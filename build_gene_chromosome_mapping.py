#!/usr/bin/env python3
"""Download the latest Ensembl human GTF and build gene → chromosome mapping.

Output TSV has three columns:
  1) gene        — HGNC symbol and Ensembl gene ID (separate rows)
  2) chromosome
  3) arm         — p / q / cen (empty for MT)
"""

from __future__ import annotations

import argparse
import gzip
import re
import ssl
import subprocess
import sys
import urllib.error
import urllib.request
from pathlib import Path

ENSEMBL_PUB = "https://ftp.ensembl.org/pub"
UCSC_CYTOBAND = {
    "GRCh38": "https://hgdownload.soe.ucsc.edu/goldenPath/hg38/database/cytoBand.txt.gz",
    "GRCh37": "https://hgdownload.soe.ucsc.edu/goldenPath/hg19/database/cytoBand.txt.gz",
}
GTF_NAME_RE = re.compile(
    r'href="(Homo_sapiens\.(GRCh\d+)\.(\d+)\.gtf\.gz)"',
    re.IGNORECASE,
)
ATTR_RE = re.compile(r'(\w+)\s+"([^"]*)"')


def _ssl_context() -> ssl.SSLContext:
    """Prefer system certs; fall back if corporate/self-signed proxies break verify."""
    try:
        return ssl.create_default_context()
    except Exception:  # noqa: BLE001
        return ssl._create_unverified_context()


def fetch_bytes(url: str, timeout: int = 120) -> bytes:
    """Fetch URL bytes via urllib, then curl fallback (handles broken SSL stores)."""
    req = urllib.request.Request(
        url,
        headers={"User-Agent": "gene-chromosome-mapping/1.0"},
    )
    try:
        with urllib.request.urlopen(req, timeout=timeout, context=_ssl_context()) as resp:
            return resp.read()
    except Exception:
        pass

    try:
        ctx = ssl._create_unverified_context()
        with urllib.request.urlopen(req, timeout=timeout, context=ctx) as resp:
            return resp.read()
    except Exception:
        pass

    result = subprocess.run(
        ["curl", "-fsSL", "--retry", "3", url],
        check=False,
        capture_output=True,
    )
    if result.returncode != 0:
        err = result.stderr.decode("utf-8", errors="replace").strip()
        raise RuntimeError(f"Failed to download {url}: {err or result.returncode}")
    return result.stdout


def fetch_text(url: str, timeout: int = 120) -> str:
    return fetch_bytes(url, timeout=timeout).decode("utf-8", errors="replace")


def download_file(url: str, dest: Path, timeout: int = 600) -> None:
    dest.parent.mkdir(parents=True, exist_ok=True)
    tmp = dest.with_suffix(dest.suffix + ".partial")
    print(f"Downloading:\n  {url}\n  -> {dest}")

    result = subprocess.run(
        [
            "curl",
            "-fL",
            "--retry",
            "3",
            "--connect-timeout",
            "30",
            "-o",
            str(tmp),
            url,
        ],
        check=False,
        capture_output=True,
    )
    if result.returncode != 0:
        data = fetch_bytes(url, timeout=timeout)
        tmp.write_bytes(data)
    tmp.replace(dest)


def resolve_latest_gtf_url() -> tuple[str, str, str]:
    """Return (filename, full_url, assembly) for the latest primary-assembly human GTF."""
    version = fetch_text(f"{ENSEMBL_PUB}/VERSION").strip()
    if not version.isdigit():
        raise RuntimeError(f"Unexpected Ensembl VERSION content: {version!r}")

    candidates = [
        f"{ENSEMBL_PUB}/current/gtf/homo_sapiens/",
        f"{ENSEMBL_PUB}/release-{version}/gtf/homo_sapiens/",
    ]

    last_err: Exception | None = None
    for index_url in candidates:
        try:
            html = fetch_text(index_url)
        except Exception as exc:  # noqa: BLE001
            last_err = exc
            continue

        matches = GTF_NAME_RE.findall(html)
        primary = [
            (name, assembly)
            for name, assembly, _rel in matches
            if ".chr." not in name
            and "abinitio" not in name
            and "hapl" not in name
            and "scaffold" not in name
        ]
        # unique by name, keep last (highest release in sorted order)
        by_name = {name: assembly for name, assembly in sorted(primary)}
        if by_name:
            filename = sorted(by_name)[-1]
            assembly = by_name[filename]
            print(f"Ensembl release {version}: {filename}")
            return filename, index_url.rstrip("/") + "/" + filename, assembly

    raise RuntimeError(
        f"Could not find Homo_sapiens GTF under {candidates}. Last error: {last_err}"
    )


def assembly_from_gtf_name(name: str) -> str:
    m = re.search(r"\.(GRCh\d+)\.", name)
    if not m:
        raise RuntimeError(f"Cannot infer assembly from GTF name: {name}")
    return m.group(1)


def parse_gtf_attributes(attr_field: str) -> dict[str, str]:
    return dict(ATTR_RE.findall(attr_field))


def normalize_chromosome(chrom: str, strip_chr_prefix: bool) -> str:
    if strip_chr_prefix and chrom.lower().startswith("chr"):
        chrom = chrom[3:]
    if chrom.upper() in {"M", "CHRM", "CHRMT"}:
        return "MT"
    if chrom.upper() == "MT":
        return "MT"
    return chrom


def _is_primary_chrom(chrom: str) -> bool:
    if chrom in {"X", "Y", "MT"}:
        return True
    return chrom.isdigit() and 1 <= int(chrom) <= 22


def load_centromeres(cytoband_path: Path) -> dict[str, tuple[int, int]]:
    """Parse UCSC cytoBand; return chrom -> (cen_start, cen_end) for acen bands."""
    open_fn = gzip.open if cytoband_path.name.endswith(".gz") else open
    spans: dict[str, list[tuple[int, int]]] = {}

    with open_fn(cytoband_path, "rt", encoding="utf-8") as fh:
        for line in fh:
            parts = line.rstrip("\n").split("\t")
            if len(parts) < 5 or parts[4] != "acen":
                continue
            chrom = normalize_chromosome(parts[0], strip_chr_prefix=True)
            start, end = int(parts[1]), int(parts[2])
            spans.setdefault(chrom, []).append((start, end))

    centromeres: dict[str, tuple[int, int]] = {}
    for chrom, intervals in spans.items():
        cen_start = min(s for s, _ in intervals)
        cen_end = max(e for _, e in intervals)
        centromeres[chrom] = (cen_start, cen_end)
    return centromeres


def assign_arm(
    chrom: str,
    start: int,
    end: int,
    centromeres: dict[str, tuple[int, int]],
) -> str:
    """Assign p / q / cen from gene midpoint vs centromere span. MT → empty."""
    if chrom == "MT":
        return ""
    cen = centromeres.get(chrom)
    if cen is None:
        return ""
    cen_start, cen_end = cen
    mid = (start + end) // 2
    if mid < cen_start:
        return "p"
    if mid > cen_end:
        return "q"
    return "cen"


def ensure_cytoband(
    assembly: str,
    cache_dir: Path,
    *,
    force: bool = False,
) -> Path:
    if assembly not in UCSC_CYTOBAND:
        raise RuntimeError(
            f"No UCSC cytoBand URL for assembly {assembly}. "
            f"Known: {', '.join(UCSC_CYTOBAND)}"
        )
    url = UCSC_CYTOBAND[assembly]
    dest = cache_dir / f"cytoBand_{assembly}.txt.gz"
    if force or not dest.exists():
        download_file(url, dest)
    else:
        print(f"Using cached cytoBand: {dest}")
    return dest


def build_mapping_from_gtf(
    gtf_path: Path,
    centromeres: dict[str, tuple[int, int]],
    *,
    include_scaffolds: bool = False,
    strip_chr_prefix: bool = True,
) -> list[tuple[str, str, str]]:
    """Parse gene features; return rows of (gene_label, chromosome, arm)."""
    open_fn = gzip.open if gtf_path.name.endswith(".gz") else open
    rows: list[tuple[str, str, str]] = []
    seen: set[tuple[str, str]] = set()

    with open_fn(gtf_path, "rt", encoding="utf-8") as fh:
        for line in fh:
            if not line or line.startswith("#"):
                continue
            parts = line.rstrip("\n").split("\t")
            if len(parts) < 9 or parts[2] != "gene":
                continue

            chrom = normalize_chromosome(parts[0], strip_chr_prefix)
            if not include_scaffolds and not _is_primary_chrom(chrom):
                continue

            start, end = int(parts[3]), int(parts[4])
            arm = assign_arm(chrom, start, end, centromeres)

            attrs = parse_gtf_attributes(parts[8])
            gene_id = attrs.get("gene_id", "").split(".")[0]
            gene_name = attrs.get("gene_name") or attrs.get("gene_symbol") or ""

            for label in (gene_id, gene_name):
                if not label:
                    continue
                key = (label, chrom)
                if key in seen:
                    continue
                seen.add(key)
                rows.append((label, chrom, arm))

    return rows


def write_mapping(rows: list[tuple[str, str, str]], out_path: Path) -> None:
    out_path.parent.mkdir(parents=True, exist_ok=True)
    with open(out_path, "w", encoding="utf-8", newline="\n") as out:
        out.write("gene\tchromosome\tarm\n")
        for gene, chrom, arm in rows:
            out.write(f"{gene}\t{chrom}\t{arm}\n")


def parse_args(argv: list[str] | None = None) -> argparse.Namespace:
    p = argparse.ArgumentParser(
        description=(
            "Download latest Ensembl Homo sapiens GTF (+ UCSC cytoBand) and write "
            "gene→chromosome→arm mapping (symbols + ENSG IDs)."
        )
    )
    p.add_argument(
        "-o",
        "--output",
        type=Path,
        default=Path("gene_chromosome_mapping.tsv"),
        help="Output TSV path (default: gene_chromosome_mapping.tsv)",
    )
    p.add_argument(
        "--cache-dir",
        type=Path,
        default=Path("annotation_cache"),
        help="Directory for downloaded files (default: annotation_cache)",
    )
    p.add_argument(
        "--gtf",
        type=Path,
        default=None,
        help="Use an existing GTF(.gz) instead of downloading",
    )
    p.add_argument(
        "--assembly",
        choices=sorted(UCSC_CYTOBAND),
        default=None,
        help="Assembly for cytoBand (default: inferred from GTF name)",
    )
    p.add_argument(
        "--include-scaffolds",
        action="store_true",
        help="Keep non-primary chromosomes/scaffolds (default: 1–22, X, Y, MT)",
    )
    p.add_argument(
        "--keep-chr-prefix",
        action="store_true",
        help="Do not strip leading 'chr' from chromosome names",
    )
    p.add_argument(
        "--force-download",
        action="store_true",
        help="Re-download remote files even if cached",
    )
    return p.parse_args(argv)


def main(argv: list[str] | None = None) -> int:
    args = parse_args(argv)

    try:
        if args.gtf is not None:
            gtf_path = args.gtf
            if not gtf_path.exists():
                print(f"GTF not found: {gtf_path}", file=sys.stderr)
                return 1
            print(f"Using local GTF: {gtf_path}")
            assembly = args.assembly or assembly_from_gtf_name(gtf_path.name)
        else:
            filename, url, assembly = resolve_latest_gtf_url()
            if args.assembly:
                assembly = args.assembly
            gtf_path = args.cache_dir / filename
            if args.force_download or not gtf_path.exists():
                download_file(url, gtf_path)
            else:
                print(f"Using cached GTF: {gtf_path}")

        cytoband_path = ensure_cytoband(
            assembly, args.cache_dir, force=args.force_download
        )
        centromeres = load_centromeres(cytoband_path)
        print(f"Loaded centromeres for {len(centromeres)} chromosomes ({assembly})")

        print("Parsing GTF gene features...")
        rows = build_mapping_from_gtf(
            gtf_path,
            centromeres,
            include_scaffolds=args.include_scaffolds,
            strip_chr_prefix=not args.keep_chr_prefix,
        )
        if not rows:
            print("No gene rows extracted from GTF.", file=sys.stderr)
            return 1

        write_mapping(rows, args.output)
        n_ensg = sum(1 for g, _, _ in rows if g.startswith("ENSG"))
        n_symbol = len(rows) - n_ensg
        arm_counts: dict[str, int] = {}
        for _, _, arm in rows:
            key = arm if arm else "NA"
            arm_counts[key] = arm_counts.get(key, 0) + 1
        chroms = sorted(
            {c for _, c, _ in rows},
            key=lambda x: (not x.isdigit(), int(x) if x.isdigit() else 99, x),
        )
        print(
            f"Wrote {len(rows)} rows (~{n_ensg} ENSG + ~{n_symbol} symbols) -> {args.output}"
        )
        print(f"Chromosomes: {', '.join(chroms)}")
        print(
            "Arms: "
            + ", ".join(f"{k}={v}" for k, v in sorted(arm_counts.items()))
        )
        return 0
    except urllib.error.URLError as exc:
        print(f"Network error: {exc}", file=sys.stderr)
        return 1
    except Exception as exc:  # noqa: BLE001 — CLI entrypoint
        print(f"Error: {exc}", file=sys.stderr)
        return 1


if __name__ == "__main__":
    raise SystemExit(main())
