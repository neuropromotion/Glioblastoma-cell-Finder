"""Variational Autoencoder for fixed gene-panel cell embeddings.

Supports high-dimensional inputs (e.g. all genes on chr7+chr10, ~7k).
Train once, then project new cells into the same latent space via the
frozen encoder (+ fitted scaler / feature list).
"""

from __future__ import annotations

import json
from dataclasses import asdict, dataclass
from pathlib import Path
from typing import Sequence

import numpy as np
import pandas as pd
import torch
import torch.nn as nn
from sklearn.preprocessing import StandardScaler
from torch.utils.data import DataLoader, TensorDataset


@dataclass
class VAEConfig:
    input_dim: int
    latent_dim: int = 3
    # Wide first layer matters for ~7k genes; 128 is too narrow.
    hidden_dims: tuple[int, ...] = (1024, 256, 64)
    beta: float = 0.01
    beta_warmup_epochs: int = 40
    lr: float = 5e-4
    weight_decay: float = 1e-5
    batch_size: int = 512
    epochs: int = 150
    seed: int = 42
    dropout: float = 0.1
    # Drop near-constant genes before training (fraction of cells with expr > 0)
    min_detection_rate: float = 0.01
    # Keep full training tensor on GPU (fast for tens of thousands of cells)
    preload_to_device: bool = True
    feature_names: list[str] | None = None


class ChromosomeVAE(nn.Module):
    def __init__(
        self,
        input_dim: int,
        latent_dim: int = 3,
        hidden_dims: Sequence[int] = (1024, 256, 64),
        dropout: float = 0.1,
    ):
        super().__init__()
        self.input_dim = input_dim
        self.latent_dim = latent_dim

        enc: list[nn.Module] = []
        prev = input_dim
        for h in hidden_dims:
            enc += [
                nn.Linear(prev, h),
                nn.BatchNorm1d(h),
                nn.ReLU(),
                nn.Dropout(dropout),
            ]
            prev = h
        self.encoder = nn.Sequential(*enc)
        self.fc_mu = nn.Linear(prev, latent_dim)
        self.fc_logvar = nn.Linear(prev, latent_dim)

        dec: list[nn.Module] = []
        prev = latent_dim
        for h in reversed(hidden_dims):
            dec += [
                nn.Linear(prev, h),
                nn.BatchNorm1d(h),
                nn.ReLU(),
                nn.Dropout(dropout),
            ]
            prev = h
        dec.append(nn.Linear(prev, input_dim))
        self.decoder = nn.Sequential(*dec)

    def encode(self, x: torch.Tensor) -> tuple[torch.Tensor, torch.Tensor]:
        h = self.encoder(x)
        mu = self.fc_mu(h)
        logvar = self.fc_logvar(h).clamp(-10.0, 10.0)
        return mu, logvar

    @staticmethod
    def reparameterize(mu: torch.Tensor, logvar: torch.Tensor) -> torch.Tensor:
        std = torch.exp(0.5 * logvar)
        eps = torch.randn_like(std)
        return mu + eps * std

    def decode(self, z: torch.Tensor) -> torch.Tensor:
        return self.decoder(z)

    def forward(
        self, x: torch.Tensor
    ) -> tuple[torch.Tensor, torch.Tensor, torch.Tensor]:
        mu, logvar = self.encode(x)
        z = self.reparameterize(mu, logvar)
        recon = self.decode(z)
        return recon, mu, logvar


def vae_loss(
    recon: torch.Tensor,
    x: torch.Tensor,
    mu: torch.Tensor,
    logvar: torch.Tensor,
    beta: float = 1.0,
) -> tuple[torch.Tensor, torch.Tensor, torch.Tensor]:
    recon_loss = nn.functional.mse_loss(recon, x, reduction="mean")
    kl = -0.5 * torch.mean(1 + logvar - mu.pow(2) - logvar.exp())
    return recon_loss + beta * kl, recon_loss, kl


def _get_device(device: str | None = None) -> torch.device:
    if device is not None:
        return torch.device(device)
    if torch.backends.mps.is_available():
        return torch.device("mps")
    if torch.cuda.is_available():
        return torch.device("cuda")
    return torch.device("cpu")


def _resolve_feature_index(chr_means: pd.DataFrame, name: str):
    """Match feature name against index that may be str or int."""
    if name in chr_means.index:
        return name
    if name.isdigit():
        as_int = int(name)
        if as_int in chr_means.index:
            return as_int
    return None


def dataframe_to_matrix(
    chr_means: pd.DataFrame,
    feature_names: Sequence[str] | None = None,
) -> tuple[np.ndarray, pd.Index, list[str]]:
    """Convert preprocessor output (features × cells) to cells × features."""
    index_as_str = [str(x) for x in chr_means.index]
    if feature_names is None:
        feature_names = index_as_str

    resolved = []
    missing = []
    for f in feature_names:
        key = _resolve_feature_index(chr_means, str(f))
        if key is None:
            missing.append(f)
        else:
            resolved.append(key)
    if missing:
        raise KeyError(
            f"Features missing in matrix index ({len(missing)}): {missing[:5]}..."
        )

    X = chr_means.loc[resolved].T.to_numpy(dtype=np.float32)
    return X, chr_means.columns, [str(f) for f in feature_names]


def filter_variable_features(
    chr_means: pd.DataFrame,
    feature_names: Sequence[str] | None = None,
    min_detection_rate: float = 0.01,
) -> list[str]:
    """Keep genes detected ( > 0 ) in at least `min_detection_rate` of cells."""
    X, _, feats = dataframe_to_matrix(chr_means, feature_names)
    det = (X > 0).mean(axis=0)
    keep = [f for f, d in zip(feats, det) if d >= min_detection_rate]
    if not keep:
        raise RuntimeError(
            f"No features passed min_detection_rate={min_detection_rate}. "
            "Lower the threshold or check preprocessing."
        )
    return keep


class ChromosomeVAETrainer:
    """Fit scaler + VAE; encode new cells into a fixed latent space."""

    def __init__(self, config: VAEConfig | None = None, device: str | None = None):
        self.config = config or VAEConfig(input_dim=0)
        self.device = _get_device(device)
        self.scaler = StandardScaler()
        self.model: ChromosomeVAE | None = None
        self.history: list[dict[str, float]] = []
        self.feature_names: list[str] | None = self.config.feature_names
        self.cell_index_: pd.Index | None = None

    def fit(
        self,
        chr_means: pd.DataFrame,
        feature_names: Sequence[str] | None = None,
        epochs: int | None = None,
        verbose: bool = True,
    ) -> ChromosomeVAETrainer:
        raw_features = list(feature_names) if feature_names is not None else None
        if self.config.min_detection_rate > 0:
            feats = filter_variable_features(
                chr_means,
                raw_features,
                min_detection_rate=self.config.min_detection_rate,
            )
            if verbose:
                n_all = (
                    len(raw_features)
                    if raw_features is not None
                    else chr_means.shape[0]
                )
                print(
                    f"Variable features: {len(feats)}/{n_all} "
                    f"(detection ≥ {self.config.min_detection_rate:.1%})"
                )
        else:
            _, _, feats = dataframe_to_matrix(chr_means, raw_features)

        X, cell_index, feats = dataframe_to_matrix(chr_means, feats)
        self.feature_names = feats
        self.cell_index_ = cell_index
        self.config.input_dim = X.shape[1]
        self.config.feature_names = feats

        if verbose:
            print(
                f"Train matrix: {X.shape[0]} cells × {X.shape[1]} genes | "
                f"device={self.device}"
            )

        torch.manual_seed(self.config.seed)
        np.random.seed(self.config.seed)

        Xs = self.scaler.fit_transform(X).astype(np.float32)
        # Guard against rare NaNs from zero-variance columns
        Xs = np.nan_to_num(Xs, nan=0.0, posinf=0.0, neginf=0.0)

        self.model = ChromosomeVAE(
            input_dim=self.config.input_dim,
            latent_dim=self.config.latent_dim,
            hidden_dims=self.config.hidden_dims,
            dropout=self.config.dropout,
        ).to(self.device)
        opt = torch.optim.Adam(
            self.model.parameters(),
            lr=self.config.lr,
            weight_decay=self.config.weight_decay,
        )

        n_cells = Xs.shape[0]
        batch_size = min(self.config.batch_size, max(n_cells, 1))
        use_gpu_preload = (
            self.config.preload_to_device and self.device.type in {"cuda", "mps"}
        )
        if use_gpu_preload:
            # Avoid per-batch CPU→GPU copies (main bottleneck on MPS)
            Xt = torch.from_numpy(Xs).to(self.device)
            if verbose:
                mb = Xt.element_size() * Xt.nelement() / 1e6
                print(f"Preloaded training tensor on {self.device} ({mb:.1f} MB)")
        else:
            loader = DataLoader(
                TensorDataset(torch.from_numpy(Xs)),
                batch_size=batch_size,
                shuffle=True,
                drop_last=n_cells > batch_size,
            )

        n_epochs = epochs if epochs is not None else self.config.epochs
        warmup = max(int(self.config.beta_warmup_epochs), 0)
        self.history = []

        for epoch in range(1, n_epochs + 1):
            self.model.train()
            if warmup <= 0:
                beta = self.config.beta
            else:
                beta = self.config.beta * min(1.0, epoch / warmup)

            total = recon_sum = kl_sum = 0.0
            n_batches = 0

            if use_gpu_preload:
                perm = torch.randperm(n_cells, device=self.device)
                if n_cells > batch_size:
                    n_use = (n_cells // batch_size) * batch_size
                    perm = perm[:n_use]
                for start in range(0, len(perm), batch_size):
                    xb = Xt.index_select(0, perm[start : start + batch_size])
                    opt.zero_grad(set_to_none=True)
                    recon, mu, logvar = self.model(xb)
                    loss, recon_l, kl_l = vae_loss(recon, xb, mu, logvar, beta=beta)
                    loss.backward()
                    nn.utils.clip_grad_norm_(self.model.parameters(), 5.0)
                    opt.step()
                    total += float(loss.detach())
                    recon_sum += float(recon_l.detach())
                    kl_sum += float(kl_l.detach())
                    n_batches += 1
            else:
                for (xb,) in loader:
                    xb = xb.to(self.device)
                    opt.zero_grad(set_to_none=True)
                    recon, mu, logvar = self.model(xb)
                    loss, recon_l, kl_l = vae_loss(recon, xb, mu, logvar, beta=beta)
                    loss.backward()
                    nn.utils.clip_grad_norm_(self.model.parameters(), 5.0)
                    opt.step()
                    total += float(loss.item())
                    recon_sum += float(recon_l.item())
                    kl_sum += float(kl_l.item())
                    n_batches += 1

            row = {
                "epoch": epoch,
                "beta": beta,
                "loss": total / max(n_batches, 1),
                "recon": recon_sum / max(n_batches, 1),
                "kl": kl_sum / max(n_batches, 1),
            }
            self.history.append(row)
            if verbose and (epoch == 1 or epoch % 10 == 0 or epoch == n_epochs):
                print(
                    f"epoch {epoch:4d}/{n_epochs}  beta={beta:.4f}  "
                    f"loss={row['loss']:.4f}  recon={row['recon']:.4f}  kl={row['kl']:.4f}"
                )
        return self

    @torch.no_grad()
    def encode(
        self,
        chr_means: pd.DataFrame,
        *,
        deterministic: bool = True,
    ) -> pd.DataFrame:
        """Project cells into latent space. Uses μ when deterministic=True."""
        if self.model is None or self.feature_names is None:
            raise RuntimeError("Model is not fitted. Call fit() or load() first.")

        X, cell_index, _ = dataframe_to_matrix(chr_means, self.feature_names)
        Xs = self.scaler.transform(X).astype(np.float32)
        Xs = np.nan_to_num(Xs, nan=0.0, posinf=0.0, neginf=0.0)
        xt = torch.from_numpy(Xs).to(self.device)

        self.model.eval()
        mu, logvar = self.model.encode(xt)
        z = mu if deterministic else self.model.reparameterize(mu, logvar)
        z_np = z.cpu().numpy()

        cols = [f"z{i + 1}" for i in range(self.config.latent_dim)]
        return pd.DataFrame(z_np, index=cell_index, columns=cols)

    def fit_encode(
        self,
        chr_means: pd.DataFrame,
        feature_names: Sequence[str] | None = None,
        **fit_kwargs,
    ) -> pd.DataFrame:
        self.fit(chr_means, feature_names=feature_names, **fit_kwargs)
        return self.encode(chr_means)

    def save(self, path: str | Path) -> None:
        if self.model is None or self.feature_names is None:
            raise RuntimeError("Nothing to save: model is not fitted.")
        path = Path(path)
        path.mkdir(parents=True, exist_ok=True)

        torch.save(self.model.state_dict(), path / "model.pt")
        with open(path / "scaler.json", "w", encoding="utf-8") as f:
            json.dump(
                {
                    "mean": self.scaler.mean_.tolist(),
                    "scale": self.scaler.scale_.tolist(),
                    "n_features_in": int(self.scaler.n_features_in_),
                },
                f,
            )
        with open(path / "config.json", "w", encoding="utf-8") as f:
            cfg = asdict(self.config)
            cfg["hidden_dims"] = list(self.config.hidden_dims)
            json.dump(cfg, f, indent=2)
        with open(path / "history.json", "w", encoding="utf-8") as f:
            json.dump(self.history, f)

    @classmethod
    def load(cls, path: str | Path, device: str | None = None) -> ChromosomeVAETrainer:
        path = Path(path)
        with open(path / "config.json", encoding="utf-8") as f:
            raw = json.load(f)
        raw["hidden_dims"] = tuple(raw["hidden_dims"])
        # backward-compatible defaults for older checkpoints
        raw.setdefault("weight_decay", 1e-5)
        raw.setdefault("dropout", 0.1)
        raw.setdefault("min_detection_rate", 0.0)
        raw.setdefault("preload_to_device", True)
        config = VAEConfig(**raw)

        trainer = cls(config=config, device=device)
        trainer.feature_names = config.feature_names

        with open(path / "scaler.json", encoding="utf-8") as f:
            sc = json.load(f)
        trainer.scaler = StandardScaler()
        trainer.scaler.mean_ = np.asarray(sc["mean"], dtype=np.float64)
        trainer.scaler.scale_ = np.asarray(sc["scale"], dtype=np.float64)
        trainer.scaler.var_ = trainer.scaler.scale_**2
        trainer.scaler.n_features_in_ = sc["n_features_in"]
        trainer.scaler.n_samples_seen_ = 1

        trainer.model = ChromosomeVAE(
            input_dim=config.input_dim,
            latent_dim=config.latent_dim,
            hidden_dims=config.hidden_dims,
            dropout=config.dropout,
        ).to(trainer.device)
        state = torch.load(path / "model.pt", map_location=trainer.device, weights_only=True)
        trainer.model.load_state_dict(state)
        trainer.model.eval()
        return trainer
