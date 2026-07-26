import pathlib
PROJ_ROOT = pathlib.Path.cwd()
HEALTHY_PATH = PROJ_ROOT / 'data' / 'healthy'
PATH_3v3_10x = PROJ_ROOT / 'data' / '3v3_10x'
PATH_GSE173278 = PROJ_ROOT / 'data' / 'GSE173278'
PATH_GSE135045 = PROJ_ROOT / 'data' / 'GSE135045'
GENE_LENGTHS_PATH = PROJ_ROOT / 'data' / 'gene_lengths.csv'
ENS_MAP_PATH = PROJ_ROOT / 'data' / 'ensembl_gene_mapping.csv'
CHR_MAPPING_PATH = PROJ_ROOT / 'data' / 'gene_chromosome_mapping.tsv'

MIN_GENES = 300
MT_THRESHOLD = 20.0
GENE_UPPER_QUANTILE = 0.95
COUNT_UPPER_QUANTILE = 0.95
MIN_CELLS = 3
