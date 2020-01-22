import os

# Data Frame Columns
COUNT = "count"
CNT_DOWN = "cnt_down"  # Count of -1s present
CNT_GROUP = "cnt_group"  # Count of members in the same group
CNT_GENE = "cnt_gene"  # Number of genes
CNT_REGULATED = "cnt_regulated"  # Count of times up- down-regulated
CNT_TERM = "cnt_term"  # Number of GO terms"
CNT_UP = "cnt_up"  # Count of +1s present
CV = "cv"  # coefficient of variation
END = "End"  # Start position of gene
GENE_ID = "GENE_ID"  # Unique string for gene
GENE_IDS = "GENE_IDs"  # List of GENE_ID
GENE_NAME = "GENE_NAME"
GROUP = "group"
GO_LINK = "GO_Link"
GO_ONTOLOGY = "GO_Ontology"
GO_TERM = "GO_Term"
GROUP = "group"
HOURS = "hours"
INDEX = "index"
LENGTH = "Length"  # Length of gene in base pairs
MEAN = "mean"  # average value
PRODUCT = "PRODUCT"  # Transcript product
SAMPLE = "sample"  # Identity of a sample
TERM = "term"  # GO term
STAGE_NAME = "name"  # Name of the stage
STAGE_COLOR = "color"  # Color for the stage
START = "Start"  # Starting position of gene
STATE = "state"
STRAND = "Strand"  # either "+" or "-"
STD = "std"  # Standard deviation
TIMEPOINT = "timepoint"  # time in the experiment

# Paths

PROJECT_DIR = os.path.abspath(__file__)
for _ in range(4):
  PROJECT_DIR = os.path.dirname(PROJECT_DIR)
DATA_DIR = os.path.join(PROJECT_DIR, 'data')
CODE_DIR = PROJECT_DIR
for directory in ["xstate", "python"]:
  CODE_DIR = os.path.join(CODE_DIR, directory)
TEST_DIR = os.path.join(CODE_DIR, "tests")
DATA_PROVIDER_PERSISTER_PATH = os.path.join(
    os.path.dirname(os.path.abspath(__file__)),
    "data_provider.pcl")
SAMPLES_DIR = os.path.join(DATA_DIR, "samples")
ENSEMBLE_PATH = os.path.join(DATA_DIR, "ensemble.pcl")

# Data Characeristics
NUM_TIMES = 26
REF_TIME = 0  # Reference time
