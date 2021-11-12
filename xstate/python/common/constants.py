import os

# Special Values
TIME_0 = "T0"  # The reference column

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
FEATURE = "feature"
GENE_DESCRIPTION = "gene_description" 
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
SCORE = "score"
SIGN = "sign"  # 1 or -1
STAGE_NAME = "name"  # Name of the stage
STAGE_COLOR = "color"  # Color for the stage
START = "Start"  # Starting position of gene
STATE = "state"
STRAND = "Strand"  # either "+" or "-"
STD = "std"  # Standard deviation
TERM = "term"  # GO term
TF = "tf"  # transcription factor
TIMEPOINT = "timepoint"  # time in the experiment
FIT_RESULT_COLUMNS = [
    STATE, GROUP, GENE_ID, SCORE, COUNT
    ]

# Paths

PROJECT_DIR = os.path.abspath(__file__)
for _ in range(4):
  PROJECT_DIR = os.path.dirname(PROJECT_DIR)
DATA_DIR = os.path.join(PROJECT_DIR, 'data')
SAMPLES_DIR = os.path.join(DATA_DIR, "samples")
TRINARY_SAMPLES_DIR = os.path.join(SAMPLES_DIR,
    "trinary_data")
CODE_DIR = PROJECT_DIR
for directory in ["xstate", "python"]:
  CODE_DIR = os.path.join(CODE_DIR, directory)
TEST_DIR = os.path.join(CODE_DIR, "tests")
DATA_PROVIDER_PERSISTER_PATH = os.path.join(
    os.path.dirname(os.path.abspath(__file__)),
    "data_provider.pcl")
ENSEMBLE_PATH = os.path.join(DATA_DIR, "ensemble.pcl")

# Data Characeristics
NUM_TIMES = 26
STATE_NORMOXIA = "Normoxia"
STATE_RESCUSCITATION = "Resuscitation"
GENE_SEPARATOR = "--"

# Data transformations
MIN_VALUE = 10e-3  # Minimum value for an average count
