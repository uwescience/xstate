"""Creates categorization data sets."""

"""
Category data sets:
  kegg_pathways: FILE_KEGG_PATHWAY FILE_KEGG_DESCRIPTION
  kegg_gene_pathway: FILE_KEGG_GENE, FILE_KEGG_PATHWAY
  kegg_gene_ec: FILE_KEGG_GENE, FILE_KEGG_EC
  kegg_gene_ko: FILE_KEGG_GENE, FILE_KEGG_KO
"""


import common.constants as cn
import common_python.constants as cpn
from common_python.bioinformatics import kegg_extractor

import os

# Files
FILE_KEGG_PATHWAYS = "mtb_kegg_pathways.csv"
FILE_KEGG_GENE_PATHWAY = "mtb_kegg_gene_pathway.csv"
FILE_KEGG_GENE_EC = "mtb_gene_ec.csv"
FILE_KEGG_GENE_KO = "mtb_gene_ko.csv"
FILES_KEGG = [
    FILE_KEGG_PATHWAYS,
    FILE_KEGG_GENE_PATHWAY,
    FILE_KEGG_GENE_EC,
    FILE_KEGG_GENE_KO,
    ]
    
# Names
MTV = "mtv"

def makePath(filename, directory):
  return os.path.join(directory, filename)

def writeFile(filename, df_base, columns, directory):
  """
  Writes the columns of the base DF to a file.
  :param str filename: filename to write
  :param pd.DataFrame df_base: DF from which to extract columns
  :param list-str columns:
  :param str directory: directory where file will be written
  """
  df = df_base[columns]
  df = df.copy()
  if cpn.KEGG_GENE in columns:
    df[cn.GENE_ID] = df[cpn.KEGG_GENE]
    df = df.drop(columns=cpn.KEGG_GENE)
  df = df.drop_duplicates()
  df.to_csv(makePath(filename, directory), index=False)

def make(max_count=-1, directory=cn.DATA_DIR):
  """
  Makes the categorization files.
  :param int max_count: maximum number of records to process.
      <0 means to process all.
  :param str directory: directory where files reside
  """
  extractor = kegg_extractor.KeggExtractor(MTV)
  # Do pathways
  df_pathways = extractor.listPathway()
  writeFile(FILE_KEGG_PATHWAYS, df_pathways, 
      df_pathways.columns, directory)
  # Other files
  df_base = extractor.getAllPathwayGenes(max_count=max_count)
  writeFile(FILE_KEGG_GENE_PATHWAY, df_base, 
      [cpn.KEGG_GENE, cpn.KEGG_PATHWAY], directory)
  writeFile(FILE_KEGG_GENE_EC, df_base, 
      [cpn.KEGG_GENE, cpn.KEGG_EC], directory)
  writeFile(FILE_KEGG_GENE_KO, df_base, 
      [cpn.KEGG_GENE, cpn.KEGG_KO], directory)


if __name__ == '__main__':
  make()
