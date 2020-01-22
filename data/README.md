# Data Files

## Data from Eliza - 2/6/2019

- "hypoxia\_timecourse\_reads.csv" - non-normalized raw read counts. The reason they are not integer values lies in how reads that straddle the gene boundary edges are treated. The DuffyNGS pipeline prorates those reads based on what fraction of the called bases actually land inside the gene. So if only 40% of the read's bases land insie the gene, then the gene receives 0.40 of a read (and the intergenic space outside the gene gets the other 0.60 of the read). For passing the data to a tool that wants raw read counts (DESeq2, etc.), we just round those READS values to integers.

- "hypoxia\_curve\_DO.csv" - dissolved oxygen (%) for the timepoints along the timecourse (T0, T1, etc)

- "stages\_matrix.csv" - the assigned "state" or "stage" for the timepoints along the timecourse (T0, T1, etc)

- "normalized\_log2\_transformed\_counts.csv" - collapsed replicates, normalized and log2 transformed. I believe they were normalized as per DESeq2...need to check with Mario on this. I believe this was the data used for PCA/t-SNE plots to identify the different states.
   - Added a header for the first column - GENE\_ID

- .doc manuscript - not completed. The sections at the end that are grayed out are probably going to be cut...I would only read what is not grayed out. I realized the methods don't include the differential expression analysis yet nor do we have the figures finalized.

## Added gene\_expression\_state
For each gene\_id, specifies its associated state. Data for the gene
can be found in other files.
This file was synthesized from files of the DE genes provided
by Mario on 3/27/2019.

## MTB.GO.All.GOterms.csv
Transformed data from Eliza into a CSV with terms

## Access KEGG pathways
- list all pathways for MTB: ``http://rest.kegg.jp/list/pathway/mtv``
- get details for a pathway: ``http://rest.kegg.jp/get/path:mtv00010``

## Categorization of genes
- mtb\_kegg\_pathway.csv - KEGG pathways and descriptions
- mtb\_kegg\_gene\_pathways.csv - gene, KEGG pathway
- mtb\_gene\_ec.csv - gene, EC number
- mtb\_gene\_ko.csv - gene, KEGG orthology
