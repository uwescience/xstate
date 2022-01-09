# Data Samples to Classify

## Description of data

- AM - alveolar macrophages (total lung) taken D20 after MTB infection
- MDM - monocyte derived macrophages (total lung) taken D20 after MTB infection
- AW pos - alveolar macrophages isolated from airway, D18 after MTB infection
- AW neg - alveolar macrophages isolated from interstitium (deeper lung tissue), D18 after MTB infection.
- galagan - time series hypoxia studies from Galagan lab
- Rustad - "a link between hypoxic microenvironments within the infected host and the latent phase of tuberculosis. Studies to test this correlation have identified the M. tuberculosis initial hypoxic response, controlled by the two-component response regulator DosR". Considers the effect of DosR (Rv2004c) deletion.

## AM and MDM data
- Eliza, 7/2/2020
  - Using differentially expressed genes and looking at whether these genes were also differentially expressed in AM vs MDM, we found a significant late hypoxia signature in MDMs. 
The AMs had a significant overlap with early/mid hypoxia differentially expressed genes. So, my expectation was that there would be a difference in terms of hypoxia state.
 - Again, using differentially expressed genes we found the AW_plus had a late hypoxia signature, but not the AW_neg.

## Rustad
- Reference: [Rustad et al.]( https://journals.plos.org/plosone/article?id=10.1371/journal.pone.0001502#s4)
- Characteristics
   * 4 replications
   * Data is already in normalized and in units of log2 w.r.t. reference
   * Dual array so that there is a reference for each sample

## GSE167232\_mtb\_transcriptome\_counts\_normalized.csv
From Eliza on 10/24/2021.
The data is from Pisu et al 2021 (attached), it is MTB from lung macrophages of infected mice, sorted on the MTB expression of a stress reporter gene (hspX). These are the TB_HIGH and TB_LOW samples. The other samples are from Pisu et al 2020, it is MTB from either alveolar macrophages (AM) or interstitial macrophages (IM) of infected mice.
The data are normalized for library size and gene size but are not
in log2 units.

## GSE167232\_mtb\_transcriptome\_counts\_normalized_reduced.csv
* Created on 1/7/2022.
* Deleted columns other than TB\_HIGH* and TB\_LOW*


``GSE167232\_mtb\_transcriptome\_counts\_normalized\_filtered.csv``
is the same data but removes all lines with "N/A" and deletes the
"Gene\_name" column.

Notes
1. Has duplicate values for gene Rv0687
