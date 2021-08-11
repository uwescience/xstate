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
