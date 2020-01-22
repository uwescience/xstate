# Analysis Plan

## Background

The analysis will focus on developing and validating algorithms for state characterization and assignment.

At present, we only have one example of ground truth for state assignments, with six states. We will extend this by proposing error models for read counts that create synthetic data to evaluate the robustness of state characterization and assignment.

## Technical Description of Tasks

A **state vector** is a vector of gene expression values. A **generalized state vector (GSV)** includes information beyond gene expression values such as: include: genes for which no value is given, whether the gene is up- or down-regulated, and variance of values. One representation of a GSV is a vector of means and a vector of variances. A variance of ``inf`` means that a gene is ignored.

The **state characterization task** specifies one or more extended state vectors for each state. The characterization of a state should adequately distinguish that state from other states.


The **state assignment task** assigns a state vector to one or more states based on the state characterizations. Assignments may be probabilistic (or weighted in some way) to indicate uncertainty as to the state assignment.



## Factors in Studies

1. Algorithm

1. Model for perturbation of read counts (error model).
   1. Random reordering of read counts.
   2. Gaussian error model with fixed variance (or coefficient of variation).

1. Units for gene expression.
   
   Units should take into account that genes may be up- or down-regulated.
   1. Normalized read counts
   2. Trinary values of read counts indicating either up- or down-regulation (w.r.t. library normalized read counts for reference state).
   2. z-score
   3. significance level
   4. Flux, change in values over times.
   5. Adjusted raw read counts (since this is what's available for a state vector)

1. Distance measure between state vectors (over time and multiple state vectors)

   The distance measure should be defined between state vectors, extended state vectors, and extended and non-extended state vectors.
   1. Eucleadian distance
   2. for significance level, infinite if below threshold and 0 if above threshold
   3. Handling distance between 1 and n state vectors (e.g., characterization). Handling distance to a set (e.g., centroid, maximum distance). Adjusting where there is variability at one point in time (e.g., weight by inverse of standard deviation plus a term $\lambda$ to handle 0 variance)?


## Initial Analysis

1. Develop a framework for analyzing the quality of state characterization.
1. Show that high quality state characterization can be accomplished using read counts.

## Notes

1. Is it possible to do state characterization without state assignment since the validation of a characterization is its effectiveness at assignment?