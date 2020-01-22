# Expression Data Normalization
Given $x_{nm}$, the read count for gene $n$ in replica $m$. What normalizations should be applied to the $x{nm}$ to account for differences in the replicas? Much of this is taken from [this reference](https://hbctraining.github.io/DGE_workshop/lessons/02_DGE_count_normalization.html).

## Concerns
   1. Read depth (library size). The number of reads may vary across replicas.
   1. Gene length. More reads will be allocated to a gene if it is longer.
   1. Very large differences in expressions. This will tend to compress values.

## Resolutions
1. Read depth. Use $x^\prime_{nm} = \frac{x_{nm}}{C_m}$, where $C_m = \sum_{n=1}^N x _{nm}$.
2. Gene length. Use TPM, *transcripts per kilobase*. Divide $x_nm$ by $L_n$, the length of gene $n$.
3. Differences in expression. Use a geometric mean.