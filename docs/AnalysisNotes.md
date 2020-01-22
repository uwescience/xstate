# Analysis Notes

1. Want a way to detect when gene expression has changed in time course data under weak assumptions. Can do a kind of non-parameteric change point detection.
   1. Consider the comparison of expression levels over 3 sequential time periods. If there is a consistent directional change, then expression is either increasing over all three periods or decreasing. Under the null hypothesis of not change, this happens with probability 0.25.
   1. Let $x_{tn}$ be the expression level of gene $n$ at time $t$. The probability of $K$ of $N$ genes have three consistent direcitonal changes is $P_K = C(N, K) 0.25^K 0.75^{N-K}$, where $C(N,K)$ is $N$ choose $K$. Thus, the significance level is $\sum_{i=K}^N P_i$.
   1. $P_K$ can be approximated with a normal.
 
1. State identification via change points and characterization via gene expression. One goal is to reduce the effect of gene expression unrelated to state transitions.
   1. The original data are $X$ which is $N \times P$, $T$ genes by $T$ times.
   1. Do PCA. Let $e_k$ be the $k-th$ largest eigenvector of $X^T X$; the eigenvector has dimension $N$.
   1. Select the set of genes $F_m \subset \{1, \cdots, N\}$, such that these genes have the $m$ largest values $E_{1,2}$, where $E_{1,2} = \{ e_{n,1}, e_{n,2}, 1 \leq n \leq N  \}$.
   1. $X^{\prime}$ is $X$ with only the rows in $F_m$.
   1. Let $\cal{B} =$ $\{ B_i \}$, biclusters on $X^{\prime}$. The $B_i$ have elements $(n, t)$, and the $B_i$ need not be disjoint.
   1. Infer the change points from $\cal{B}$ for each replication and hence the states. Also, use $\cal{B}$ to infer the genes that characterize each state.