# Analyzing Transcriptome State Dynamics

## Problem Description

1. By *gene expression state (state)* is meant that for a time period or condition the set of gene expressions within the state are readily distinguished from gene expressions outside the state.
2. There is existing work on characterizing the differences in gene expression where different gene expressions are in distinctive states.
3. However, there is an additional challenge with dynamics in that the state transitions are not provided a priori. Thus, states must first be identified. Further, state characterization in a way that makes sense for a sequence of states.
4. To summarize, we see the following challenges in the analysis of transcriptome state dynamics:
    1. **Identification**. Identify the time boundaries of a state.
    2. **Characterization**. Characterize each state in a parsimonious and distinctive way in terms of the genes that are differentially expressed. The genes constituting the state are referred to as the **state vector**.
    3. **Assignment**. Given a new state vector, determine the probability that the state vector belongs to each of a set of previously characterized states.

## State Identification For Time Course Expression Data

### Statistical Test
This approach to state identification performs a test at each time instant to detect if a state transition has occurred. Two assumptions underlie this approach:

1. It should be rare that a differential expression is sustained across a state boundary.
2. The null distribution for differential expression is that gene expression is independent and identically distributed at each time period for each gene with fixed probability $p$.
3. Genes are expressed independently of each other.

Assumption (3) is violated if genes are in the same operon or controlled by the same TF. Thus, there is a need to group genes under the same controls. This can be done with a priori knowledge (gene networks) or statistically by finding correlations.

Define $G$ as the set of grouped, differentially expressed genes. Define $q$ to be the fraction of time that a gene in $G$ is differentially expressed. Let $t$ be a time at which we are evaluating if there is a state transition. That is, we want to determine if there is a significantly low number of genes that are differentially expressed at time $t-1$ and $t$.

Suppose we want to determine if there is a state transition at time $t$.
That is, are there significantly fewer genes that are differentially expressed at both $t-1$ and $t$.
For a single gene, this happens with probability $q^2$ under the null hypothesis.
For $M$ gene groups, the probability of $m$ genes being differentially expressed at $t-1$ and $t$ is
$B(M, m, q^2)$, where this is the binomial probabililty.
Thus, the probability that there are too few contiguous differential expressions is $\alpha_t = \sum_{k=0}^m B(M, k, q^2)$. 

There is an additional consideration since we will repeatedly apply the above test to all $T$. We want to choose a cut-off value of $\alpha^\star$ such that the probability of having $n$ times with $\alpha_t < \alpha^\star$ is less than $p$, a desired significance level. This is a multiple hypothesis testing problem.

We renumber the time points as $n_1 , \cdots, n_T$ such that $\alpha_1 \leq \cdots \leq \alpha_T$ computed as above for the $T$ time periods, where $\alpha_1 = 0$.
We choose the smallest $k$ such that
$B(T, k+1, \alpha_{k+1}) > p$, or $T$ if no such $k$ exists. Then, state transitions occur at times $n_1 , \cdots, n_k$.


### Contig Analysis
A differential expression contig (**DEC**) is specified by a gene and a time period over which that gene is differentially expressed in the same direction. DECs are used to evaluate the quality of a state identification algorithm.

Consider a set of genes $G$ and times over which these genes may be differentially expressed $T = \{t_1 , \cdots, t_N \}$. A DEC $c$ specifies a $g \in G$ and a range of times from $t_i$ to $t_j$ ($t_i \leq t_j$) such that $g$ is differentially expressed in the same direction for $[t_i, t_j]$.

Let $S$ be a state description in terms of times. That is, $S = \{t_{i_1}, \cdots, t_{i_n} \}$ where the times indicate the state boundaries.

A quality measure is the area under the curve where the x-axis is fraction of DE contig in a stage and the y-axis is the fraction of the stage covered by the DE-contig.

1. Ideally want genes that are DE in just one stage.
2. Terms
   1. State characterization. Time based id of states.
   2. State identification. Gene based id of states
   3. State assignment. Given new data, determine its state.
1. Measure of quality for state characterization.
   1. Area of DEC within each state divided by the total area of the contig.
   2. Area of DEC with the states divided by the total area of the state for the associated genes within the state.