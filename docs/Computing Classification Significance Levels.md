# Computing Significance Levels For Case Based Classification
This is about classifying gene expression data. The data consist of columns of genes (features) and rows of instances. Values are in the set $\{-1, 0, 1\}$ that indicate under expression, no differential expression, and over expression. The classification problem is to find the state that best matches a new feature vector.

Case based classification involves counting the (positive) cases for a state $s$ in a feature vector. Typically, an experiment is replicated and so the matching involves multiple feature vectors. The result is a count of $n_{sc}$ occurrences of case $c$ in state $s$. The statistical questions are:

1. If we consider a single state and case, is $n_{sc}$ statistical significant?
2. If we consider a single state, are the cases for $s$ statistically significant?
3. Considering all states, are the results statistically significant?

## Notation
### Data
1. $N$ is the number of replications.
1. $S$ is the set of classification states, where $s \in S$.

### Cases
2. $F_s$ is the set feature sets for state $s$, where $f \in F_s$.
2. $C_s$ is the set of positive cases for state $s$, where $c \in C_s$. The feature set for $c$ is $f_c$.

### Probabilities
1. $p_c$ is the probability that case $c$ occurs in a feature vector under the null hypothesis.
2. $q_{cn}$ is the probability of case $c$ occurring in at least $n$ replications under the null hypothesis.
3. $q_s$ is the probability of $s$ having at least one positive case under the null hypothesis given that there are $\{n_c\}$ cases.
4. $q$ is the probability of at least one state having at least one positive case under the null hypothesis.


## Null Hypothesis
We assume that: (a) features are independent; (b) feature sets are independent (e.g., non-overlapping features); and (c) values in the feature vector are uniformly distributed across $\{-1, 0, 1\}$.

## Derivation of Distribution
First, note that under the null hypothesis,  $p_c =\left (\frac{1}{3}\right )^{|f_c|}$.

Next, we derive $q_{cn}$. This is binomial with parameters $N$ and $p_c$. That is, $q_{cn} = \sum_{j=n_c}^{N} \binom{N}{j} {p_c}^j (1-p_c)^{N - j}$.

Since feature sets are independent, $q_s = \pi_c 1 - \Pi_c (1 - q_{cn})$ given $\{n_c\}$. And, since states are independent, $q = 1 - \pi_s (1 - q_s)$.