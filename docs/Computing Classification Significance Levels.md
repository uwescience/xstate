# Computing Significance Levels For Case Based Classification
This is about classifying gene expression data. The data consist of columns of genes (features) and rows of instances. Values are in the set $\{-1, 0, 1\}$ that indicate under expression, no differential expression, and over expression. The classification problem is to find the state that best matches a new feature vector.

Case based classification involves counting the (positive) cases for a state $s$ in a feature vector. Typically, an experiment is replicated and so the matching involves multiple feature vectors. The result is a count of $n_{sc}$ occurrences of case $c$ in state $s$. The statistical questions are:

1. If we consider a single state and case, is $n_{sc}$ statistical significant?
2. If we consider a single state, are the cases for $s$ statistically significant?
   1. What is the probability of at least one of the positive cases occurring?
   2. What is the probability of at least $n$ positive cases occurring?
3. Considering all states, are the results statistically significant?

## Notation


### Cases
1. $F^{\star}$ is the set of all features $f$.
2. $\cal{F}_s$ is the set feature sets for state $s$, where $F^{\star} \supseteq F \in \cal{F}_s$.
3. A case $c$ has an associated feature set $F_s \in \cal{F}_s$.
4. $v_{cf} \in \{-1, 0, 1\}$ is the value of feature $f$ for case $c$.
5. $s \in S$ is a classification state.
3. $C_s$ is the set of cases for state $s$. $C_s^P$ is the set of positive cases for state $s$, where $c \in C_s$. $C_s^N$ is the set of negative cases for state $s$. $C_s^P, C_s^N$ partition $C_s$.

### Data
1. $M$ is the number of replications.
2. $n^P_s$ is the number of occurrences of the positive cases in state $s$.
3. $n^N_s$ is defined analogously.

### Probabilities for Significance Levels
2. $E_{c}$ is the event that case $c$ occurs in a replication.
2. $q_c = P(E_{c})$ under the null hypothesis.
4. $Q_{sk}^P$ is the probability that there are at least $k$ positive case occurrences in a replication.
5. $Q_{sm}^N$ is the probability that there are at least $m$ negatuve case occurrences in a replication.
6. $\rho_c$ is the significance level of case $c$ across replications.
7. $\rho_s$ is the significance level of all cases for a state.
6. $\rho_s^P$ is the significance level for positive cases for state $s$.
7. $\rho_s^N$ is the significance level for negative cases for state $s$.


## Null Hypothesis
We assume that: (a) features are independent; (b) feature sets are independent (e.g., non-overlapping features); and (c) values in the feature vector are uniformly distributed across $\{-1, 0, 1\}$.

## Derivation of Probabilities
First, note that under the null hypothesis,  $q_c =\left (\frac{1}{3}\right )^{|F_c|}$.

Note that $\rho_c$ is binomial with parameters $M$, $m$, and $q_c$. That is, $\rho_c = \sum_{j=n_c}^{M} \binom{M}{j} {q_c}^j (1-q_c)^{M - j}$.

Now consider $Q_{sk}^P$. Clearly, $Q_{sk}^P = 0$ if $k > |C_s^P|$.
Otherwise, we calculate the event combinatorics.

Denote the power set of $C$ as $2^C$, and let $2_k^C \subseteq 2^C$
such that all subsets are of size $k$.
Let $C^{\prime} \subseteq C$ be a set of cases.
Then, the probability that exactly $k$ of these cases occur in a single replication is
$P(s,k,C^{\prime}) = \sum_{C \in 2_k^{C^{\prime}}} \Pi_{c \in C} \rho_c \Pi_{c \in C^{\prime} - C} (1 - \rho_c)$.
$Q^P_{sk} = \sum_{j=k}^{|C_s^P|} P(s, k, C_s^P)$.
$Q^N_{sk}$ can be calculated analogously.

The significance level for positive cases for state $s$ is
the probability that each replication has at least as many positive cases
as what was observed under the null hypothesis.
That is, $\rho_s^P = \Pi_m Q^P_{s, n_s^P}$.
Similarly, $\rho_s^P = \Pi_m Q^N_{s, n_s^N}$.

## Issues