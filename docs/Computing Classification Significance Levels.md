# Computing Significance Levels For Case Based Classification
This is about classifying gene expression data. The data consist of columns of genes (features) and rows of instances. Values are in the set $\{-1, 0, 1\}$ that indicate under expression, no differential expression, and over expression. The classification problem is to find the state that best matches a new feature vector.

Case based classification involves counting the (positive) cases for a state $s$ in a feature vector. Typically, an experiment is replicated and so the matching involves multiple feature vectors. The result is a count of $n_{sc}$ occurrences of case $c$ in state $s$. The statistical questions are:

1. If we consider a single state and case, is $n_{sc}$ statistical significant?
2. If we consider a single state, are the cases for $s$ statistically significant?
3. Considering all states, are the results statistically significant?

## Notation


### Cases
2. $F_s$ is the set feature sets for state $s$, where $f \in F_s$.
2. $C_s$ is the set of positive cases for state $s$, where $c \in C_s$. The feature set for $c$ is $f_c$.

### Feature Data for Which Significance Level is Evaluated
1. $N$ is the number of replications.
1. $s \in S$ is a classification state.
2. $C_s^O$ is the observed positive cases for state $s$.
3. $S^O$ is the set of states such that $|C_s^O| >0$ for $s\in S^O$.
2. $n_{sc}$ is the number of occurrences of the positive case $c$. For $c\in C_s^O$, $n_{sc} \geq 1$.

### Probabilities for Significance Levels
1. $p_c$ is the probability that case $c$ occurs in a feature vector under the null hypothesis.
2. $E_{sc}$ is the event that case $c$ in state $s$ occurs in at least $n_{sc}$ replications. Note that for $c \notin C_s^O$, then $n_c=0$ and so $E_{sc}$ is always true.
2. $q_{sc} = P(E_{sc})$ under the null hypothesis.
3. $q_s^O$ is the probability that $E_{sc}$ occurs for at least one $c \in C_s^O$ under the null hypothesis for $s \in S^O$.
4. $q$ is the probability that $E_{sc}$ occurs for at least one $s \in S^O$.


## Null Hypothesis
We assume that: (a) features are independent; (b) feature sets are independent (e.g., non-overlapping features); and (c) values in the feature vector are uniformly distributed across $\{-1, 0, 1\}$.

## Derivation of Probabilities
First, note that under the null hypothesis,  $p_c =\left (\frac{1}{3}\right )^{|f_c|}$.

Next, we derive $q_{sc}$. This is binomial with parameters $N$ and $p_c$. That is, $q_{sc} = \sum_{j=n_c}^{N} \binom{N}{j} {p_c}^j (1-p_c)^{N - j}$.

Since feature sets are independent, $q_s^O = 1-\Pi_{c \in C_s^O} (1 - q_{cn})$.

Since states are independent, $q = 1 - \Pi_{s \in S^O} (1 - q_{C_s^O})$.

## Issues
1. Should *all* positive cases be considered, not just those that are observed?
2. Should there be companion statistics for observed negative cases?