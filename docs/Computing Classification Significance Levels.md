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
3. $C_s$ is the set of cases for state $s$. $C_s^P$ is the set of positive cases for state $s$, where $c \in C_s$. $C_s^N$ is the set of negative cases for state $s$. $C_s^P, C_s^N$ partition $C_s$.
3. A case $c$ has an associated feature set $F_s \in \cal{F}_s$. 
4. $v_{c,f}$ is the value of $f$ for $c$.

### Feature Data for Which Significance Level is Evaluated
1. $M$ is the number of replications.
1. $s \in S$ is a classification state.
2. $C_s^O$ is the observed positive cases for state $s$.
3. $S^O$ is the set of states such that $|C_s^O| >0$ for $s\in S^O$.
2. $n_{sc}$ is the number of occurrences of the positive case $c$. For $c\in C_s^O$, $n_{sc} \geq 1$.

### Probabilities for Significance Levels
1. $p_c$ is the probability that case $c$ occurs in a feature vector under the null hypothesis.
2. $E_{sc}$ is the event that case $c$ in state $s$ occurs in at least $n_{sc}$ replications. Note that for $c \notin C_s^O$, then $n_c=0$ and so $E_{sc}$ is always true.
2. $q_{sc} = P(E_{sc})$ under the null hypothesis.
3. $q_s^O$ is the probability that $E_{sc}$ occurs for at least one $c \in C_s^O$ under the null hypothesis for $s \in S^O$.
4. $Q_s$ is the probability that there are at least $n_s = \sum_{c \in C_s^P} n_{sc}$ positive case occurrences.
4. $q$ is the probability that $E_{sc}$ occurs for at least one $s \in S^O$.


## Null Hypothesis
We assume that: (a) features are independent; (b) feature sets are independent (e.g., non-overlapping features); and (c) values in the feature vector are uniformly distributed across $\{-1, 0, 1\}$.

## Derivation of Probabilities
First, note that under the null hypothesis,  $p_c =\left (\frac{1}{3}\right )^{|f_c|}$.

Next, we derive $q_{sc}$. This is binomial with parameters $M$ and $p_c$. That is, $q_{sc} = \sum_{j=n_c}^{M} \binom{M}{j} {p_c}^j (1-p_c)^{M - j}$.

Since feature sets are independent, $q_s^O = 1-\Pi_{c \in C_s^O} (1 - q_{cn})$.

Since states are independent, $q = 1 - \Pi_{s \in S^O} (1 - q_{C_s^O})$.

Let $g(C, n)$ be a randomly chosen subset of size $n$ and $\cal{C}_n$ be the set of all such subsets. Then, the probability of exactly $n$ positive cases in a replication is $Q_{sn} = \sum_{C \in g(C, n)} \Pi_{c \in C} p_{sc} \Pi_{c \in g(C, n) - C} (1 - p_{sc})$, and so the probability of at least $n$ is 
$\sum_{n = {n_s}}^{|C_s^P|} Q_{sn}$.

A partition of an integer $K$ is a set of non-negative numbers $\{k_i\}$ such that $\sum_i k_i = K$. Let $\cal{P}(K, m)$ be the set of partitions of $K$ of size $j$.
Then, $Q_s = \sum_{m=1}^M \sum_{P \in \cal{P}(n_c, m)$ FINISH 

## Issues
1. Should *all* positive cases be considered, not just those that are observed?
2. Should there be companion statistics for observed negative cases?