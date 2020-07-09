# Computing Classification Significance Levels
Classification involves matching cases against replications of an experiment. The result of this matching is that for state $s$ there are $n_{fs}$ occurrences of features set $f$. Are the $\{n_{fs}\}$ statistically significant?

## Null Hypothesis
We assume that the feature vector is randomly distributed so that for the $k$-th feature set $f_1,\cdots,f_{N_k}$, the probability of a feature vector satisfying a case is $q=\left (\frac{1}{3}\right )^{N_k}$.

## Derivation of Distribution
We want to find $Q = P(n_1 \geq n_{c_1}, \cdots, n_C \geq n_{c_C})$, the significance level of the occurrence of positive cases.

Given the positive cases $c_1,\cdot, c_{N_C}$ that occur with $n_c$ times and have the feature set $f_c$. Let $N_R$ be the number of replications. Then, the probability of having at least $n_c$ cases is $Q_c = \sum_{j=n_c}^{N_R} \binom{N_R}{j} q_{f_c}^j (1-q_{f_c})^{N_R - j}$.
$Q_c$ is the significance level for the occurrence of a case.

Since feature sets are disjoint, $Q = \Pi_c Q_c$.