# Scoring the Results of Conditions-Based Classification and Assessing its Statistical Significance

This procedure is motivated by a need to assess the significance of outcomes of multiple methods for constructing classifiers and preparation of test data.

## Scoring

An ideal result is that the predicted stage (class) is the same within a condition, and the predicted stage under distinct conditions differ. We can construct a score that quantifies this idea as follows.
Let $y_{cr}$ be the class predicted for the $r^{th}$ replication of the $c^{th}$ class.
Suppose there are $C$ conditions, each with $R$ replications.
We want the mean sum of squares *within* conditions,
$W = \frac{1}{R(C-1} \sum_c \sum_r (y_{rc} - \bar{y}_c)^2$ (where $\bar{y}_c$ 
is the mean value within condition $c$), to be small.
Also, we want the mean sum of square *between* conditions, 
$B = \frac {1}{C - 1}\sum_c (\bar{y}_c - \bar{y})^2$,
to be large. So, we score a plot by the ratio $\frac{B}{W}$.
To handle cases in which $W = 0$, we use $\tilde{W}$ in the denominator, where $\tilde{W} = max(0.1, W)$.

## Statistical Significance

Now consider the statistical significance of obtaining a large score.
We consider the following situation.
An experiment produces one of $K$ outcomes.
An experiment is conducted under a
set of conditions;
there are $R$ replications of each
condition, and $C$ conditions in
a study.
We consider $M$ study variations,
essentially other conditions that
apply to all conditions in the study variation.

An **event** for a condition of a study occurs when
there are at least $n$ of the
$R$ replications have the same outcome.
We want to find the probability of at least one study having
such an event for all of its conditions.

Define the following:

  * $K$ number of distinct outcomes for an experiment
  * $R$ number of replications for a condition
  * $C$: number of conditions in a study
  * $M$: number of study variations
  * $p_k$: probility of outcome $k$ under the null hypothesis
  * $BT(N, p, n)$: tail probability of binomial, $P(\tilde{n} \geq n; N, p)$.

Let $P_k (n)$ be the probability of 
an event for outcome $k$;
that is, there is
at least $n$ outcomes of type $k$ for a single condition within a study. $P_k (n) = BT(p_k, R, n)$.

Let $\tilde{P} (n)$ be the probability of an event for a condition.
We assume that $n$ is large enough so that there can be an event for at most
one outcome.
Thus, $\tilde{P} (n) \approx \sum_k P_k (n) \prod_{l \neq k} (1 - P_l)$.
Also, let $\cal{K}$ be a subset of the classes.
We consider this probability if only a subset of classes are considered.
That is, $\tilde{P} (n, \cal{K}) \approx \sum_{k \in \cal{K}} P_k (n) \prod_{l \neq k, l \in \cal{K}} (1 - P_l)$.

Let $\tilde{Q} (n)$ be the probability that an event occurs for
all conditions in a study.
Then, $\tilde{Q} (n) = \left( \tilde{P} (n) \right)^C$.
Let $Q (n)$ be the probability that there is at least one study of the
study variations such that an event occurs for all of its conditions.
Then $Q(n) = 1 - \left( 1 - \tilde{Q} \right)^M$.
