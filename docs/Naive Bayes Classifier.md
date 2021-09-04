# Probabilistic Classifiers

This writeup explores the use of Naive Bayes classification and some extensions.

## Naive Bayes

Let $c_i$ be a class, and $e_k$ be the $k$-th categorical evidence (e.g., the value of a feature)
with ${\bf e} = (e_1, \cdots, e_K)$.
Naive Bayes operates by calculating
We want to find $P(c_i | {\bf E} )$ for each
$c_i$.
The chosen class is the one with the largest probability
of occurrence.

Using the naive Bayes probability formula,
$P(c_i | {\bf e}) = \frac{\Pi_k P(e_k | c_i) P(c_i)}{P({\bf e})}$.

## Choosing Evidence

By *feature*, we mean an attribute of the data.
By *evidence*, we mean some transformation of
feature values so that defines a boolean function
of a feature vector (e.g., the feature vector
satisfies a case).

We choose cases based on their ability to discriminate between classes.
Each case provides evidence in that
the presence (or absence) of the case in
feature data ${\bf x}$
supports the class.
Define
$e({\bf x})$ (hereafter just denoted by $e$)
to be boolean, with a value of $1$ if the case
is supportive of the class in ${\bf x}$.
Since ${\bf x}$ is assumed to be selected randomly (from a stationary distribution),
$e$ is a Bernoulli variable.
In the following, $c$ be a binary variable that is 1 if the
sample is from the class that $e$ supports and is 0 otherwise.

We want to find cases such that either
$r = \frac{P(e = 1 | c = 1)}{P(e = 0 | c = 0)}$
is large.
Note that $0 \leq r \leq \infty$.

## Extended Binomial Naive Bayes
Binomial Naive Bayes conditions on the occurrence of a set
of $K$ events.
A generalization is to condition on the occurrence of at
least $M$ of the $K$ events.
That is, we calculate
$P(c_i | \mbox{at least }$M$ \mbox{ of } {\bf e})$.

## Implementation Considerations
There are two considerations:
1. Select cases
1. Calculate $P(c_i | \mbox{at least }$M$ \mbox{ of } {\bf e})$.

### Selecting Cases
The $e_i$ must be conditionally independent and so their
underlying cases should have disjoint features.
Let $v_k = 1$ if case $k$ is present;
otherwise $v_k = 0$.

Below is a
"greedy'' algorithm for selecting $K$ cases
that have no feature in common
from an initial set of cases.
1. Calculate $r_k = max \left[ 
\frac{P(v_k = 1 | c = 1)}{P(v_k = 0 | c = 0)},
\frac{P(v_k = 0 | c = 1)}{Pv_k = 1 | c = 0)}
\right]$
1. Sort by descending $r_k$,
and let $n$ index this ordering.
1. Delete $e_j$ if its associated case has a feature in common
with the case for any $e_{j ^{\prime}}$, $j^{\prime} < j$.
1. Choose the top $K$ cases.

### Calculating Probabilities
Since the probability of events have different probabilities,
probabilities are calculated numerically using
combinatorics.
We expect that a typical value for $M$ will be 2 or 3.
So, it will be easier to calculate
$1 - P(c_i | \mbox{fewer than }$M$ \mbox{ of } {\bf e} \mbox{ are not present})$.

## Plan
1. Try off-the-shelf naive Bayes with all features.

1. Use cases
    1. ``Case``: Calculate $P(e_k | c_i)$ and save as ``- log``.
    1. ``CaseCollection``
       1. ``createIndependentCases``. Use the algorithm above to create cases for each class that do not share features.
    1. Ensure that by default there is no filtering by significance level. 
       1. Calculate $log ( P(c_i))$
       1. ``calcLogJointProbability`` calculates
       $\sum_k log v(e_k)$, where $v(e_k) = -log (P(e_k | c_i))$ if $e_k$ is present and $-log (1 - P(e_k | c_i))$ otherwise.
       Provide the actual probability as well.
       1. New bar plots that plot the probability of each class.
   
## Notes
1. ``sklearn`` has an implementation of [Bernoulli Naive Bayes](https://scikit-learn.org/stable/modules/generated/sklearn.naive_bayes.BernoulliNB.html).

## Issues
1. Potentially dramatic impact of small probabilities because they result in large negative values.