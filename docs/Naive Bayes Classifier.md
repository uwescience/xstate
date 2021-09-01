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

## Implementation
$P(e_k | c_i)$ is the probability of the feature values occurring in the class.
This is a bit different from the current calculation in that the denominator should be the number of occurrences of the class.

The most profound implication is that the $e_i$ must be conditionally independent and so they should not share features.
An easier way to ensure this is to do the following for each ``case\_collection`` ($c_i$):
1. Calculate $- log(P(e_k | c_i))$ and $- log \left( 1 - P(e_k | c_i) \right)$
1. Sort by descending value so that $k$ indicates this ordering.
1. Delete $e_k$ if it has a feature in commong
with $e_n$, $n < k$.

## Plan
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
1. ``sklearn`` has an implementation of naive bayes.

## Issues
1. Potentially dramatic impact of small probabilities because they result in large negative values.