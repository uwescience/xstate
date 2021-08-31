# Naive Bayes Classifier

# Theory

Let $c_i$ be classes and $e_k$ be the $k$-th categorical evidence (e.g., the value of a feature)
and ${\bf e} = (e_1, \cdots, e_K)$
We want to find $P(c_i | {\bf E} )$ for each
$c_i$.

Using the naive Bayes assumption
$P(c_i | {\bf e}) = \frac{\Pi_k P(e_k | c_i) P(c_i)}{P({\bf e})}$.

# Implementation
$P(e_k | c_i)$ is the probability of the feature values occurring in the class.
This is a bit different from the current calculation in that the denominator should be the number of occurrences of the class.

The most profound implication is that the $e_i$ must be conditionally independent and so they should not share features.
An easier way to ensure this is to do the following for each ``case\_collection`` ($c_i$):
1. Calculate $- log(P(e_k | c_i))$ and $- log \left( 1 - P(e_k | c_i) \right)$
1. Sort by descending value so that $k$ indicates this ordering.
1. Delete $e_k$ if it has a feature in commong
with $e_n$, $n < k$.

# Plan
1. ``Case``: Calculate $P(e_k | c_i)$ and save as ``- log``.
1. ``CaseCollection``
   1. ``createIndependentCases``. Use the algorithm above to create cases for each class that do not share features.
1. Ensure that by default there is no filtering by significance level. 
   1. Calculate $log ( P(c_i))$
   1. ``calcLogJointProbability`` calculates
   $\sum_k log v(e_k)$, where $v(e_k) = -log (P(e_k | c_i))$ if $e_k$ is present and $-log (1 - P(e_k | c_i))$ otherwise.
   Provide the actual probability as well.
   1. New bar plots that plot the probability of each class.
   
# Issues
1. Potentially dramatic impact of small probabilities because they result in large negative values.