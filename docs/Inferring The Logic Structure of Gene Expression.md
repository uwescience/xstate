# Inferring The Logic Structure of Gene Expression State

## Motivation
1. Simplify the description of state classifications.

## Approach
1. Do SVM classification for a set of genes whose expression states are -1 (underexpresed), 0 (normal), 1 (overexpression).
2. The parameters parameters of the SVM allow for inferring a ternary truth table.
3. Evaluate which ternary logic is most approprite.
4. Can use observational counts to rank how well different logic expressions match the data.
5. To implement matching, consider explicitly evaluating all possible expressions. For binary expressions:
  1. There are two possibilities for each variable, $A$ and $\bar{A}$.
  2. There are four conjuctions of two variables: $A \land B$, $\bar{A} \land B$, $A \land \bar{B}$, and $\bar{A} \land \bar{B}$. So, if there are $n$ boolean variables, then there are $2n$ possible singleton terms in disjunctive normal form (DNF), $n \choose 2$ pairs of variables, each of which has 4 variants. Hence, $4^{n \choose 2}$ binary terms. Thus, the base 2 logorithm of the size of the set possible DNF is no larger than $2n + 4^{n \choose 2}$. For $n=1$, this is 2. For $n=2$, this is 256. For $n=3$, this is $\approx 10^{20}$.