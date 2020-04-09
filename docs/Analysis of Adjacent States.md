# Analysis of Adjacent States

This is an analysis of state assignment prediction for
classification of time serial observations.
We consider misclassifications, and investigate if misclassifications are for a state that is adjacent in time.
Such confusions might be expected in a dynamical system.

Let $t$ be the time index of observations, $s\in S$
be the set of states, and $c(t)\in S$
be the state at $t$.
Let $t_b = max_{u < t} \{ c(u) \neq c(t) \}$, and define the **before state**
as $b(t) = s(t_b)$.
Similarly, let $t_a = min_{u > t} \{ c(u) \neq c(t) \}$,
and define the **after state** as
$a(t) = c(t_a)$.
We define $A_t = \{b(t), a(t)\}$ to be the adjacent states of $t$.

Let $p(t)$ be the predicated state at time $t$ for a random classifier that only uses
the distribution of states to do class assignment.
That is, the probability that of assigning $s$ to observation $t$ is
$q^S_s = \frac{\sum_t c(t) }{|S|}$.
We want to calculate $q^A_t$, the probability that $p(t)$
is in $A_t$ given that $p(t) \neq c(t)$.
Since we use a random classifier, 
$q^A_t = P(c(t) \in A_t| p(t) \neq c(t)) = \frac{q_{b(t)} + q_{a(t)} }{1 - q_{p(t)}}$.

We are given observations at times $t_1, \cdots, t_N$ such that $c(t) \neq p^{\prime}(t)$ for
some classifier.
Suppose that for $M$ of these events we have that $p(t) \in {b(t), a(t)}$. We want to
assess statistical significance.

Our null hypothesis is that that before and after states are assigned at random.
Thus, the probability of $M$ such events is the sum of the probabilities of the
probability of all combinations of $M$ of these events occurring with probability
$q_t$ and $N-M$ events
not occurring with probability $1 - q_t$.
Suppose that misclassifications at times $t_1, \cdots, t_M$ are adjacent and
those at $t_{M+1}, \cdots t_N$ are not.
Let $T_1 =  \{ t_1, \cdots, t_M \}$ and $T_2 = \{ t_{M+1}, \cdots, t_N \}$.
The probability of these joint events is
$\prod_{t \in T_1} q_t \prod_{t \in T_2} (1 - q_t)$.
