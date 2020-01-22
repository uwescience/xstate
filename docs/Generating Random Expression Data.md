# Generating Random Expression Data
Consider  $X = \{x_{m,n} \}$ where $m \leq M$ is the number of instances and $n \leq N$ is the number of features. In our case, there are $N$ genes and $M$ measurements. We want to generate a new matrix is which $f$ of the values are changed in the following way. $x_{m,n}$ is unchanged with probability $1-f$. With probability $f$, $x_{m,n}$ is drawn uniformly from $\{x_{1,n}, \cdots, x_{M,n}\}$.

## Algorithm
The algorithm must have the following characteristics:

1. preseve the marginals for each feature;
2. preserve covariances for $f = 0$, and eliminate covariances when $f=1$;
3. be computationally efficient.

Given $X$, $K \times N$, the approach is:

1. $X^{\prime}_0$ is obtained by sampling $M$ rows from $X$ with replacement.
1. $Y$ is $M \times N$. It contains values drawn from a Bernoulli process with parameter $f$.
1. $Z$ is $M \times N$, where $x_{m,n}$ is draw from the marginal distribution of feature $n$.
1. $X^{\prime} = X^{\prime}_0 - X^{\prime}_0 * Y + Z * Y$, where $*$ means element wise multiplication.
        
The following tests should be used to evaluate $X^{\prime}$:
 
 - Distribution of the marginals is the same as the original empirical distribution.
 - Correlations should decrease with $f$.
 - Dataframes have the correct shape.       
  

## Visualizations

1. To visualize marginal distributions for a large number of features, do a heatmap with:
   - x-axis is value (-1, 0, 1)
   - y-axis is gene, sorted by $P(-1)$ then $P(0)$ and then $P(1)$
   - Value is $P(x)$

  The function should return the ordering of the genes and there should be an
option to force the same order.

1. To compare marginal distributions $A$ and $B$ for $X^A$, $X^B$, do a line plot
   - x axis is $P(a)$ y-axis is $P(b)$
   - points are $(P_A, P_B)$ for $x^A_{m,n}$, $x^B_{m,n}$.

1. Correlations are depicted with a heatmap that only plot correlations within the Barlett limit. Returns the genes plotted.


## Implementation Notes

1. Want to have a total of about 1500 instances balanced across the sates ([reference says the number of instances should be about the same as the number of uncorrelated features](https://academic.oup.com/bioinformatics/article/21/8/1509/249540)).

1. Create a notebook with analysis of data generation.


