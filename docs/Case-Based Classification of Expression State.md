# Classification of Gene Expression State Using Case Distributions
## Motivation
1. More reliable to do interpolation than extrapolation
2. Want visibility as to the data used in decision making
3. Key terms
   1. Feature set
   2. Case

## Approach
1. Find cases
   1. Find feature sets with disjoint features (so have disjoint cases)
   2. Select cases with high significance levels
1. evaluate(feature_vector)
   1. for each state
      1. Find the significance levels for all applicable cases
      1. Plot distribution
   1. Choose state with mass of distribution at higher significance levels

## Considerations
1. Bogus cases. With a large number of features and small number of instances, can have cases with high significance level due to chance.
   1. Consider $m$ instances, $n$ features, $c$ class instances.
   2. Assume only have over expression or no differential expression. 0.5 probability of each.
   3. Probability of a random feature vector matching a positive a single case is $c(0.5^n)(1-0.5^n)^{m-1}$.