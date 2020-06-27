# Classification of Gene Expression State Using Case Distributions
## Motivation
1. More reliable to do interpolation than extrapolation
2. Want visibility as to the data used in decision making
3. Key terms
   1. Feature set
   2. Casex

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
   1. Consider $N_I$ instances of which $N_C$ are the positive class, $N_F$ features, $p$ the probability of two feature values matching.
   2. The probability that all features match the feature vector is $q = p^{N_C}$.
   3. Consider the probability that at least one of the instances of the positive class matches the feature vector and no negative class matches the feature vector. This would result in a significance level of 0.
      1. The probability of at least one of the positive class instances matching is $1 - (1 - q)^{N_C}$.
      2. The probability of also having no negative instance matching is $P = [1 - (1 - q)^{N_C}](1-q)^{N_I-N_C} = (1 -q)^{N_I}[1 - (1-q)^{N_C}]$.
      3. As $N_F$ gets large, $q \rightarrow 0$