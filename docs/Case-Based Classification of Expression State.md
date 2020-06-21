# Case-Based Classification of Expression State
## Motivation
1. More reliable to do interpolation than extrapolation
2. Want visibility as to the data used in decision making

## Approach
1. Separate classifier for each state.
2. Feature selection based on accuracy achieved from classification
3. Use trinary valued features
4. Rank feature sets for each classifier
5. For each feature set, calculate the following for each combination of trinary valued features:
   1. Number of positive cases
   2. Prob of at least that many positive cases if 0.5 probability of class given the assignment of trinary values to features.
   3. Prob of fewer than this many negative cases under the same conditions.
1. Classifying a new feature vector
   1. For each binary classifier
      1. For each feature set
         1. Find significance level for positive and negative cases
      1. Aggregate significance levels for positive and negative class
   1. Report significance levels by classifier