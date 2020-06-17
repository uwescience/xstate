# Multi-Classifier
This proposal describes a multi-state classifier that constructs separate classifiers for each state. The problem is motivated by classifying expression states in cell biology. The challenges are:

1. generalizing from one set of conditions to another
2. correlations between features as a result of transcription regulation

## Technical Approach
We propose using several levels of classifier ensembles.
1. A separate ensemble for each state with separate feature selection for each state.
2. Within each state, an ensemble whereby a classifier is described in terms of feature groups, where one gene is chosen from each group to construct a classifier.
3. Classifiers are constructed by training with one holdout for each state. Predictions are probabilistic based on the number of ensemble instances.

## Work Plan
1. Construct basic statistics for
   - single feature classifier. report the accuracy of each gene when it is the sole feature in a classifier.
   - classifier prediction correlation. For pairs of single gene classifiers, report the correlation between their predictions.
   - incremental classifier accuracy. Report the increase in accuracy of using a pair of features instead of the large accuracy of individual features.

1. Evaluate a classifier using regulators as features.
1. Validate the gene equivalence idea.
