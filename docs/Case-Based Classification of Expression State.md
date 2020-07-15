# Ensemble Case Based Classification (eCBC) Robustly Identifies Hypoxia in Expression Data for Mycobacterium Tuberculosis
## Motivation
1. Use of classification of expression state
   1. Reproduced experimental conditions
   2. Assess treatments intended to drive pathogen into a desired state
   3. Compare data from different sources
1. By robust, we mean: generalizable, accurate, indicate when don't know.
1. Features that increase robustness:
   1. interpolation than extrapolation
   2. Know which data are used so can evaluate the reliability of those data (e.g., document when used and outcomes)
   3. Augment machine automation with human insight by judgements (e.g., choice of gene combinations used to draw conclusions).

## Classification and Cases
1. Features, instances, class.
1. A classifier is an algorithm that inputs feature data for one or more instances and outputs a class assignment for each instance.
2. Training data. There is typically a pre-step in which the classifier is **trained** on a set of data with known classes so as to estimate parameters of the classification algorithm.
1. Feature set
2. Case. A set of features and instances with the same values of the features that are obtained from training data. Ideally, we want cases where all instances are from the same class. A decision tree can be viewed as constructing cases based on the predicates used to reach a leaf from the root. Cases can also be constructed if features have few values. In expression data, there typically trinary values: -1, 0, and 1.

## Approach
1. Training
   1. Train an ensemble of accurate classifiers using different features.
   2. For each classifier in the ensemble, select cases.
1. Prediction for a feature vector
   1. for each class
      1. for each case for the class
         1. If the case applies to the feature vector, count it
   1. Evaluate the statistical significance of
      1. Each case
      2. Each state

## Considerations
1. Feature Set Filtering. Want to use feature sets that a generalizable in that data taken from one set of experimental conditions are reflective of hypoxia in another.
   1. We selected feature sets that enrich GO term for hypoxia, lipid metabolism, and fatty acid production.
1. Bogus cases. With a large number of features and small number of instances, can have cases with high significance level due to chance.
   1. Consider $N_I$ instances of which $N_C$ are the positive class, $N_F$ features, $p$ the probability of two feature values matching.
   2. The probability that all features match the feature vector is $q = p^{N_C}$.
   3. Consider the probability that at least one of the instances of the positive class matches the feature vector and no negative class matches the feature vector. This would result in a significance level of 0.
      1. The probability of at least one of the positive class instances matching is $1 - (1 - q)^{N_C}$.
      2. The probability of also having no negative instance matching is $P = [1 - (1 - q)^{N_C}](1-q)^{N_I-N_C} = (1 -q)^{N_I}[1 - (1-q)^{N_C}]$.
      3. As $N_F$ gets large, $q \rightarrow 0$

## Contributions of eCBC
1. Novel approach to constructing cases using an ensemble of classifiers.
2. Demonstrate generalizing to other datasets (Galagan). Existing studies do not consider the challenge of generalizing from biological data collected under one set of conditions to classifying data collected under different conditions. We address this by: (a) using regulators as features since their role may be more consistent across conditions; (b) using an ensemble of classifiers.
3. Evaluating the sigificance of a class based on the positive cases for the class.
4. Selecting the ensemble
   1. Focus on performance characteristics of classifiers, not similarity of feature values in training data so ensure an impact on prediction.

## Related Work
1. Case Based Reasoning for Medical Informatics, Periklis Andritsos, Igor Jurisica, Janice I. Glasgow.
   1. Examples of use CBR in medicine
   2. Elements: Cases (features, text), retrieval engine, distance algorithm, case mainenance
1. Case-Based Retrieval Framework for Gene Expression Data, Annassi.
   1. Connects classification with CBR. But objective is to retrieve the "best case", not to find cases that support class prediction.
   2. Unclear that the data contain class information. Rather, this is a case retrieval system.
   3. Discusses a separate step for training (construction of the case base).
1. A balanced iterative random forest for gene selection from microarray data. Annassi.
   1. Uses a "balanced" variation of random forests on microarray data to train a classifier.
1. Evaluation of Gene Expression Classification Studies: Factors Associated with Classification Performance. Putri W. Novianti.
   1. Proposes that best choice depends on the data. This raises a question about how to generalize from one dataset to another.
1. Galagan and Sherman have data that provide a basis for evaluating generalizability.
2. Decision trees.
   1. Do leaves of Tree specify a case?
   2. No because there could be many values assigned to features.
   4. Could argue that leaves are cases if the minimum number of instances for a node is 1 so that all leaves are class homogeneous. Then, a node can be viewed as an aggregation of cases. But still need to identify "don't know" cases, those with feature values that are not present in the data but would be classified in the Tree.