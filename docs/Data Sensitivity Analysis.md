# Data Sensitivity Analysis
This write-up addresses concerns about the extent to which a feature vector may be inappropriate for a classifier We use the term **sensitive data items** to refer to samples whose presence in the test data to changes predictions to a different class from that expected by the training data. We consider two cases for sensitive data items: (a) the test data are drawn from a different distribution from the training data and (b) the classification large rests on a small subset of features in the test data. We refer to (a) as **distribution sensitivity** and (b) as **feature sensitivity**.

## Feature Sensitivity

Item (b) is largely addressed by the ensemble approach since each classifier in the ensemble is trained on a subset of data. If the classification depends on a single value, then a subset without this value predicts a different class from a subset with the value.

We can quantify the foregoing as follows. Let $n$ be the number of data samples, of which $n-k$ are used in training. We assume that the number of classifiers is considerably larger than $n$. Assume that there is a single sensitive data item. Then, the probability that a training set does *not* contain this data item is $p=(1-\frac{1}{n})^k$. Further assume that the non-sensitive data predict class $c$ and the presence of the sensitive data item yields a prediction of $c\prime$. So, the presence of a sensitive data item is indicated if predictions are split between two classes with probabilities $p$ and $1-p$. In our analysis, $n=24$ and $k=1$, and so the probabilities are approximately $p=0.95$ and $1-p=0.05$.

We can generalize this by considering that there are $m$ sensitive data items where $m << n$. We assume that the presence of any of the sensitive data items changes the classification. So, now the probability of *not* including a sensitive data item is $p=1-(1-\frac{m}{n})^k$ (assuming that $m<<n$). With this formula and given a split between two classes with probabilities $q_1$ and $q_2$, we can find a value of $m$ that would be required to produce this split.

## Distribution Sensitivity

1. Distribution sensitivity should be weighted by feature importance.
2. How do feature sensitivity and distribution sensitivity differ from each other?