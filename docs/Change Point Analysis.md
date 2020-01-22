# Change Point Analysis
This procedure describes how to detect change points in time course expression data and identify the genes that characterize the change points.

1. Normalize the data.
   1. Convert to transcript per kilobase.
   2. Convert to fraction of reads at time point.
   3. Calculate average values across replicas at the reference time.
   4. Compute log 2 ratio w.r.t. reference value. (How handle 0 references?)
1. Discretize the data. Using a threshold (e.g., 1, -1), convert expression values into the {-1, 0, 1}.
1. Aggregate replicas, computing: (a) fraction upregulated; (b) fraction down regulated. Check for instances of non-zero values of both (a) and (b).
1. Cluster rows based on absolute values of aggregated discretized expression levels.
1. Look at differences in values to find change points.