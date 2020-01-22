# State Characterization Algorithms

## Considerations

## Approach 1: Minimum Intra-State Distance, Maximum Inter-State Distinction (MIDMID)

**Construct Candidate Groupings**

1. For state in states

   1. geneClusters = hierarchicalCluster(vectors(state))
   2. breakpoint = distance at which there is a "knee" in number clusters
   3. groupings = groups(geneClusters, breakpoint)


**Refine Groupings**

1. for pairs state1, state2
   1. for group1 in groupings(state1), group2 in groupings(state2)
      1. If intersectionScore(group1, group2) is large
         1. Discard group1 from groupings(state1), group2 from groupings(state2)
      