# Optimal Stress Seriation 

## Instances:

- *easy* instances: square matrices nxn designed to contain structures that facilitate ordering. To generate them, we draw $n$ points in [0,100] and output their distance matrix (pairwise absolute differences). Folder *easy_instances*. Five random instances for each n in [10,20,30,40,50,100].
- structured *square* instances: Square matrices (nxn) but contain more complex structured patterns. They are generated also as pairwise Euclidean distances between randomly generated points in [0,100]x [0,100]. Folder *easy_instances*.
- *nonsquare* instances: For n > m, we generate rectangular matrices (nxm) created with structured patterns that introduce complexity in ordering tasks. They are generated also as pairwise Euclidean distances between two sets of (n and m) randomly generated points in [0,100]x [0,100]. Folder *nonsqr_instances*.
- *binary square* instances: We randomly generate binary square instances (nxn), but with certain densities (percent of elements that are 1 in th whole matrix). We chose densities d in {25%,50%,75%}. Folder *bin_instances* with pattern *bin{d}\_{n}\_instance.txt*.
- *binary nonsquare* instances: Rectangular matrices (n x m) and contain only binary values, and with the same density parameters than the binary square instances.