# Optimal Stress Seriation 

## Usage

**seriation_mathopt.py** contains all the functions to compute the optimal stress seriation with different methods.

In the class, _MatrixSeriation_, the parameters are set:

seriation_solver = MatrixSeriation(matrix, file=file, method=method, symmetric_ordering=sym, eps_neigh=eps)

where:

    - _matrix_: is the matrix to be seriated
    - file (optional): name of the file (for the expriments to track the obtained results).
    - _method_: 
        * tsp: Hamiltonian path approaches.
        * spp: shortest path approaches
        * s: models based on the "s"-variables.
        * general: four index formulation (default)
    - _symmetric_ordering_; True: if both rows and columns are sorted with the same permutations; False: if rows and columns are differently permuted.
    - eps_neigh: 1 (Default): if the von Neumman neighbor is used, 1.5: is the Moore neighborhood is used for the stress seriation.

the function:

optimal_order = seriation_solver.solve()

computes the re-arranged matrix with the optimal seriation (with the indicated parameters). The information can be extracted with _optimal_order.ATTRIBUTE_ where ATTRIBUTE is one of the following:

    - nrows: number of rows of the matrix
    - ncols: number of columns of the matrix
    - order_rows: list of sorted rows with the optimal permutation
    - order_columns: list of sorted columns with the optimal permutation.
    - transformed_matrix: re-arranged matrix with the optimal permutations.

It is possible to extract extra information about the solution procedure by calling: seriation_solver.solution_info. Specifically, apart from the above one obtains:

    - CPUTime
    - ObjectiveValue
    - MIPGap
    - LPRelaxation at Root
    - Work
    - NodeCount
    - Memory
    - NumBinVars
    - NumCtrs
    - Homogeneity Original
    - Homogeneity Transformed
    - Error Original
    - Error Transformed

The heatmap matrices (original and transformed) are ploted with the following command: **seriation_solver.plot_matrices**




## Instances:

- *easy* instances: square matrices nxn designed to contain structures that facilitate ordering. To generate them, we draw $n$ points in [0,100] and output their distance matrix (pairwise absolute differences). Folder *easy_instances*. Five random instances for each n in [10,20,30,40,50,100].
- structured *square* instances: Square matrices (nxn) but contain more complex structured patterns. They are generated also as pairwise Euclidean distances between randomly generated points in [0,100]x [0,100]. Folder *easy_instances*.
- *nonsquare* instances: For n > m, we generate rectangular matrices (nxm) created with structured patterns that introduce complexity in ordering tasks. They are generated also as pairwise Euclidean distances between two sets of (n and m) randomly generated points in [0,100]x [0,100]. Folder *nonsqr_instances*.
- *binary square* instances: We randomly generate binary square instances (nxn), but with certain densities (percent of elements that are 1 in th whole matrix). We chose densities d in {25%,50%,75%}. Folder *bin_instances* with pattern *bin{d}\_{n}\_instance.txt*.
- *binary nonsquare* instances: Rectangular matrices (n x m) and contain only binary values, and with the same density parameters than the binary square instances.