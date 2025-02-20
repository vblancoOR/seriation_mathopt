from seriation_mathopt import *
from sklearn.preprocessing import MinMaxScaler
scaler = MinMaxScaler(feature_range=(0, 1))

results_df = pd.DataFrame()
NN=[10, 20, 30, 40, 50, 100]
MM=[10, 20, 30, 40, 50, 100]
for method in ["tsp", "spp", "s", "general"]:
    if method=="s":
        EPS=[1]
    else:
        EPS =[1, 1.5]
    for eps in EPS:
        for n in NN:
            for m in MM:
                if n>=m:
                    for it in range(5):
                        if n==m:
                            SYM=[True, False]
                        else:
                            SYM=[False]
                        for sym in SYM:
                            for type in ["bin", "easy", "sqr", "nonsqr"]:
                                if type=="bin":
                                    D=[25,50,75]
                                else:
                                    D=[1]
                                for dens in D:
                                    if n==m:
                                        if type=="bin":
                                            file=f"{type}{dens}_{n}_{it+1}.txt"
                                        else:
                                            file=f"{type}_{n}_{it+1}.txt"
                                    else:
                                        if type=="bin":
                                            file=f"{type}{dens}_{n}_{m}_{it+1}.txt"
                                        else:
                                            file=f"{type}_{n}_{m}_{it+1}.txt"
                                    try:
                                        test_matrix=np.loadtxt(f"{type}_instances/"+file, delimiter=",")
                                        if type!="bin":
                                            test_matrix=scaler.fit_transform(test_matrix)
                                        seriation_solver = MatrixSeriation(test_matrix, file=file, method=method, symmetric_ordering=sym, eps_neigh=eps)
                                        ##
                                        print(f"type: {type}, n: {n}, m: {m}, Method: {method}, Sym: {sym}, eps: {eps}, it: {it+1}")
                                        optimal_order = seriation_solver.solve()
                                        print("Solution Info:", seriation_solver.solution_info)
                                    # Convert result to DataFrame
                                        new_row = pd.DataFrame([seriation_solver.solution_info])
                                        # Append result to the main DataFrame
                                        results_df = pd.concat([results_df, new_row], ignore_index=True)
                                        # Check if file exists to control the header writing
                                        file_exists = os.path.isfile("test_seriation.csv")
                                        # Save to CSV iteratively, only write header if file does not exist
                                        new_row.to_csv("test_seriation.csv", mode="a", index=False, header=not file_exists)
                                    except:
                                        print(f"{file} does not exist.")
