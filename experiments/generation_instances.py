from seriation_mathopt import *

#Generate random matrices
NN=[10, 20, 30, 40, 50, 100]
MM=[10, 20, 30, 40, 50, 100]
for n in NN:
    for it in range(5):
        test_matrix = generate_easy_instances(n)#np.round(np.random.rand(n, m),decimals=0)  # Example similarity/distance matrix
        file=f"easy_{test_matrix.shape[0]}_{it+1}.txt"
        np.savetxt("easy_instances/"+file, test_matrix, fmt="%.2f", delimiter=",")

        test_matrix = generate_sqr_instances(n)#np.round(np.random.rand(n, m),decimals=0)  # Example similarity/distance matrix
        file=f"sqr_{test_matrix.shape[0]}_{it+1}.txt"
        np.savetxt("sqr_instances/"+file, test_matrix, fmt="%.2f", delimiter=",")

        m=n
        test_matrix =generate_binary_instances(n, m, density=0.25)
        file=f"bin25_{test_matrix.shape[0]}_{it+1}.txt"
        np.savetxt("bin_instances/"+file, test_matrix, fmt="%d", delimiter=",")

        test_matrix =generate_binary_instances(n, m, density=0.5)
        file=f"bin50_{test_matrix.shape[0]}_{it+1}.txt"
        np.savetxt("bin_instances/"+file, test_matrix, fmt="%d", delimiter=",")

        test_matrix =generate_binary_instances(n, m, density=0.75)
        file=f"bin75_{test_matrix.shape[0]}_{it+1}.txt"
        np.savetxt("bin_instances/"+file, test_matrix, fmt="%d", delimiter=",")

for n in NN:
    for m in MM:
            if n>m:
                for it in range(5):
                    test_matrix = generate_nonsqr_instances(n,m)#np.round(np.random.rand(n, m),decimals=0)  # Example similarity/distance matrix
                    file=f"nonsqr_{n}_{m}_{it+1}.txt"
                    np.savetxt("nonsqr_instances/"+file, test_matrix, fmt="%f.2f", delimiter=",")

                    test_matrix =generate_binary_instances(n, m, density=0.25)
                    file=f"bin25_{n}_{m}_{it+1}.txt"
                    np.savetxt("bin_instances/"+file, test_matrix, fmt="%d", delimiter=",")

                    test_matrix =generate_binary_instances(n, m, density=0.5)
                    file=f"bin50_{n}_{m}_{it+1}.txt"
                    np.savetxt("bin_instances/"+file, test_matrix, fmt="%d", delimiter=",")

                    test_matrix =generate_binary_instances(n, m, density=0.75)
                    file=f"bin75_{n}_{m}_{it+1}.txt"
                    np.savetxt("bin_instances/"+file, test_matrix, fmt="%d", delimiter=",")