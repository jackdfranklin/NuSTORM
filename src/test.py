import numpy as np

def chi_squared(N_exp, N_theo):
        χ_squared = np.empty(10,10)

        for i in range(10):
            for j in range(10):
                χ_squared[i, j] = (N_exp[i, j] - N_theo[i, j])**2 /(N_exp[i,j])

        χ_squared = χ_squared.flatten()

        return χ_squared

chi_squared(4, 4)