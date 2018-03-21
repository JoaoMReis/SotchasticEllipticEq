import numpy as np
from scipy.linalg import eigh

def modesKL(Nel, var, l):
    # INPUTS: Nel -> number of elements
    # .       var -> variance of the parameter
    # .       l   -> correlation lenght of the parameter
    
    # OUPUTS: matrix with \lambda_i and \phi_i of the KL
    #             expansion K = \sum_i \lambda_i \phi_i(x) \eta_i(\theta)

    # START
    
    # number of finite elements
    Omega = np.linspace(0,1,Nel)
    h = Omega[1] - Omega[0]

    # Define the Kernel
    K = np.zeros((Nel, Nel))
    M = np.zeros((Nel, Nel))
    for i in range(Nel):
        # xi = Omega[i]
        xi = h*(i + 0.5)
        for j in range(i,Nel):
            xj = h*(j + 0.5)
            xij = xi - xj
            K[i,j] = var*np.exp(-xij*xij/(2*l**2))*h**2
            if i != j:
                K[j,i] = K[i,j]
        M[i,i] = h
    
    Eign = eigh(K,M,eigvals_only=False)
    
    return Eign