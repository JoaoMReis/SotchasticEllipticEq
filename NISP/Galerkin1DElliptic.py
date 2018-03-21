import numpy as np
from sympy import *
from scipy.interpolate import interp1d
from scipy import integrate





def Galerkin1DElliptic(N, f, alpha, beta, u_alpha, u_beta, k_field):
    
    
    # N = 10             number of spacial points
    # f = lambda x: -1   source function
    # u_alpha            LHS boundary value
    # u_beta             RHS boundary value
    # k_field            convective field term
     
    
    # problem set up
    Omega = np.linspace(alpha, beta, N)    # spacial domain
    u0 = u_alpha                           # for convineance
    uN = u_beta                            # for convineance
    dx = Omega[1] - Omega[0]               # step size
    
    # interpolating conductive term over Omega
    kfield_interp = interp1d(Omega, k_field, kind='cubic') 
    
    # setting orthogonal basis functions
    x, x_km1, x_k = symbols('x x_{k-1} x_k')
    phi_km1 = lambda z, x_km1, x_k: (z - x_km1)/(x_k - x_km1)
    phi_k = lambda z, x_km1, x_k: (x_k - z)/(x_k - x_km1)


    # filling mass matrix
    A_aux = np.zeros((N,N))
    for k in range(N-1):
        a = Omega[k]
        b = Omega[k+1]
        M_k = np.zeros((N,N))
        M_k[k,k] = integrate.quad(kfield_interp, a, b)[0]
        M_k[k,k+1] = -M_k[k,k]
        M_k[k+1,k] = -M_k[k,k]
        M_k[k+1,k+1] = M_k[k,k]
        A_aux += 1/(dx)**2*M_k
    # ------------------------------


    # filling source vector
    F_aux = np.zeros(N)
    for k in range(N-1):
        a = Omega[k]
        b = Omega[k+1]
        F_k = np.zeros(N)
        
        int_km1 = lambda x: f(x)*phi_km1(x,a,b)
        F_k[k] = integrate.quad(int_km1, a, b)[0]
        
        int_k = lambda x: f(x)*phi_k(x,a,b)
        F_k[k+1] = integrate.quad(int_k, a, b)[0]
        
        F_aux += F_k
    # ------------------------------
    

    # dirichlet boundary conditions
    A = A_aux[1:-1,1:-1]
    F = F_aux[1:-1]
    F[0] += -A_aux[1,0]*u0
    F[-1] += -A_aux[-2, -1]*uN
    # ------------------------------



    # constructing solution
    u_aux = np.linalg.solve(A,F)
    u = np.zeros(N)
    u[0] = u0
    u[-1] = uN
    u[1:-1] = u_aux
    
    return u