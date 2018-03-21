
import numpy as np

from generateSamples import normalSamples, gaussHermPts
from truncatedKL import truncatedKL
from modesKL import modesKL

# This routine computes the exact solution of the 
#    1D stochastic elliptic equation

def solution(f, Nel, alpha, beta, u0, uN, xi, mu, sigma, l):
	# INPUTS: f     -> source vector
	#         u0    -> LHS boundary condition
	#         uN    -> RHS boundary condition
	#         mu    -> mean of k field
	#         sigma -> standard deviation of k field

	#OUTPUT:  exact solution foda-se!


	Nkl = len(xi)

	var = sigma**2
	Eign = modesKL(Nel, var, l) 
	xi = np.reshape(xi,(Nkl,1))
	k_field = truncatedKL(Nel, xi, mu, Eign)

	Omega = np.linspace(alpha, beta, Nel)    # spacial domain

	# vectorise source vector and boundary conditions
	f_vector = np.zeros(Nel)
	for i in range(Nel):
		f_vector[i] = f(Omega[i])
	# u0_vector = u0*np.ones(Nel)
	# uN_vector = uN*np.ones(Nel)

	# #vectorise means
	# mean_uN = np.mean(uN_vector)
	# mean_u0 = np.mean(u0_vector)

	#exact solution
	exact_solution =  np.zeros(Nel)
	for i in range(Nel):
		exact_solution[i] = -(f_vector[i]*Omega[i]**2)/(2*k_field[i]) + \
		(uN + f_vector[i]/(2*k_field[i]) - u0)*Omega[i] + u0



	return exact_solution