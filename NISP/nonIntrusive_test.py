
import numpy as np
import math
import itertools
import matplotlib.pyplot as plt

from generateSamples import normalSamples, gaussHermPts
from modesKL import modesKL
from truncatedKL import truncatedKL
from Galerkin1DElliptic import Galerkin1DElliptic
from hermPoly import pherm
from multiIndex import multiIndex
from solEllipticEq import solution 

# # some_file.py
# import sys
# sys.path.insert(0, '/Users/joaoreis/Documents/Codes/Python Notebooks/Galerkin Method/Routines')

# from massMatrix import galerkinMassMatrix
# from galerkinSourceVector import galerkinSource
# from galerkin1D import galerkinSolution


def nisp(Nel, f_source, alpha, beta, u_alpha, u_beta, Nkl, No, mu, sigma, l, M):


	## Generate M samples of $\underline{\xi}^{(j)}$ that are normal distributed. 
	#    We choose this distribution because we it's a nice one;
	# xi_samples = normalSamples(Nkl, M)
	xi_samples, weights_N, Nq = gaussHermPts(Nkl,No)
	print("Total of " + str(Nq) + " Gauss-Hermite points")





	## Compute the KL modes $\sqrt{\lambda_i}$ and eigenvectors $\phi_i(x)$ of $\hat{k}$;
	print("Computing KL modes")
	var = sigma**2
	Eign = modesKL(Nel, var, l) 
	print(np.sqrt(np.flip(Eign[0],axis=0)))



	## Compute the samples $\hat{k}^{(j)}$ using truncated KL-expansion;
	print("Computing truncated KL expansion")
	k_samples = truncatedKL(Nel, xi_samples, mu, Eign)


	## Solve $M$ deterministic problems $\partial_x [ \hat{k}^{(j)} \partial_x u(x, k(\theta))] = f(x)$;
	
	

	print("Computing solution samples")
	solution_samples = np.zeros((Nq,Nel))

	# print(".   building source")
	# F = galerkinSource(Nel, f_source, alpha, beta)



	print(".   building mass matrices and solutions")
	for q in range(Nq):
		if q % 1000 == 0:
			print(".         computing solution " + str(q))
		# A = galerkinMassMatrix(Nel, alpha, beta, k_samples[:,q])
		solution_samples[q,:] = Galerkin1DElliptic(Nel, f_source, alpha, beta, u_alpha, u_beta,\
		                                           k_samples[:,q])
		# solution_samples[q,:] = galerkinSolution(Nel, A, F, u_alpha, u_beta)


	###############################################################################
	## Generate the coefficients of a spectral expansion of the solution 
	##     by an orthogonal projection approach


	# multi-indix set to generate the spectral expansion
	print("Computing coefficients of spectral expansion")
	index = multiIndex(Nkl,No)
	


	# number of terms in the expansion
	P = int(math.factorial(No + Nkl)/(math.factorial(No)*math.factorial(Nkl))) - 1


	# creatres vector to store coefficients and matrix to store polynomials
	uCoeff = np.zeros((Nel, P+1))

	print(".     generating and evaluating " + str(P+1)+ " polynomials")
	psi_matrix = np.zeros((P+1, Nq))
	for s in range(P+1):
		for q in range(Nq):
			psi_matrix[s,q] = pherm(xi_samples[:,q], No, s, index)

	# for s in range(P+1):
	# 	for m in range(M):
	# 		p = p_N[m]
	# 		psi_matrix[s,m] = pherm(xi_samples[:,m], No, s, index)


	print(".     fill in " + str(P+1)+ " coefficients")
	# fill the coeffiients vector 
	for s in range(P+1):
		if s % 100 == 0:
			print(".           Computing coefficient: " + str(s))

		# Gauss-Hermite quadrature approximation of the integral
		integral = np.zeros(Nel)
		for q in range(Nq):
			w = np.prod(weights_N[q])
			integral += w*solution_samples[q,:]*psi_matrix[s,q]

		uCoeff[:,s] = integral

	###############################################################################


	## Computes the solution using spectral expansion

	print("Computing exact solution")
	## Computes exact solution using truncated KL expansion
	xi_normal = normalSamples(Nkl, M) # new samples for the error
	u0 = u_alpha
	uN = u_beta
	exactSolution_samples = np.zeros((Nel, M))
	for m in range(M):
		exactSolution_samples[:,m] = solution(f_source, Nel, alpha, beta, u0, uN, xi_normal[:,m], mu, sigma, l)


	#computing hermite polynomials for new samples
	print("Compute Hermite polynomials for new samples")
	newPsi_matrix = np.zeros((P+1, M))
	for s in range(P+1):
		for m in range(M):
			newPsi_matrix[s,m] = pherm(xi_normal[:,m], No, s, index)


	print("Computing PC solution for each normal distributed new sample")
	num_sol = np.zeros((Nel, M))

	for i in range(Nel):
		for m in range(M):
			num_sol[i,m] = sum(uCoeff[i,:]*newPsi_matrix[:,m])
	

	Omega = np.linspace(alpha, beta, Nel)
	
	image = plt.plot(Omega, num_sol[:,0], '-0', Omega, exactSolution_samples[:,0])
	plt.legend(['numerical solution', 'exact solution'])
	plt.savefig('image.png')




	print("Computing truncation error")
	## ERROR	
	error = 0
	ind = int(Nel/2) 
	for m in range(M):
		# exact_solution = solution(f_source, Nel, alpha, beta, u0, uN, sigma, k_samples[:,m])
		error += np.linalg.norm(num_sol[ind,m] - exactSolution_samples[ind,m])/M


	return num_sol, error
