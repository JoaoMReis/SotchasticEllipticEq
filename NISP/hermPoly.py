# This routine computes the multivariate Hermite polynomials 
#   Psi_k as in (2.50) in  the book

import numpy as np 
import scipy.special
import math
from fodase import pherm1D


# This routine returns the nomalised multivariate Hermite polynomials 
#    with weight function exp(-x^2/2)
def pherm(x0, No, k, index):
	#  INPUTS: x0    -> point to be evaluated
	#          No    -> order of spectral expansion
	#          k     -> index position
	#          index -> index that gives the correct order/point evaluations

	
	N = len(x0)
	gamma = index[k]
	
	# order_vector = np.ones(No+1)*range(No+1)
	# psi_vector = list()
	# for j in range(No+1):
	# 	psi_vector.append(scipy.special.hermite(order_vector[j]))
	
	# order_vector = np.ones(No+1)*range(No+1)
	# psi_vector = np.zeros((No+1,N))
	# for j in range(No+1):
	# 	psi_vector[j,:] = scipy.special.eval_hermitenorm(j,x0)

	psi = 1
	# psi_normalised = 1
	
	if k != 0:
		for i in range(N):
			order = gamma[i]
			xi = x0[i]
			# need to divide by 1/sqrt(2) to get the right type of polynomial
			# psi *= scipy.special.eval_hermitenorm(order,xi)
			psi *= pherm1D(order, xi)
			# psi *= scipy.special.hermitenorm(order)(xi)			
			# # normalisation
			# psi_normalised *=(2**order + math.factorial(order))  

		# psi = psi/psi_normalised

	return psi