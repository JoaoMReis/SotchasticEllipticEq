import numpy as np
import itertools

# This routine generates M standard normal distributed samples 
#    with dimension Nkl



def normalSamples(Nkl, M):
	# Nkl   -> dimension of the sample vector
	# M     -> number of samples


	normal_samples = np.zeros((Nkl, M))
	for i in range(Nkl):
	    normal_samples[i, :] = np.random.normal(0, 1, M) 

	return normal_samples

def gaussHermPts(Nkl,No):
	# Nkl   -> dimension of the sample vector
	# M     -> number of samples


	pw = np.polynomial.hermite.hermgauss(No)
	points = pw[0]*np.sqrt(2)
	weights = pw[1]/np.sqrt(np.pi)

	points_N = list(itertools.product(points, repeat = Nkl))
	weights_N = list(itertools.product(weights, repeat = Nkl))
	Nq = len(weights_N)

	gaussHermPts = np.zeros((Nkl, Nq))
	for q in range(Nq):
		gaussHermPts[:,q] = points_N[q][:]
	

	return gaussHermPts, weights_N, Nq









