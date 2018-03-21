import numpy as np
import math

def pherm1D(No, x0):


	p0 = 1
	p=p0

	p1 = x0

	if No == 1: p = p1

	for i in range(2,No+1):
		p = p1*x0 - (i-1)*p0
		p0 = p1
		p1 = p


	p *= 1/np.sqrt(math.factorial(No))

	return p

