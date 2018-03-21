import numpy as np

# This routine computes the truncated KL expansion of a set
#  of M samples of a field K given a set of normal distributed samples

def truncatedKL(Nel, xi_samples, mu_vector, Eign):

	# INPUTS: Nel        -> Number of points in space
	#         xi_samples -> Normal distributed N(0,1) samples
	#         mu         -> Mean of the field K
	#         Eign       -> Matrix with KL modes and eignvectors of field K

	# OUTPUT: k_samples  -> (Nel, M) matrix with M samples of a field defined 
	#                          in Nel space points



	Nkl = len(xi_samples[:,0])
	M = len(xi_samples[0,:])
	k_samples = np.zeros( (Nel, M) ) # k(x, xi, j), j = 1,..., M

	# here we keeping the lardest Nkl modes corresponding eigenvectors 
	# modes = np.sqrt(Eign[0])[Nel - Nkl: Nel]
	modes = np.sqrt(np.flip(Eign[0],axis=0)[0: Nkl])
	# eigenvectors[j,i] gives the eigenvector on point xi_j for mode i 
	eigenvectors = np.flip(Eign[1],axis=1)[0: Nkl,:]
	# print("engenvectors" + str(eigenvectors))
	

	# mu_vector = mu*np.ones(Nel)
	for j in range(M):
	    prod_aux = np.zeros( (Nkl,Nel) )
	    for i in range(Nkl):
	        prod_aux[i,:] = modes[i]*xi_samples[i,j]*eigenvectors[i,:]
	    k_samples[:,j] = np.exp(mu_vector[:] + np.sum(prod_aux))

	return k_samples