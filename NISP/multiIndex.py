import itertools
import numpy as np



def multichoose(n,k):
    if k < 0 or n < 0: return "Error"
    if not k: return [[0]*n]
    if not n: return []
    if n == 1: return [[k]]

    return [[0]+val for val in multichoose(n-1,k)] + \
        [[val[0]+1]+val[1:] for val in multichoose(n,k-1)]

def multiIndex(N, No):

    multi_index = []
    for n in range(No+1):
        index = multichoose(N,n)
        multi_index.extend(index)

    return multi_index

# def multiIndex(N, No):

# 	# This routine returns the set al multi_indeces \lambda(p) as in
# 	#    2.48 page 32 of the book

# 	# INPUTS: N -> dimension of xi
#     #         No -> order of the polynomial

#     # OUTPUT: the set of multi-indices lambda(p) where each set
#     #             corresponds to an order up to p sets. See reference 
#     #             for more details

	

# 	S = np.ones(No+1, dtype='int')*range(No+1) # all orders that we need
# 	# combinat = list(itertools.combinations_with_replacement(S,N))
# 	# # S = np.flip(S,0)

# 	# permutations1 = list(itertools.product(S,repeat = N))

# 	# permutations of the orders in the number of variables
# 	permutations = []
# 	permutations.append(list(itertools.product(S,repeat = N)))
	      
# 	# # multi_index is a vector with elements [[lambda(1)], lambda(2), ..., lambda(No)]

# 	multi_index = []

# 	for n in range(No+1):
# 	    for l in range(len(permutations)):
# 	        myset = permutations[l]
# 	        for i in range(len(myset)):
# 	            mytuple = myset[i]
# 	            if sum(mytuple) == n:
# 	                multi_index.append(mytuple)

# 	# combinat = list(itertools.product(S,repeat = N))

	
# 	# for gamma in combinat:
# 	# 	if np.sum(gamma) > No:
# 	# 		combinat.remove(gamma)
	

# 	return multi_index
