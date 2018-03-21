from nonIntrusive_test import nisp
import numpy as np
import itertools

from IPython.display import Markdown, display
import matplotlib.pyplot as plt


import math
import scipy.special


from generateSamples import normalSamples
from modesKL import modesKL
from truncatedKL import truncatedKL
from Galerkin1DElliptic import Galerkin1DElliptic
from hermPoly import pherm
from multiIndex import multiIndex

Nel = 10
f_source = lambda x: x
alpha = 0
beta = 1
u_alpha = 0
u_beta = 0
Nkl = 1
No = 10
mu = np.ones(Nel)
sigma = 0.1 
l = 0.1
M = 100000

# index = multiIndex(Nkl, No)

# print(index)
solution, error = nisp(Nel, f_source, alpha, beta, u_alpha, u_beta, Nkl, No, mu, sigma, l, M)

Omega = np.linspace(alpha, beta, Nel)

# plt.plot(Omega, solution, '-o', Omega, exact_solution,)
# plt.savefig("solution with NISP.png")
print(error)


