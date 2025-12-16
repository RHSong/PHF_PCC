import math
import numpy as np
"""
calculate Wigner (small) d-matrix, en.wikipedia.org/wiki/Wigner_D-matrix
d^j_nm = <jn|exp(-i beta Sy)|jm>
D^j_nm = exp(-in alpha) d^j_nm exp(-im gamma)
"""
def WignerMat(j,n,m,alpha,beta,gamma):
	prefac = 1
	prefac *= math.factorial(j+n) * math.factorial(j-n)
	prefac *= math.factorial(j+m) * math.factorial(j-m)
	prefac = math.sqrt(prefac)
	smax = min(j+m,j-n)
	smin = max(0,m-n)
	Sum = 0
	for s in range(smin,smax+1):
		p1 = n - m + s
		p2 = 2*j + m - n -2*s
		p3 = n - m + 2*s
		numer = (-1)**p1 * (math.cos(beta/2))**p2 * (math.sin(beta/2))**p3
		denom = 1
		denom *= math.factorial(j+m-s) * math.factorial(s)
		denom *= math.factorial(n-m+s) * math.factorial(j-n-s)
		Sum += numer / denom
	djnm = prefac * Sum
	Djnm = np.exp(-1j * n*alpha) * djnm * np.exp(-1j * m*gamma)
	return Djnm
