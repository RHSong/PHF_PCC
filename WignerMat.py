import math
import numpy as np
"""
calculate Wigner (small) d-matrix, en.wikipedia.org/wiki/Wigner_D-matrix
d^j_nm = <jn|exp(-i beta Sy)|jm>
D^j_nm = exp(-in alpha) d^j_nm exp(-im gamma)
j, n, m may be integers or half-integers (n, m = -j, ..., j); factorials via the Gamma function.
"""
def _fact(x):
	return math.gamma(round(2 * x) / 2 + 1)
def WignerMat(j,n,m,alpha,beta,gamma):
	prefac = _fact(j+n) * _fact(j-n) * _fact(j+m) * _fact(j-m)
	prefac = math.sqrt(prefac)
	smax = int(round(min(j+m,j-n)))
	smin = int(round(max(0,m-n)))
	Sum = 0
	for s in range(smin,smax+1):
		p1 = int(round(n - m + s))
		p2 = int(round(2*j + m - n -2*s))
		p3 = int(round(n - m + 2*s))
		numer = (-1)**p1 * (math.cos(beta/2))**p2 * (math.sin(beta/2))**p3
		denom = _fact(j+m-s) * _fact(s) * _fact(n-m+s) * _fact(j-n-s)
		Sum += numer / denom
	djnm = prefac * Sum
	Djnm = np.exp(-1j * n*alpha) * djnm * np.exp(-1j * m*gamma)
	return Djnm
