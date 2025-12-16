import numpy as np
from math import *
from scipy import linalg
def ao2mo(mat,MOs,dim):
# transfer integral from ao to mo
	if (dim==2):
		newmat=np.einsum('ia,ab,bj -> ij',MOs.T.conj(),mat,MOs,optimize='optimal')
	if (dim==4):
		newmat=np.einsum('ia,jb,abcd,ck,dl -> ijkl',MOs.T.conj(),MOs.T.conj(),mat,MOs,MOs,optimize='optimal')
	return newmat

def getTrans(Ovlp):
	tol = 1e-10
	evals,evecs = linalg.eigh(Ovlp,lower=False)
	invS = [np.sqrt(abs(i))/(i+tol) for i in evals]
	Sinv = np.diag(invS)
	X = np.matmul(evecs,Sinv)
	return X
