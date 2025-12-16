import numpy as np
from parameter_fc import *
from ao2mo import *
from itertools import product
def BuildSz(NAO):
	NSO = 2 * NAO
	Sz = np.zeros([NSO,NSO])
	for i in range(NAO):
		Sz[i,i] = 1
		Sz[i+NAO,i+NAO] = -1
	return Sz / 2

def BuildSy(NAO):
	NSO = 2 * NAO
	Sy = np.zeros([NSO,NSO], dtype=complex)
	for i in range(NAO):
		Sy[i,i+NAO] = -1j
		Sy[i+NAO,i] = 1j
	return Sy / 2

def BuildSx(NAO):
	NSO = 2 * NAO
	Sx = np.zeros([NSO,NSO])
	for i in range(NAO):
		Sx[i,i+NAO] = 1
		Sx[i+NAO,i] = 1
	return Sx / 2

def BuildS2(NAO):
	NSO = NAO*2
	S2One = np.zeros([NSO,NSO])
	S2Two = np.zeros([NSO,NSO,NSO,NSO])
	for i in range(NAO):
		S2One[i,i] = 0.75
		S2One[i+NAO,i+NAO] = 0.75
	for i in range(NAO):
		for j in range(NAO):
			S2Two[i+NAO,j+NAO,i+NAO,j+NAO] = 0.25
			S2Two[i,j,i,j] = 0.25
			S2Two[i,j+NAO,i+NAO,j] = 1
			S2Two[i,j+NAO,i,j+NAO] = -0.5
	S2Two = antisymm(S2Two)
	return S2One, S2Two

def calcS(Ref,NOcc,NAO):
	NSO = 2 * NAO
	Sx = BuildSx(NAO)
	Sy = BuildSz(NAO)
	Sz = BuildSz(NAO)
	rdm = np.matmul(Ref[:,:NOcc],Ref[:,:NOcc].T.conj())
	ExpSx = np.einsum('ij,ji',Sx,rdm)
	ExpSy = np.einsum('ij,ji',Sy,rdm)
	ExpSz = np.einsum('ij,ji',Sz,rdm)
	return ExpSx, ExpSy, ExpSz

def calcS2(Ref,NOcc,NAO):
	NSO = NAO*2
	S2One, S2Two = BuildS2(NAO)
	rdm = np.matmul(Ref[:,:NOcc],Ref[:,:NOcc].T.conj())
	ExpS2 = np.einsum('ij,ji',S2One,rdm)
	ExpS2 += np.einsum('ijkl,ki,lj',S2Two,rdm,rdm)
	ExpS2 += -np.einsum('ijkl,li,kj',S2Two,rdm,rdm)
	return ExpS2

def SpinDen(MO):
	mo = orth @ MO
	mo1 = mo[:NAO,:NAO]
	mo2 = mo[NAO:,NAO:]
	Sa = mo1[:,:NOccAO] @ mo1[:,:NOccAO].T.conj()
	Sb = mo2[:,:NOccAO] @ mo2[:,:NOccAO].T.conj()
	return np.diag(Sa-Sb)

def SSHG(MO,NAO,NOccSO):
	Sx = BuildSx(NAO)
	Sy = BuildSy(NAO)
	Sz = BuildSz(NAO)
	Ox = MO.T.conj() @ Sx @ MO
	Oy = MO.T.conj() @ Sy @ MO
	Oz = MO.T.conj() @ Sz @ MO
	Olist = [Ox, Oy, Oz]
	A = np.eye(3, dtype=complex) * NOccSO / 4
	for i in range(3):
		for j in range(3):
			A[i,j] -= np.trace(Olist[i] @ Olist[j])
	return A

#def SGHFS2(Ref,f,NOcc,NAO):
#	NSO = NAO*2
#	Sz = BuildSz(NAO)
#	S2One, S2Two = BuildS2(NAO)
#	Sz = ao2mo(Sz,Ref,2)
#	S2One = ao2mo(S2One,Ref,2)
#	S2Two = ao2mo(S2Two,Ref,4)
#	Ovlp = 0
#	S1 = 0
#	S2 = 0
#	MOs = np.zeros([NSO,NOcc])
#	for i in range(NOcc):
#		MOs[i,i] = 1
#	for l,y in product(range(ngrid[0]),range(ngrid[1])):
#		R = BuildRotMat(l,l,y)
#		R = ao2mo(R,Ref,2)
#		M = MOs.T.conj() @ R @ MOs
#		Minv = np.linalg.inv(M)
#		det = np.linalg.det(M)
#		rdm = R @ MOs @ Minv @ MOs.T.conj()
#		ExpS1 = np.einsum('ij,ji',Sz,rdm)
#		ExpS2 = np.einsum('ij,ji',S2One,rdm)
#		ExpS2 += np.einsum('ijkl,ki,lj',S2Two,rdm,rdm)
#		ExpS2 += -np.einsum('ijkl,li,kj',S2Two,rdm,rdm)
#		ExpS1 = ExpS1 * weights[l,y] * det
#		ExpS2 = ExpS2 * weights[l,y] * det
#		ovlp = weights[l,y] * det
#		for ni in range(-j,j+1):
#			for nj in range(-j,j+1):
#				wig = WignerMat(j,ni,nj,roota[l],rootb[l],rooty[y])
#				S1 += wig * f[ni+j] * f[nj+j] * ExpS1
#				S2 += wig * f[ni+j] * f[nj+j] * ExpS2
#				Ovlp += wig * f[ni+j] * f[nj+j] * ovlp
#	return S1/Ovlp, S2/Ovlp

