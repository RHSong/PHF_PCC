import numpy as np
from math import *
from pyscf import gto, scf
from ao2mo import *
def BuildHamitonian(nx,ny,u,pbcx=True,pbcy=True):
	n = nx * ny
	T0 = np.zeros((n,n))
	U0 = np.zeros((n,n,n,n))
	for i in range(n):
		U0[i,i,i,i] = u
	for i in range(nx):
		for j in range(ny-1):
			ni = ny * i + j
			nj = ny * i + j + 1
			T0[ni,nj] = T0[nj,ni] = -1
	for i in range(nx-1):
		for j in range(ny):
			ni = ny * i + j
			nj = ny * (i+1) + j
			T0[ni,nj] = T0[nj,ni] = -1
	if (pbcy and ny > 1):
		for i in range(nx):
			ni = ny * i
			nj = ny * i + ny - 1
			T0[ni,nj] += -1
			T0[nj,ni] += -1
	if (pbcx and nx > 1):
		for j in range(ny):
			ni = j
			nj = ny * (nx-1) + j
			T0[ni,nj] += -1
			T0[nj,ni] += -1
	return T0, U0

def BuildHamitonianSO(NAO,u,t):
# build H for 1D Hubbard model in spin orbitals 
	nso = 2*NAO
	T0=np.zeros((nso,nso))
	U0=np.zeros((nso,nso,nso,nso))
	for i in range(NAO-1):
		T0[i,i+1] = T0[i+1,i] = -t
		T0[NAO+i,NAO+i+1] = T0[NAO+i+1,NAO+i] = -t
		U0[i,i,i+NAO,i+NAO] = U0[i+NAO,i+NAO,i,i] = u
	T0[NAO-1,0] = T0[0,NAO-1] = T0[nso-1,NAO] = T0[NAO,nso-1] =-t
	U0[NAO-1,NAO-1,nso-1,nso-1] = U0[nso-1,nso-1,NAO-1,NAO-1] = u
	return T0,U0

# here the 2ele part is (ij|kl), which is equivalant to <ik|jl>,
# <i!j> = <j|b>delta(ab)<a|i> = Cja (C_dag)ai = (CC_dag)ji = rho(ji)
# <ik|jl>i!k!lj = (ij|kl) rho(lk)i!j - (ij|kl) rho(jk)i!l

# def 2-ele integral
def Jint(dm,eleint):
#	J=np.einsum('ij,klji',dm,eleint,optimize='greedy')
	J=np.einsum('ijkl,lk',eleint,dm,optimize='greedy')
	return J

def Kint(dm,eleint):
#	K=np.einsum('ij,kijl',dm,eleint,optimize='greedy')
	K=np.einsum('ijkl,jk',eleint,dm,optimize='greedy')
	return K

def Mulliken2Dirac(H2):
# from Mulliken to Dirac notation and antisymmetrize it
	tmp = np.einsum('ijkl->ikjl',H2)
	U = tmp - np.einsum('ijkl->ijlk',tmp)
	return U

def antisymm(H2):
	tmp = H2 - np.einsum('ijkl->ijlk',H2)
	H2new = tmp - np.einsum('ijkl->jikl',tmp)
	return H2new / 4

def MakeFock(H1,H2,MO_Occ):
	rho = MO_Occ @ MO_Occ.T.conj()
	Fock = H1 + np.einsum('ijkl,lj->ik', H2, rho)
	return Fock

def MakeEne(H1,H2,MO_Occ):
	rho = MO_Occ @ MO_Occ.T.conj()
	Ene = np.einsum('ij,ji', H1, rho)
	Ene += 0.25 * np.einsum('ijkl,ki,lj', H2, rho, rho)
	Ene -= 0.25 * np.einsum('ijkl,kj,li', H2, rho, rho)
	return Ene

def SemiCanon(H1,H2,MO,NOcc):
	F = MakeFock(H1,H2,MO[:,:NOcc])
	F = MO.T.conj() @ F @ MO
	Foo = F[:NOcc,:NOcc]
	Fvv = F[NOcc:,NOcc:]
	e, U = np.linalg.eigh(Foo)
	mo_o = MO[:,:NOcc] @ U
	e, U = np.linalg.eigh(Fvv)
	mo_v = MO[:,NOcc:] @ U
	return np.hstack((mo_o, mo_v))
