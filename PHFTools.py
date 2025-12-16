import numpy as np
import sys
from parameter_fc import *
from solve_ci import *
from copy import deepcopy
Tol = 1e-6
Id = np.identity(NSO)
def EigenSolver(H,S):
	evals,evecs = np.linalg.eigh(S)
	Dim = len(evals)
	invS = [1/np.sqrt(abs(i)) for i in evals if (np.abs(i) > Tol)]
	dim = len(invS)
	Sinv = np.diag(invS)
	X = np.matmul(evecs[:,Dim-dim:],Sinv)
	h = X.T.conj() @ H @ X
	evals,evecs = np.linalg.eigh(h)
	evecs = X @ evecs
	return evals, evecs

def VecDecomp(vec):
	M = vec.reshape((ncisp*ncipg,ncik),order='F')
	fk = M[0,:] / M[0,0]
	M = M[:,0].reshape((ncisp,ncipg),order='F')
	fpg = M[0,:] / M[0,0]
	fsp = M[:,0]
	return fsp, fpg, fk

def EvalOvlp(r1,r2,rpg,rk,fsp,fpg,fk):
	Ovlp = 0
	newMOs = np.identity(NSO,dtype=complex)
	for il in range(ngrid[0]):
		for iy in range(ngrid[1]):
			for ipg in range(npg):
				for ik in range(ncik):
					R = r1[il,:,:] @ r2[iy,:,:] @ rpg[ipg,:,:] @ rk[ik,:,:]
					rho0, s0 = Rdm(R,newMOs[:,:NOccSO])
					s0 = weightsp[il,iy] * s0
					for ni in range(ncisp):
						for nj in range(ncisp):
							if (SP == 2):
								wig = WignerMat(J,ni-J,nj-J,roota[il],rootb[iy],rooty[il])
							else:
								wig = WignerMat(J,ni-J,nj-J,roota[il],rootb[il],rooty[iy])
							for p in range(ncipg):
								for q in range(ncipg):
									wpg = weightpg[ipg,p,q]
									w = wig * wpg * fsp[ni].conj() * fsp[nj] * fpg[p].conj() * fpg[q] * fk[ik]
									Ovlp += s0 * w
	return Ovlp


#def BuildHSG(HOne,HTwo,MOs,Z0):
#	Z = 1 * Z0
##	if (np.max(np.abs(Z.real))) > 10 or (np.max(np.abs(Z.imag))) > 10:
##		sys.exit("Z too Large")
#	if (SP == 2):
#		Z = makeUHF(Z)
#	if (CmplxConj == 0):
#		Z.imag = 0
#	Hmat = np.zeros([nci,nci], dtype=complex)
#	Smat = np.zeros([nci,nci], dtype=complex)
#	Gradmat = np.zeros([nci,nci,NSO,NSO], dtype=complex)
#	Rdmmat = np.zeros([nci,nci,NSO,NSO], dtype=complex)
#	newMOs = np.identity(NSO,dtype=complex)
#	newMOs[NOccSO:,:NOccSO] = Z
#	for il in range(ngrid[0]):
#		for iy in range(ngrid[1]):
#			for ipg in range(npg):
#				R = R1[il,:,:] @ R2[iy,:,:] @ Rpg[ipg,:,:]
#				R = ao2mo(R,MOs,2)
#				rho0, s0 = Rdm(R,newMOs[:,:NOccSO])
#				h0 = Kernels(HOne,HTwo,rho0,s0)
#				g0 = Gradient(HOne,HTwo,rho0,s0)
#				if (CmplxConj == 2):
#					rho0c, s0c = Rdm(R,newMOs[:,:NOccSO].conj())
#					g0c = Gradient(HOne,HTwo,rho0c,s0c)
#					rho1, s1 = RdmConj(R,newMOs[:,:NOccSO])
#					h1 = Kernels(HOne,HTwo,rho1,s1)
#					g1 = Gradient(HOne,HTwo,rho1,s1)
#					rho1c, s1c = RdmConj(R,newMOs[:,:NOccSO].conj())
#					g1c = Gradient(HOne,HTwo,rho1c,s1c)
#				else:
#					h1 = h0
#					g1 = g0
#					s1 = s0
#					g0c = g0
#					s0c = s0
#					g1c = g0
#					s1c = s0
#					rho0c = rho0
#					rho1 = rho0
#					rho1c = rho0
#				rho0 = rho0 * weightsp[il,iy] * s0
#				rho1 = rho1 * weightsp[il,iy] * s1
#				rho0c = rho0c * weightsp[il,iy] * s0c
#				rho1c = rho1c * weightsp[il,iy] * s1c
#				h0 = weightsp[il,iy] * h0
#				s0 = weightsp[il,iy] * s0
#				h1 = weightsp[il,iy] * h1
#				s1 = weightsp[il,iy] * s1
#				g0 = g0 * weightsp[il,iy]
#				g1 = g1 * weightsp[il,iy]
#				g0c = g0c * weightsp[il,iy]
#				s0c = s0c * weightsp[il,iy]
#				g1c = g1c * weightsp[il,iy]
#				s1c = s1c * weightsp[il,iy]
#				for ni in range(ncisp):
#					for nj in range(ncisp):
#						if (SP == 2):
#							wig = WignerMat(J,ni-J,nj-J,roota[il],rootb[iy],rooty[il])
#						else:
#							wig = WignerMat(J,ni-J,nj-J,roota[il],rootb[il],rooty[iy])
#						for p in range(ncipg):
#							for q in range(ncipg):
#								wpg = weightpg[ipg,p,q]
#								m = ni + p * ncisp
#								n = nj + q * ncisp
#								w = wig * wpg
#								Hmat[m,n] += h0 * w
#								Smat[m,n] += s0 * w
#								Gradmat[m,n,:,:] += g0 * w
#								Rdmmat[m,n,:,:] += rho0 * w
#								if (CmplxConj == 2):
#									d = ncisp * ncipg
#									Hmat[m,n+d] += h1 * w
#									Smat[m,n+d] += s1 * w
#									Gradmat[m,n+d,:,:] += g1 * w
#									Rdmmat[m,n+d,:,:] += rho1 * w
#									Gradmat[m+d,n,:,:] += g1c * w
#									Rdmmat[m+d,n,:,:] += rho1c * w
#									Gradmat[m+d,n+d,:,:] += g0c * w
#									Rdmmat[m+d,n+d,:,:] += rho0c * w
#	norm = np.max(np.abs(Smat.real))
#	Hmat = Hmat / norm
#	Smat = Smat / norm
#	Gradmat = Gradmat / norm
#	Rdmmat = Rdmmat / norm
#	if (CmplxConj == 2):
#		Hmat[ncisp*ncipg:,:ncisp*ncipg] = Hmat[:ncisp*ncipg,ncisp*ncipg:].T.conj()
#		Hmat[ncisp*ncipg:,ncisp*ncipg:] = Hmat[:ncisp*ncipg,:ncisp*ncipg]
#		Smat[ncisp*ncipg:,:ncisp*ncipg] = Smat[:ncisp*ncipg,ncisp*ncipg:].T.conj()
#		Smat[ncisp*ncipg:,ncisp*ncipg:] = Smat[:ncisp*ncipg,:ncisp*ncipg]
#	if (np.max(np.abs(Hmat.T.conj()-Hmat)) > Tol or np.max(np.abs(Smat.T.conj()-Smat)) > Tol):
#		print(Hmat)
#		print(Smat)
#		sys.exit("Not Hermitian")
#	Hmat = 0.5 * (Hmat + Hmat.T.conj())
#	Smat = 0.5 * (Smat + Smat.T.conj())
#	return Hmat, Smat, Gradmat, Rdmmat

def LocalG(Gradmat,Rdmmat,Smat,E0,fsp,fpg,fk):
	G = 0
	Ovlp = 0
	Rdm = 0
	for ni in range(ncisp):
		for nj in range(ncisp):
			for p in range(ncipg):
				for q in range(ncipg):
					m = ni + p * ncisp
					n = nj + q * ncisp
					G += fsp[ni].conj() * fsp[nj] * fpg[p].conj() * fpg[q] * fk[0].conj() * fk[0] * Gradmat[m,n,:,:]
					Ovlp += fsp[ni].conj() * fsp[nj] * fpg[p].conj() * fpg[q] * fk[0].conj() * fk[0] * Smat[m,n]
					Rdm += fsp[ni].conj() * fsp[nj] * fpg[p].conj() * fpg[q] * fk[0].conj() * fk[0] * Rdmmat[m,n,:,:]
					if (CmplxConj == 2):
						d = ncisp * ncipg
						G += fsp[ni].conj() * fsp[nj] * fpg[p].conj() * fpg[q] * fk[1].conj() * fk[1] * Gradmat[m+d,n+d,:,:].conj()
						Ovlp += fsp[ni].conj() * fsp[nj] * fpg[p].conj() * fpg[q] * fk[1].conj() * fk[1] * Smat[m+d,n+d]
						Rdm += fsp[ni].conj() * fsp[nj] * fpg[p].conj() * fpg[q] * fk[1].conj() * fk[1] * Rdmmat[m+d,n+d,:,:].conj()
						G += fsp[ni].conj() * fsp[nj] * fpg[p].conj() * fpg[q] * fk[0].conj() * fk[1] * Gradmat[m,n+d,:,:]
						Ovlp += fsp[ni].conj() * fsp[nj] * fpg[p].conj() * fpg[q] * fk[0].conj() * fk[1] * Smat[m,n+d]
						Rdm += fsp[ni].conj() * fsp[nj] * fpg[p].conj() * fpg[q] * fk[0].conj() * fk[1] * Rdmmat[m,n+d,:,:]
						G += fsp[ni].conj() * fsp[nj] * fpg[p].conj() * fpg[q] * fk[0].conj() * fk[1] * Gradmat[m+d,n,:,:].conj()
						Ovlp += fsp[ni].conj() * fsp[nj] * fpg[p].conj() * fpg[q] * fk[1].conj() * fk[0] * Smat[m+d,n]
						Rdm += fsp[ni].conj() * fsp[nj] * fpg[p].conj() * fpg[q] * fk[0].conj() * fk[1] * Rdmmat[m+d,n,:,:].conj()
	G -= E0 * Rdm
	return G[NOccSO:,:NOccSO] / Ovlp

def Rdm(R,MOs):
	M = np.matmul(MOs.T.conj(), R)
	M = np.matmul(M, MOs)
	Minv = np.linalg.inv(M)
	det = np.linalg.det(M)
	rdm = np.matmul(R, MOs)
	rdm = np.matmul(rdm, Minv)
	rdm = np.matmul(rdm, MOs.T.conj())
	return rdm, det

def RdmConj(R,MOs):
	M = np.matmul(MOs.T.conj(), R)
	M = np.matmul(M, MOs.conj())
	Minv = np.linalg.inv(M)
	det = np.linalg.det(M)
	rdm = np.matmul(R, MOs.conj())
	rdm = np.matmul(rdm, Minv)
	rdm = np.matmul(rdm, MOs.T.conj())
	return rdm, det

def Kernels(HOne,HTwo,rdm,det):
	h = np.einsum('ij,ji',HOne,rdm)
	h += 0.5*np.einsum('ijkl,ki,lj',HTwo,rdm,rdm)
	h = det * h
	return h

def Gradient(HOne,HTwo,rho,det):
	h = Makeh(HOne,HTwo,rho)
	gamma = MakeGamma(HTwo,rho)
	g = h * rho
	g += (Id-rho) @ (HOne + gamma) @ rho
	g = g * det
	return g

def MakeGamma(HTwo,rho):
	return np.einsum('ijkl,lj->ik',HTwo,rho)

def Makeh(HOne,HTwo,rho):
	h = np.trace(HOne @ rho)
	gamma = MakeGamma(HTwo,rho)
	h += 0.5 * np.trace(gamma @ rho)
	return h

def MLdecompose(Z):
	Id = np.identity(NOccSO)
	mat = Id + Z.T @ Z.conj()
	L = np.linalg.cholesky(mat)
	Id = np.identity(NSO-NOccSO)
	mat = Id + Z.conj() @ Z.T
	M = np.linalg.cholesky(mat)
	return L, M

def SortOrb(MO,NAO,NOccAO,Type):
# if type = 1
# rearrange orbitals from OA VA  0  0 to OA  0 VA  0
#                          0  0 OB VB     0 OB  0 VB
# if type = 0, do the inverse
	if (Type == 1):
		NVrtAO = NAO-NOccAO
		newMO = np.zeros([NSO,NSO])
		newMO[:,:NOccAO] = MO[:,:NOccAO]
		newMO[:,NOccAO:NOccSO] = MO[:,NAO:NAO+NOccAO]
		newMO[:,NOccSO:NOccSO+NVrtAO] = MO[:,NOccAO:NOccAO+NVrtAO]
		newMO[:,NAO+NOccAO:NSO] = MO[:,NAO+NOccAO:NSO]
		return newMO
	elif (Type == 0):
		newMO = np.zeros([NSO,NSO])
		newMO[:,:NOccAO] = MO[:,:NOccAO]
		newMO[:,NOccAO:NAO] = MO[:,NOccSO:NOccSO+NVrtA]
		newMO[:,NAO:NAO+NOccAO] = MO[:,NOccAO:NOccSO]
		newMO[:,NAO+NOccAO:NSO] = MO[:,NAO+NOccAO:NSO]
		return newMO

def makeUHF(Z):
# for parameter Zai, zero out the spin-flip part
	newZ = Z
	newZ[NVrtA:,:NOccA] = 0
	newZ[:NVrtA,NOccA:] = 0
	return newZ

def MFEne(HOne,HTwo,Z):
	newMOs = np.eye(NSO,NOccSO,dtype=complex)
	newMOs[NOccSO:,:NOccSO] = Z
	rdm, ovlp = Rdm(np.eye(NSO),newMOs)
	h = Kernels(HOne,HTwo,rdm,ovlp)
	return h / ovlp

def MFGrad(HOne,HTwo,Z, E0):
	Id = np.identity(NSO)
	newMOs = np.eye(NSO,NOccSO,dtype=complex)
	newMOs[NOccSO:,:NOccSO] = Z
	rdm, ovlp = Rdm(np.eye(NSO),newMOs)
	gk = Gradient(HOne,HTwo,rdm,ovlp) - E0 * rdm * ovlp
	return gk[NOccSO:,:NOccSO] / ovlp

def GetThouless(MO):
	M = MO[:NOccSO,:NOccSO]
	ovlp = np.linalg.det(M)
	Minv = np.linalg.inv(M)
	Z = MO[NOccSO:,:NOccSO] @ Minv
	return ovlp, Z

def Thouless2MOs(Z):
	MOs = np.identity(NSO, dtype=complex)
	MOs[NOccSO:,:NOccSO] = Z
	MOs[:NOccSO,NOccSO:] = -Z.T.conj()
	L, M = MLdecompose(Z)
	Linv = np.linalg.inv(L)
	Minv = np.linalg.inv(M)
	newMOs = np.zeros([NSO,NSO], dtype=complex)
	newMOs[:,:NOccSO] = MOs[:,:NOccSO] @ Linv.T
	newMOs[:,NOccSO:] = MOs[:,NOccSO:] @ Minv.T
	return newMOs

def AddFrozenMO(MO,NFO):
	newMO = np.zeros([NSO+NFO*2,NSO+NFO*2])
	newMO[:NFO,:NFO] = np.eye(NFO)
	newMO[NFO:NFO+NAO,NFO:NFO+NOccA] = MO[:NAO,:NOccA]
	newMO[NFO*2+NAO:,NFO:NFO+NOccA] = MO[NAO:,:NOccA]
	newMO[NAO+NFO:NAO+2*NFO,NOccA+NFO:NOccA+2*NFO] = np.eye(NFO)
	newMO[NFO:NFO+NAO,NOccA+2*NFO:NOccSO+2*NFO] = MO[:NAO,NOccA:NOccSO]
	newMO[NFO*2+NAO:,NOccA+2*NFO:NOccSO+2*NFO] = MO[NAO:,NOccA:NOccSO]
	newMO[NFO:NFO+NAO,NOccSO+2*NFO:] = MO[:NAO,NOccSO:]
	newMO[NFO*2+NAO:,NOccSO+2*NFO:] = MO[NAO:,NOccSO:]
	return newMO

def RmvFrozenMO(MO,NFO):
	newMO = np.zeros([NSO-NFO*2,NSO-NFO*2])
	newMO[:NAO-NFO,:NOccA-NFO] = MO[NFO:NAO,NFO:NOccA]
	newMO[NAO-NFO:,:NOccA-NFO] = MO[NFO+NAO:NSO,NFO:NOccA]
	newMO[:NAO-NFO,NOccA-NFO:NOccSO-2*NFO] = MO[NFO:NAO,NOccA+NFO:NOccSO]
	newMO[NAO-NFO:,NOccA-NFO:NOccSO-2*NFO] = MO[NFO+NAO:NSO,NOccA+NFO:NOccSO]
	newMO[:NAO-NFO,NOccSO-2*NFO:] = MO[NFO:NAO,NOccSO:]
	newMO[NAO-NFO:,NOccSO-2*NFO:] = MO[NFO+NAO:NSO,NOccSO:]
	return newMO
