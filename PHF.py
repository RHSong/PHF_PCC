import numpy as np
import os
import pickle
from Spin import *
from parameter_fc import *
from PHFTools import *
from FixGauge import *
from scipy.optimize import minimize, basinhopping
#from ipopt import minimize_ipopt
from PHFTools_F import phftools
Sx = BuildSx(NAO)
Sy = BuildSy(NAO)
Sz = BuildSz(NAO)
Zeros = np.zeros([NSO,NSO,NSO,NSO])
fixgauge = False
NFO = 8
def optPHF(H1, H2, MOs, comm):
	if (os.path.exists("hubsghf.p")):
		z0 = pickle.load(open( "hubsghf.p", "rb" ))
	else:
		z0 = np.zeros(2 * (NSO-NOccSO)*NOccSO)
	minimizer_kwargs = {"method":"BFGS", "jac":True, "args":(H1,H2,MOs,comm),"options":{'gtol': 1e-4}}
	res = basinhopping(EandG, z0, minimizer_kwargs=minimizer_kwargs, niter=0, T=0.2, stepsize=0.2)
#	res = minimize(EandG, z0, method='BFGS', args=(H1,H2,MOs,comm), jac=True, options={'gtol': 5e-5})
#	cons = ({'type':'eq','fun': lambda x: SxCons(x,MOs), 'jac': lambda x: SxConsGrad(x,MOs)})
##	{'type':'eq','fun': lambda x: SyCons(x,MOs), 'jac': lambda x: SyConsGrad(x,MOs)})
#	res = minimize_ipopt(EandG, z0, args=(H1,H2,MOs,comm), jac=True, constraints=cons, tol=5e-5)
	print(res)
	sol = res.x[:NVrtSO*NOccSO] + res.x[NVrtSO*NOccSO:] * 1j
	pickle.dump(res.x, open( "hubsghf.p", "wb" ))
	z = sol.reshape([NSO-NOccSO,NOccSO])
	if (NFO > 0):
		z = FrozenProj(z)
	if (CmplxConj == 0):
		z.imag = 0
	if (SP == 2):
		z = SzProj(z)
	if (is_RHF):
		z = RHFProj(z)
	G2S = Thouless2MOs(z)
	newMOs = MOs @ G2S
	if (fixgauge):
		newMOs = FixGauge(newMOs)
	Hmat, Smat, Gradmat, Rdmmat = phftools.buildhsg(H1,H2,MOs,z,ngrid,nci,CmplxConj,ncisp,NOccSO, \
	R1,R2,Rpg,weightsp,weightpg,roota,rootb,rooty,J,SP,comm,npg,ncipg,NSO)
	evals, evecs = EigenSolver(Hmat,Smat)
	E0 = evals[0]
	fsp, fpg, fk = VecDecomp(evecs[:,0])
	print("CI Ene=", evals+Enuc)
	if (CmplxConj == 2):
		newMOs = FixCmplxGauge(fk[1],newMOs)
		fk[1] = 1.0
	return E0, newMOs, fsp, fpg, fk

def EandG(Z, *args):
	lamb = 0
	pickle.dump(Z, open( "newsol.p", "wb" ))
	z0 = Z[:NVrtSO*NOccSO] + Z[NVrtSO*NOccSO:] * 1j
	sol = z0.reshape(NSO-NOccSO,NOccSO)
	if (NFO > 0):
		sol = FrozenProj(sol)
	if (CmplxConj == 0):
		sol.imag = 0
	if (SP == 2):
		sol = SzProj(sol)
	if (is_RHF):
		sol = RHFProj(sol)
	HOne,HTwo,MOs, comm = args
#	Hmat, Smat, Gradmat, Rdmmat = BuildHSG(HOne,HTwo,MOs,sol)
	Hmat, Smat, Gradmat, Rdmmat = phftools.buildhsg(HOne,HTwo,MOs,sol,ngrid,nci,CmplxConj,ncisp,NOccSO, \
	R1,R2,Rpg,weightsp,weightpg,roota,rootb,rooty,J,SP,comm,npg,ncipg,NSO)
	evals, evecs = EigenSolver(Hmat,Smat)
	E0 = evals[0]
	fsp, fpg, fk = VecDecomp(evecs[:,0])
#	G = LocalG(Gradmat,Rdmmat,Smat,E0,fsp,fpg,fk)
	G = phftools.localg(Gradmat,Rdmmat,Smat,E0,fsp,fpg,fk,CmplxConj,NSO,nci,ncisp,ncipg)
	E0 = evals[0] + lamb * z0.T.conj() @ z0
	G = G[NOccSO:,:NOccSO] + lamb * z0.reshape(NSO-NOccSO,NOccSO)
	if (NFO > 0):
		G = FrozenProj(G)
	if (CmplxConj == 0):
		G.imag = 0
	if (SP == 2):
		G = SzProj(G)
	if (is_RHF):
		G = RHFProj(G)
	G = G.reshape((NSO-NOccSO)*NOccSO)
	output = open('PHFout', 'a')
	print("E and max |G|", evals[0]+Enuc, np.max(np.abs(np.imag(np.diag(Hmat/Smat)))),np.max(np.abs(G.real)), np.max(np.abs(G.imag)), file=output)
	output.close()
	return E0.real, np.concatenate((G.real,G.imag))

def SzCons(Z,MO):
	z0 = Z[:NVrtSO*NOccSO] + Z[NVrtSO*NOccSO:] * 1j
	z0 = z0.reshape(NSO-NOccSO,NOccSO)
	if (CmplxConj == 0):
		z0.imag = 0
	sz = ao2mo(Sz,MO,2)
	s_z = MFEne(sz,Zeros,z0)
	return s_z.real

def SzConsGrad(Z,MO):
	z0 = Z[:NVrtSO*NOccSO] + Z[NVrtSO*NOccSO:] * 1j
	z0 = z0.reshape(NSO-NOccSO,NOccSO)
	if (CmplxConj == 0):
		z0.imag = 0
	sz = ao2mo(Sz,MO,2)
	s_z = MFEne(sz,Zeros,z0)
	szg = MFGrad(sz,Zeros,z0,s_z)
	print("Sz and max |SG|", s_z.real, np.max(np.abs(szg.real)), np.max(np.abs(szg.imag)))
	return np.concatenate((szg.real,szg.imag))

def SxCons(Z,MO):
	z0 = Z[:NVrtSO*NOccSO] + Z[NVrtSO*NOccSO:] * 1j
	z0 = z0.reshape(NSO-NOccSO,NOccSO)
	if (CmplxConj == 0):
		z0.imag = 0
	sx = ao2mo(Sx,MO,2)
	s_x = MFEne(sx,Zeros,z0)
	return s_x.real

def SxConsGrad(Z,MO):
	z0 = Z[:NVrtSO*NOccSO] + Z[NVrtSO*NOccSO:] * 1j
	z0 = z0.reshape(NSO-NOccSO,NOccSO)
	if (CmplxConj == 0):
		z0.imag = 0
	sx = ao2mo(Sx,MO,2)
	s_x = MFEne(sx,Zeros,z0)
	sxg = MFGrad(sx,Zeros,z0,s_x)
	print("Sx and max |SG|", s_x.real, np.max(np.abs(sxg.real)), np.max(np.abs(sxg.imag)))
	return np.concatenate((sxg.real,sxg.imag))

def SzProj(Z):
	Znew = np.zeros([NVrtSO,NOccSO], dtype=complex)
	Znew[:NVrtA,:NOccA] = Z[:NVrtA,:NOccA]
	Znew[NVrtA:,NOccA:] = Z[NVrtA:,NOccA:]
	return Znew

def RHFProj(Z):
	Znew = np.zeros([NVrtSO,NOccSO], dtype=complex)
	Znew[:NVrtA,:NOccA] = Z[:NVrtA,:NOccA]
	Znew[NVrtA:,NOccA:] = Z[:NVrtA,:NOccA]
	return Znew

# The input MO should be RHF type
def FrozenProj(Z):
	Znew = np.zeros([NVrtSO,NOccSO], dtype=complex)
	Znew[:,:] = Z
	Znew[:,:NFO] = 0
	Znew[:,NOccAO:NOccAO+NFO] = 0
	return Znew

def PseudoDiag(HOne,HTwo,A2G):
	from makeH import MakeFock
	newMO = np.zeros_like(A2G)
	F = MakeFock(HOne,HTwo,A2G[:,:NOccSO])
	e, V = np.linalg.eigh(F)
	uocc = V[:,:NOccSO].T.conj() @ A2G[:,:NOccSO]
	R = subspace_eigh(np.diag(e[:NOccSO]), uocc)[1]
	newMO[:,:NOccSO] = V[:,:NOccSO] @ R
	uvir = V[:,NOccSO:].T.conj() @ A2G[:,NOccSO:]
	R = subspace_eigh(np.diag(e[NOccSO:]), uvir)[1]
	newMO[:,NOccSO:] = V[:,NOccSO:] @ R
	return newMO

def subspace_eigh(Fock, MO_Occ):
	f = MO_Occ.T.conj() @ Fock @ MO_Occ
	moe, u = np.linalg.eigh(f)
	return moe, MO_Occ @ u
