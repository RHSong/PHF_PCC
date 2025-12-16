import numpy as np
import scipy.sparse.linalg as lin
import time
from parameter_fc import *
from PHFTools import *
from PHFTools_F import phftools, buildsc_dr
dim = 1 + NOccSO * NVrtSO + NOccSO * (NOccSO-1) * NVrtSO * (NVrtSO-1) / 4
#dim = 1 + OV + OV**2
dim = int(dim)
def Lanczos(Rsp1, Rsp2, Rpg1, Rk1, fsp, fpg, fk):
	print("Total SD dimension is",dim)
	N = min(280,dim)
	nb = min(280,dim)
	N = N // nb
	print("Total block:", N, ", block dimension:", nb)
	tol = 1e-6
	Kspace = np.zeros([dim,N*nb], dtype=complex)
	qold = np.zeros([dim,nb], dtype=complex)
	qnew = np.eye(dim,nb, dtype=complex)
	A = np.zeros([nb,nb], dtype=complex)
	B = np.zeros([nb,nb], dtype=complex)
	T = np.zeros([(N+1)*nb,(N+1)*nb], dtype=complex)
	tmp = np.zeros([dim,nb], dtype=complex)
	for i in range(N):
		t1 = time.time()
		for j in range(nb):
			tmp[:,j] = BuildSC(Rsp1, Rsp2, Rpg1, Rk1, qnew[:,j], fsp, fpg, fk)
		tmp = tmp - qold @ B.T.conj()
		A = qnew.T.conj() @ tmp
		R = tmp - qnew @ A
		qold = qnew
		Kspace[:,i*nb:(i+1)*nb] = qold
		qnew, B = np.linalg.qr(R)
		T[i*nb:(i+1)*nb,i*nb:(i+1)*nb] = A
		T[(i+1)*nb:(i+2)*nb,i*nb:(i+1)*nb] = B
		T[i*nb:(i+1)*nb,(i+1)*nb:(i+2)*nb] = B.T.conj()
		t2 = time.time()
		print("Block time", i, t2-t1)
	evals, evecs = np.linalg.eigh(T[:N*nb,:N*nb])
	m = sum([1 for i in range(N*nb) if (abs(evals[i])<tol)])
	print("Linear Dependency number is",m)
	return evals[:m], Kspace @ evecs[:,:m]

def Davidson(Rsp1, Rsp2, Rpg1, Rk1, fsp, fpg, fk):
	print("Total SD dimension is",dim)
	N = min(200,dim)
	nb = min(50,dim)
	print("Krylov dimension:", N, ", block dimension:", nb)
	tol = 1e-6
	t = np.eye(dim,nb, dtype=complex)
	V = np.zeros([dim,N], dtype=complex)
	for m in range(nb,N+1,nb):
		t1 = time.time()
		V[:,m-nb:m] = t
		tmp = np.zeros([dim,m], dtype=complex)
		V[:,:m],R = np.linalg.qr(V[:,:m])
		for j in range(m):
			tmp[:,j] = BuildSC(Rsp1, Rsp2, Rpg1, Rk1, V[:,j], fsp, fpg, fk)
		T = V[:,:m].T.conj() @ tmp
		evals, evecs = np.linalg.eigh(T)
		for j in range(0,nb):
			w = tmp @ evecs[:,j]
			w -= evals[j] * V[:,:m] @ evecs[:,j]
			t[:,j] = w
		t2 = time.time()
		print("Block time", m//nb,t2-t1)
	m = sum([1 for i in evals if (abs(i)<tol)])
	print("Linear Dependency number is",m)
	return evals[:m], V @ evecs[:,:m]

def BuildSC(Rsp1, Rsp2, Rpg1, Rk1, Vecin, fsp, fpg, fk):
	v0, v1, v2 = Vec2Mat(Vecin)
	c0, c1, c2 = phftools.buildsc(v0,v1,v2,NOccSO,ngrid,J,nci,Rsp1,Rsp2,Rpg1,Rk1, \
	roota,rootb,rooty,weightsp,weightpg,fsp,fpg,fk,NSO,npg,ncisp,ncipg,ncik)
	CVec = Mat2Vec(c0,c1,c2)
	return CVec

def ArpackSolver(Rsp1, Rsp2, Rpg1, Rk1, fsp, fpg, fk):
	tol = 1e-4
	Av = lambda x: BuildSC(Rsp1, Rsp2, Rpg1, Rk1, x, fsp, fpg, fk)
	Ahv = lambda x: BuildSC(Rsp1, Rsp2, Rpg1, Rk1, x, fsp, fpg, fk)
	A = lin.LinearOperator((dim,dim), matvec=Av, rmatvec=Ahv, dtype=complex)
	t1 = time.time()
#	evals, evecs = lin.eigsh(A, k=100,sigma=0.0,which='SM',mode='normal',ncv=200)
#	evecs,evals,vh = lin.svds(A,k=100,which='SM',solver='lobpcg',ncv=250)
	t2 = time.time()
	print("LinDep time:",t2-t1)
	m = sum([1 for i in evals if (abs(i)<tol)])
	print("Linear Dependency number is",m)
	return evals[:m], evecs[:,:m]

def ArpackSolver2(Rsp1, Rsp2, Rpg1, Rk1, fsp, fpg, fk):
	tol = 1e-5
	Av = lambda x: BuildSC(Rsp1, Rsp2, Rpg1, Rk1, x, fsp, fpg, fk)
	Ahv = lambda x: BuildSC(Rsp1, Rsp2, Rpg1, Rk1, x, fsp, fpg, fk)
	A = lin.LinearOperator((dim,dim), matvec=Av, rmatvec=Ahv, dtype=complex)
	b = np.zeros(dim, dtype=complex)
	NGuess = 110
	V0 = np.eye(dim,NGuess, dtype=complex)
	evecs = np.zeros([dim,NGuess], dtype=complex)
	t1 = time.time()
	for i in range(NGuess):
		if (i % 10 == 0):
			print("i-th guess",i)
		evecs[:,i], exit_code = lin.cg(A, b, x0=V0[:,i])
		if (exit_code != 0):
			print("not converge")
	S = evecs.T.conj() @ evecs
	e,c = np.linalg.eigh(S)
	invS = [1/np.sqrt(abs(i)) for i in e if (np.abs(i) > tol)]
	Sinv = np.diag(invS)
	X = np.matmul(c[:,len(e)-len(invS):],Sinv)
	evecs = evecs @ X
	t2 = time.time()
	print("LinDep time:",t2-t1)
	print("Linear Dependency number is",len(evecs[0,:]))
	return len(evecs[0,:]), evecs

#def BuildSC(Rsp1, Rsp2, Rpg1, Rk1, Vecin, fsp, fpg, fk):
#	v0, v1, v2 = Vec2Mat(Vecin)
#	C0 = 0
#	C1 = 0
#	C2 = 0
#	Ovlp = 0
#	newMOs = np.eye(NSO,NOccSO)
#	for il in range(ngrid[0]):
#		for iy in range(ngrid[1]):
#			for ipg in range(npg):
#				for ik in range(nk):
#					R = Rsp1[il,:,:] @ Rsp2[iy,:,:]
#					R = R @ Rpg1[ipg,:,:]
#					R = Rk1[ik,:,:] @ R
#					a0 = v0
#					a1 = ao2mo(v1,R.T.conj(),2)
#					a2 = ao2mo(v2,R.T.conj(),4)
#					rho0, s0 = Rdm(R,newMOs[:,:NOccSO])
#					c0, c1, c2 = buildsc_dr.buildsc_vec(a0,a1,a2,rho0,NSO)
#					s0 = s0 * weightsp[il,iy]
#					c0 = c0 * s0
#					c1 = c1 * s0
#					c2 = c2 * s0
#					for ni in range(ncisp):
#						for nj in range(ncisp):
#							wig = WignerMat(J,ni-J,nj-J,roota[il],rootb[il],rooty[iy])
#							for p in range(ncipg):
#								for q in range(ncipg):
#									wpg = weightpg[ipg,p,q]
#									fac = fsp[ni].conj() * fsp[nj] * fpg[p].conj() * fpg[q] * fk[ik].conj()
#									C0 += c0 * wig * wpg * fac
#									C1 += c1 * wig * wpg * fac
#									C2 += c2 * wig * wpg * fac
#	CVec = Mat2Vec(C0,C1,C2) / Ovlp
#	return CVec


def Vec2Mat(Vec):
	v0 = Vec[0]
	v1 = np.zeros([NSO,NSO], dtype=complex)
	v2 = np.zeros([NSO,NSO,NSO,NSO], dtype=complex)
	tmp = 0
	for i in range(NOccSO):
		for a in range(NOccSO,NSO):
			tmp += 1
			v1[a,i] = Vec[tmp]
	for j in range(1,NOccSO):
		for i in range(j):
			for b in range(NOccSO+1,NSO):
				for a in range(NOccSO,b):
					tmp += 1
					v2[a,b,i,j] = Vec[tmp]
	return v0, v1, v2

def Mat2Vec(v0,v1,v2):
	Vec = np.zeros(dim, dtype=complex)
	Vec[0] = v0
	tmp = 0
	for i in range(NOccSO):
		for a in range(NOccSO,NSO):
			tmp += 1
			Vec[tmp] = v1[a,i]
	for j in range(1,NOccSO):
		for i in range(j):
			for b in range(NOccSO+1,NSO):
				for a in range(NOccSO,b):
					tmp += 1
					Vec[tmp] = v2[a,b,i,j]
	return Vec
