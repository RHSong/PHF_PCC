import numpy as np
from copy import deepcopy
from scipy.optimize import minimize
from parameter import *
from ao2mo import *
from PHFTools import *
from Spin import *
def FixGauge(MOs):
	S2One, S2Two = BuildS2(NAO)
	S1 = ao2mo(S2One,MOs,2)
	S2 = ao2mo(S2Two,MOs,4)
	guess = np.zeros(3)
	res = minimize(SandG, guess, args=(S1,S2,MOs), method='BFGS', jac=True)
#	a,b,ib,y,iy = res.x
	a,ib,y = res.x
	b = 0
	iy = 0
	A = complex(a,0)
	B = complex(b,ib)
	Y = complex(y,iy)
	R = BuildRealRotMat(A,B,Y)
	R = ao2mo(R,MOs,2)
	ovlp, Z = GetThouless(R)
	newMOs = Thouless2MOs(Z)
	return MOs @ newMOs

def SandG(x0, *args):
	SOne, STwo, MOs = args
	SS = SVal(x0,SOne,STwo,MOs)
	glist = np.zeros([3,2])
	dx = 1e-3
	s2 = J*(J+1)
	for n in range(len(glist[0,:])):
		for i in range(len(glist[:,0])):
			xp = deepcopy(x0)
			xn = deepcopy(x0)
			xp[i] += (n+1) * dx
			xn[i] -= (n+1) * dx
			fp = SVal(xp,SOne,STwo,MOs)
			fn = SVal(xn,SOne,STwo,MOs)
			glist[i,n] = (fp - fn) / (2 * (n+1) * dx)
	Sg = (4 *  glist[:,0] - glist[:,1]) / 3
	print("S and max |G|", SS-s2, np.max(np.abs(Sg)))
	return (SS-s2)**2, 2*(SS-s2)*Sg


def SVal(x0,SOne,STwo,MOs):
#	a,b,ib,y,iy = x0
	a,ib,y = x0
	b = 0
	iy = 0
	A = complex(a,0)
	B = complex(b,ib)
	Y = complex(y,iy)
	R = BuildRealRotMat(A,B,Y)
	R = ao2mo(R,MOs,2)
	ovlp, Z = GetThouless(R)
	SS = MFEne(SOne,4*STwo,Z)
	return SS.real

def BuildRealRotMat(a,b,y):
	Ap = np.exp(a/2)
	An = np.exp(-a/2)
	Bp = (np.exp(b/2) + np.exp(-b/2)) / 2
	Bn = (np.exp(b/2) - np.exp(-b/2)) / 2j
	Yp = np.exp(y/2)
	Yn = np.exp(-y/2)
	SaRotMat = np.zeros([NSO,NSO],dtype=complex)
	SbRotMat = np.zeros([NSO,NSO],dtype=complex)
	SyRotMat = np.zeros([NSO,NSO],dtype=complex)
	for i in range(NAO):
		SaRotMat[i][i] = Ap
		SaRotMat[i+NAO][i+NAO] = An
		SbRotMat[i][i] = SbRotMat[i+NAO][i+NAO] = Bp
		SbRotMat[i][i+NAO] = Bn
		SbRotMat[i+NAO][i] = -Bn
		SyRotMat[i][i] = Yp
		SyRotMat[i+NAO][i+NAO] = Yn
	RotMat = SaRotMat @ SbRotMat @ SyRotMat
	return RotMat

def FixCmplxGauge(fk,MO):
	phi = -np.angle(fk) / (2 * NOccSO)
	newMO = np.exp(1j * phi) * np.eye(NSO)
	return newMO @ MO


def FixSz(MOs, tol=1e-10):
	"""Gauge fix using the Sz generator only: apply the non-unitary boost exp(a*Sz)
	(real a; BuildRealRotMat(a,0,0)), which leaves the J=0 projected state and the PHF
	energy invariant, and choose a such that <Sz> of the re-orthonormalised determinant
	vanishes. <Sz>(a) is monotonic in a (its derivative is the Sz variance), so a is
	found by bracketing + brentq. Only one-body quantities are needed."""
	from scipy.optimize import brentq
	sz = ao2mo(BuildSz(NAO), MOs, 2)
	def SzOf(a):
		R = ao2mo(BuildRealRotMat(complex(a, 0), 0, 0), MOs, 2)
		ovlp, Z = GetThouless(R)
		newMOs = np.eye(NSO, NOccSO, dtype=complex)
		newMOs[NOccSO:, :NOccSO] = Z
		rdm, o = Rdm(np.eye(NSO), newMOs)
		return np.einsum('ij,ji', sz, rdm).real
	s0 = SzOf(0.0)
	if (abs(s0) < tol):
		return MOs
	lo, hi = -0.5, 0.5
	while SzOf(lo) > 0: lo *= 2
	while SzOf(hi) < 0: hi *= 2
	a = brentq(SzOf, lo, hi, xtol=1e-12)
	print("FixSz: <Sz> %.6f -> %.2e  (a = %.6f)" % (s0, SzOf(a), a))
	R = ao2mo(BuildRealRotMat(complex(a, 0), 0, 0), MOs, 2)
	ovlp, Z = GetThouless(R)
	return MOs @ Thouless2MOs(Z)
