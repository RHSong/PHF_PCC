import numpy as np
import sys
from scipy.optimize import minimize, basinhopping
from parameter import *
from SUHF import Thouless2MOs
from SUHFTools import MFEne, MFGrad
def optGHF(HOne, HTwo, MOs):
	z0 = np.zeros((NSO-NOccSO)*NOccSO, dtype=complex)
	guess = np.concatenate((z0.real,z0.imag))
	minimizer_kwargs = {"method":"BFGS", "jac":True, "args":(HOne,HTwo)}
	res = basinhopping(EandG, guess, minimizer_kwargs=minimizer_kwargs, niter=2, T=0.5)
	sol = res.x[:NVrtSO*NOccSO] + res.x[NVrtSO*NOccSO:] * 1j
	z = sol.reshape([NSO-NOccSO,NOccSO])
	newMOs = Thouless2MOs(z)
	newMOs = MOs @ newMOs
	return res.fun, newMOs

def EandG(Z, *args):
	z0 = Z[:NVrtSO*NOccSO] + Z[NVrtSO*NOccSO:] * 1j
	sol = z0.reshape(NSO-NOccSO,NOccSO)
	HOne,HTwo = args
	E0 = MFEne(HOne,HTwo,sol)
	if (abs(E0.imag) > 1e-8):
		sys.exit("Ene is imaginary")
	G = MFGrad(HOne,HTwo,sol,E0)
	G = G.reshape((NSO-NOccSO)*NOccSO)
	return E0.real, np.concatenate((G.real,G.imag))
