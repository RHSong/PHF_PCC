from parameter_fc import *
from Spin import *
from PHFTools import *
import pickle
import time
import os
import pyscf
from scipy.linalg import block_diag
	
## RHF
RHF = scf.RHF(mol)
RHF.get_hcore = lambda *args: h1
RHF.get_ovlp = lambda *args: Ovlp
RHF._eri = pyscf.ao2mo.restore(8, eri, NAO)
RHF.diis_space = 10
if (os.path.exists("RHFdm.p")):
	dm = pickle.load(open( "RHFdm.p", "rb" ))
	print("load RHF DM")
	RHF.kernel(dm)
else:
	RHF.kernel()
# check RHF stab
RHFstab = RHF.stability(external=True)
# save and print results
print("E(RHF)=",RHF.e_tot + Enuc)
pickle.dump(RHF.make_rdm1(), open( "RHFdm.p", "wb" ))
A2R = np.zeros([NSO,NSO])
A2R[:NAO,:NAO] = A2R[NAO:NSO,NAO:NSO] = RHF.mo_coeff

# FCI
#mci = fci.direct_spin1.FCI()
#mci = fci.addons.fix_spin_(mci, ss=J*(J+1))
#fcie, fcivec = mci.kernel(h1, eri, NAO, NOccSO, nroots=4)
#print("E(FCI1)= %8.8f" %(fcie[0]+Enuc))
#print("E(FCI2)= %8.8f" %(fcie[1]+Enuc))

## UHF
UHF = scf.UHF(mol)
UHF.get_hcore = lambda *args: h1
UHF.get_ovlp = lambda *args: Ovlp
UHF._eri = pyscf.ao2mo.restore(8, eri, NAO)
UHF.diis_space = 10
UHF.level_shift = 0.1
if (os.path.exists("UHFdm.p")):
	dm = pickle.load(open( "UHFdm.p", "rb" ))
	print("load UHF DM")
	UHF.kernel(dm)
else:	
	guessUHF = RHFstab[1]
	dm1 = np.dot(guessUHF[0][:,:NOccAO],guessUHF[0][:,:NOccAO].T)
	dm2 = np.dot(guessUHF[1][:,:NOccAO],guessUHF[1][:,:NOccAO].T)
	dm = np.array([dm1,dm2])
	UHF.kernel(dm)
# check UHF stab
UHFstab = UHF.stability(external=True)
#guessUHF = UHFstab[0]
#dm1 = np.dot(guessUHF[0][:,:NOccAO],guessUHF[0][:,:NOccAO].T)
#dm2 = np.dot(guessUHF[1][:,:NOccAO],guessUHF[1][:,:NOccAO].T)
#dm = np.array([dm1,dm2])
#UHF.kernel(dm)
# Print results
print("E(UHF)=",UHF.e_tot + Enuc)
pickle.dump(UHF.make_rdm1(), open( "UHFdm.p", "wb" ))
A2U = np.zeros([NSO,NSO])
A2U[:NAO,:NAO] = UHF.mo_coeff[0]
A2U[NAO:NSO,NAO:NSO] = UHF.mo_coeff[1]
R2U = A2R.T @ A2U
newA2R = SortOrb(A2R,NAO,NOccAO,1)
newR2U = SortOrb(R2U,NAO,NOccAO,1)
newA2U = SortOrb(A2U,NAO,NOccAO,1)

## GHF
GHF = scf.GHF(mol)
GHF.get_hcore = lambda *args: block_diag(h1,h1)
GHF.get_ovlp = lambda *args: block_diag(Ovlp,Ovlp)
GHF._eri = pyscf.ao2mo.restore(8, eri, NAO)
GHF.diis_space = 10
GHF.level_shift = 0.1
if (os.path.exists("GHFdm.p")):
	dm = pickle.load(open( "GHFdm.p", "rb" ))
	print("load GHF DM")
	GHF.kernel(dm)
else:
	guessGHF = UHFstab[1]
	guessGHF = SortOrb(guessGHF,NAO,NOccAO,1)
	dm = np.dot(guessGHF[:,:NOccSO],guessGHF[:,:NOccSO].T)
	GHF.kernel(dm)
GHFstab = GHF.stability(external=True)
print("E(GHF)=",GHF.e_tot + Enuc)
pickle.dump(GHF.make_rdm1(), open( "GHFdm.p", "wb" ))
#A2G = GHF.mo_coeff
A2G = newA2U
#A2G = MatchOrb(newA2U)
if (is_RHF):
	pickle.dump(newA2R, open( "GHFMO.p", "wb" ))
if (SP == 2):
	pickle.dump(newA2U, open( "GHFMO.p", "wb" ))
else:
	pickle.dump(GHF.mo_coeff, open( "GHFMO.p", "wb" ))

## GCCSD
#gcc = cc.UCCSD(UHF)
#if (os.path.exists("t1t2.p")):
#	t1, t2 = pickle.load(open( "t1t2.p", "rb" ))
#	print("load GCCSD T")
#	gcc.kernel(t1,t2)
#else:
#	gcc.kernel()
#pickle.dump([gcc.t1,gcc.t2],open( "t1t2.p", "wb" ))
#print('E(GCCSD)= %8.8f' %(gcc.e_tot + Enuc))
#et = gcc.ccsd_t()
#print('E(GCCSD-T)= %8.8f' %(et + gcc.e_tot + Enuc))

