from parameter_fc import *
from Spin import *
#from LinDep import *
#from PUCC_Exact import gcc
#from PUCC_X3 import pgcc
from PUCC_X2 import pgcc
import pickle
import time
import os
import pyscf
from mpi4py import MPI
from scipy.linalg import block_diag
from lean_ints import unpack_eri, ghf_fock_ao, semicanon_from_fock, dirac_mo_eri_F
fcomm = MPI.COMM_WORLD.py2f()

SGHFMO, fsp, fpg, fk = pickle.load(open( "SGHF.p", "rb" ))[:4]   # 5th entry (PHFFock moe) ignored here
#SGHFMO = pickle.load(open( "GHFMO.p", "rb" ))
HOne = np.zeros([NSO,NSO])
HOne[:NAO,:NAO] = h1[:,:]
HOne[NAO:,NAO:] = h1[:,:]
X = block_diag(OrthAO, OrthAO)
Xinv = block_diag(Xinv, Xinv)
# Lean integrals (lean_ints.py): no spin-orbital AO tensor, one complex NSO^4 array
# (the final <pq||rs>, Fortran order -> passed to PGCC without an f2py copy).
eri4 = unpack_eri(eri, NAO)
SGHFMO = np.array(SGHFMO, dtype=complex)
SGHFMO = semicanon_from_fock(ghf_fock_ao(h1, eri4, SGHFMO[:,:NOccSO], NAO), SGHFMO, NOccSO, SP)

# PCC
H1 = ao2mo(HOne,SGHFMO,2)
H2 = dirac_mo_eri_F(eri4, SGHFMO, NAO)
HOne = None
eri4 = None
eri = None
for i in range(ngrid[0]):
	R1[i,:,:] = ao2mo(R1[i,:,:],SGHFMO,2)
for i in range(ngrid[1]):
	R2[i,:,:] = ao2mo(R2[i,:,:],SGHFMO,2)
for i in range(npg):
	Rpg[i,:,:] = ao2mo(Rpg[i,:,:],SGHFMO,2)
Rk = CmplxProj(SGHFMO,NSO,ncik)
#print("Ovlp=", EvalOvlp(R1,R2,Rpg,Rk,fsp,fpg,fk))
if (ncik == 1):
	fk = np.ones(1)
nBroyVec = 20      # Broyden history (was 30): saves ~1.4 GB per 10 vectors for Cr2 NFC=10
PCC = pgcc(SGHFMO,H1,H2,NAO,NOccSO,J,CmplxConj,SP,ngrid,R1,R2,Rpg,Rk,roota,rootb,rooty,weightsp, \
			weightpg,fsp,fpg,fk,nBroyVec,X,Xinv,Enuc,fcomm,NSO,npg,ncisp,ncipg,ncik)
print("E(PCC)=",PCC[0].real)

