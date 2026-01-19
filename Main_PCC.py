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
fcomm = MPI.COMM_WORLD.py2f()

SGHFMO, fsp, fpg, fk = pickle.load(open( "SGHF.p", "rb" ))
#SGHFMO = pickle.load(open( "GHFMO.p", "rb" ))
HOne = np.zeros([NSO,NSO])
HTwo = np.zeros([NSO,NSO,NSO,NSO])
HOne[:NAO,:NAO] = h1[:,:]
HOne[NAO:,NAO:] = h1[:,:]
HTwo[:NAO,:NAO,:NAO,:NAO] = eri[:,:,:,:]
HTwo[NAO:,NAO:,NAO:,NAO:] = eri[:,:,:,:]
HTwo[NAO:,NAO:,:NAO,:NAO] = eri[:,:,:,:]
HTwo[:NAO,:NAO,NAO:,NAO:] = eri[:,:,:,:]
HTwo = Mulliken2Dirac(HTwo)
X = block_diag(OrthAO, OrthAO)
Xinv = block_diag(Xinv, Xinv)
SGHFMO = SemiCanon(HOne, HTwo, SGHFMO, NOccSO, SP)

# PCC
H1 = ao2mo(HOne,SGHFMO,2)
H2 = ao2mo(HTwo,SGHFMO,4)
HOne = None
HTwo = None
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
PCC = pgcc(SGHFMO,H1,H2,NAO,NOccSO,J,CmplxConj,SP,ngrid,R1,R2,Rpg,Rk,roota,rootb,rooty,weightsp, \
			weightpg,fsp,fpg,fk,nBroyVec,X,Xinv,Enuc,fcomm,NSO,npg,ncisp,ncipg,ncik)
print("E(PCC)=",PCC[0].real)

