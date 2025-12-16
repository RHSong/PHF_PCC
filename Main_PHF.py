from parameter_fc import *
from Spin import *
from PHF import *
from LinDep import *
from makeH import SemiCanon
import pickle
import time
import os
import pyscf
from mpi4py import MPI
from scipy.linalg import block_diag
fcomm = MPI.COMM_WORLD.py2f()

# PHF
t1 = time.time()
if (os.path.exists("SGHF.p")):
	A2G, fsp, fpg, fk = pickle.load(open( "SGHF.p", "rb" ))
else:
	A2G = pickle.load(open( "GHFMO.p", "rb" ))
#A2G = np.eye(NSO)
#A2G = SortOrb(A2G,NAO,NOccAO,1)
HOne = np.zeros([NSO,NSO])
HTwo = np.zeros([NSO,NSO,NSO,NSO])
HOne[:NAO,:NAO] = h1[:,:]
HOne[NAO:,NAO:] = h1[:,:]
HTwo[:NAO,:NAO,:NAO,:NAO] = eri[:,:,:,:]
HTwo[NAO:,NAO:,NAO:,NAO:] = eri[:,:,:,:]
HTwo[NAO:,NAO:,:NAO,:NAO] = eri[:,:,:,:]
HTwo[:NAO,:NAO,NAO:,NAO:] = eri[:,:,:,:]
HTwo = Mulliken2Dirac(HTwo)

H1 = ao2mo(HOne,A2G,2)
H2 = ao2mo(HTwo,A2G,4)
EPHF, SGHFMO, fsp, fpg, fk = optPHF(H1,H2,A2G,fcomm)
#SGHFSz,SGHFSS = SGHFS2(SGHFMO,f0,NOccSO,NAO)
GHFS = calcS(SGHFMO,NOccSO,NAO)
GHFSS = calcS2(SGHFMO,NOccSO,NAO)
#####S2ghf = SGHFS2(SGHFMO,NOccSO,NAO)
#####print("S^2(GHF)=",S2ghf)
t2 = time.time()
print("E(PHF)=",EPHF + Enuc)
print("GHF S^2, Sxyz=", GHFSS,GHFS)
#print("SGHF S^2, Sz=", SGHFSS.real,SGHFSz.real)
print("PHF time=",t2-t1)
pickle.dump([SGHFMO, fsp, fpg, fk], open( "SGHF.p", "wb" ))
G2S = A2G.T @ SGHFMO
ovlp, Z = GetThouless(G2S)
print("ovlp=", ovlp)
for i in range(ngrid[0]):
	R1[i,:,:] = ao2mo(R1[i,:,:],SGHFMO,2)
for i in range(ngrid[1]):
	R2[i,:,:] = ao2mo(R2[i,:,:],SGHFMO,2)
for i in range(npg):
	Rpg[i,:,:] = ao2mo(Rpg[i,:,:],SGHFMO,2)
Rk = CmplxProj(SGHFMO,NSO,ncik)
print("Ovlp=", EvalOvlp(R1,R2,Rpg,Rk,fsp,fpg,fk))

