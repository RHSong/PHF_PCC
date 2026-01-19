import pickle
from pyscf import gto, scf, fci, cc, tdscf
from ao2mo import *
from makeH import *
from Spin_Proj import *
from PG_Proj import *
from Cmplx_Proj import *
mol = gto.Mole()
HamType = 'Hub'
if (HamType == 'Mol'):
	dista = 4.0 # unit Bohr or Ang
	a0 = 0.529177210903 # Bohr to Ang
	mol.atom = [['N',(0, 0, 0)], ['N',(0, 0, dista)]]
	mol.spin = 0
	mol.basis = 'cc-pvdz'
#	mol.symmetry = 'Dooh'
	mol.build()
	ovlp = mol.intor('int1e_ovlp')
	h1 = mol.intor('int1e_kin') + mol.intor('int1e_nuc')
	eri = mol.intor('int2e')
	NAO = mol.nao
	NOccSO = sum(mol.nelec)
	Enuc = mol.energy_nuc()
	nx = 0
	ny = 0
if (HamType == 'Hub'):
	nele = 6
	nx = 1
	ny = 6
	mol.nelectron = nele
	nsite = nx * ny
	U0 = 4
	h1, eri = BuildHamitonian(nx,ny,U0,pbcx=True,pbcy=True)
	ovlp = np.eye(nsite)
	Enuc = 0
	NAO = nsite
	NOccSO = nele
# redefine parameters
Enuc,h1,eri,ncipg, weightpg,Rpg, NAO,NOccSO = pickle.load(open( "FC.p", "rb" ))
mol.nelectron = NOccSO
ovlp = np.eye(NAO)
#ncipg = 1
#weightpg = np.ones([1,1,1])
#Rpg = np.zeros([1,2*NAO,2*NAO])
#Rpg[0,:,:] = np.eye(2*NAO)

# Do RHF
OrthAO = np.eye(NAO)
Xinv = np.linalg.inv(OrthAO)
Ovlp = ao2mo(ovlp,OrthAO,2)
h1 = ao2mo(h1,OrthAO,2)
eri = ao2mo(eri,OrthAO,4)
NSO = 2*NAO
NVrtSO = NSO - NOccSO
NOccAO = NOccSO // 2
OV = NOccSO * NVrtSO
NOccA = NOccAO
NVrtA = NAO-NOccA
''' 
Spin Projection, SP = 0, No projection
				 SP = 1, SGHF
				 SP = 2, SUHF
				 SP = 3, SzPHF
for SUHF & SzPHF, set grid pts as [1,n]
'''
SP = 2
ngrid = [1,8]
J, M = [0, 0]
ncisp, roota, rootb, rooty, weightsp, R1, R2 = Spin_Proj(SP, ngrid, NAO, NSO, J)
npg = len(weightpg)

'''
Complex Conj = 0, real
			 = 1, complex
			 = 2, complex projection
'''
is_RHF = False
CmplxConj = 0
if (CmplxConj == 2):
	ncik = 2
else:
	ncik = 1
nci = ncisp * ncipg * ncik
npoints = ngrid[0] * ngrid[1] * npg
# Broyden Space
nBroyVec = 30
