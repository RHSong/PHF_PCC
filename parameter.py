import pickle
from pyscf import gto, scf, fci, cc, tdscf
from ao2mo import *
from makeH import *
from Spin_Proj import *
from PG_Proj import *
from Cmplx_Proj import *
mol = gto.Mole()
HamType = 'Mol'
if (HamType == 'Mol'):
	dista = 1.3 # unit Bohr or Ang
	a0 = 0.529177210903 # Bohr to Ang
	mol.atom = [['N',(0, 0, 0)], ['N',(0, 0, dista)]]
	mol.spin = 0
	mol.basis = 'cc-pvdz'
	mol.symmetry = 'D2h'
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
# FC parameters
NFC = 0
NFV = 0
# Do RHF
OrthAO = getTrans(ovlp)
Xinv = np.linalg.inv(OrthAO)
NSO = 2*NAO
NVrtSO = NSO - NOccSO
NOccAO = NOccSO // 2
OV = NOccSO * NVrtSO
NOccA = NOccAO
NVrtA = NAO-NOccA
is_RHF = True
''' 
Spin Projection, SP = 0, No projection
				 SP = 1, SGHF
				 SP = 2, SUHF
				 SP = 3, SzPHF
for SUHF & SzPHF, set grid pts as [1,n]
'''
SP = 1
ngrid = [6,6]
J, M = [0, 0]
ncisp, roota, rootb, rooty, weightsp, R1, R2 = Spin_Proj(SP, ngrid, NAO, NSO, J)
'''
Point Group Projection

'''
PG = 'None'
Irrep = 'Ag'
kx = 1
ky = 1
ncipg, weightpg, Rpg = PG_Proj(mol, PG, Irrep, ovlp, OrthAO, NAO, nx=nx, ny=ny, kx=kx, ky=ky)
npg = len(weightpg)

'''
Complex Conj = 0, real
			 = 1, complex
			 = 2, complex projection
'''
CmplxConj = 0
if (CmplxConj == 2):
	ncik = 2
else:
	ncik = 1
nci = ncisp * ncipg * ncik
npoints = ngrid[0] * ngrid[1] * npg
# Broyden Space
nBroyVec = 30
