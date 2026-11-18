"""Molecule, frozen core and point group for Main_FC.py (step 1).  N2 / 6-31G at R = 3.0 bohr, N 1s frozen, C4 projection."""
import pickle
from pyscf import gto, scf
from ao2mo import *
from makeH import *
from Spin_Proj import *
from PG_Proj import *
from Cmplx_Proj import *

# molecule
mol = gto.Mole()
dista = 3.0                      # bond length (bohr)
mol.atom = [['N', (0, 0, 0)], ['N', (0, 0, dista)]]
mol.unit = 'Bohr'
mol.basis = '6-31g'
mol.spin = 0                     # 2S = n_alpha - n_beta (open shells: ROHF reference)
mol.symmetry = True              # pyscf Dooh: needed for the C4v / C4 point-group projection
mol.build()
ovlp = mol.intor('int1e_ovlp')
h1 = mol.intor('int1e_kin') + mol.intor('int1e_nuc')
eri = mol.intor('int2e')
NAO = mol.nao
NOccSO = sum(mol.nelec)
Enuc = mol.energy_nuc()
nx = 0
ny = 0

# frozen core: NFC lowest RHF orbitals (here the two 1s combinations), NFV highest virtuals
NFC = 2
NFV = 0

OrthAO = getTrans(ovlp)
Xinv = np.linalg.inv(OrthAO)
NSO = 2 * NAO
NVrtSO = NSO - NOccSO
Spin2 = mol.spin
OV = NOccSO * NVrtSO
NOccA = (NOccSO + Spin2) // 2
NOccB = (NOccSO - Spin2) // 2
NVrtA = NAO - NOccA
NVrtB = NAO - NOccB
NOccAO = NOccA
is_RHF = (Spin2 == 0)

# spin projection is set in parameter_fc.py; here only a placeholder for the imports of Main_FC.py
SP = 0
ngrid = [1, 1]
J, M = [Spin2 / 2, Spin2 / 2]
ncisp, roota, rootb, rooty, weightsp, R1, R2 = Spin_Proj(SP, ngrid, NAO, NSO, J)

# point-group projection: 'None', or for linear molecules (mol.symmetry = True) 'C4v' / 'D4h' with a C4v / D4h irrep
PG = 'C4v'
Irrep = 'A1'                      # Sigma_g+ ground state -> A1 of C4v -> A of C4 (cut in step 1 of run.sh)
kx = 1
ky = 1
ncipg, weightpg, Rpg = PG_Proj(mol, PG, Irrep, ovlp, OrthAO, NAO, nx=nx, ny=ny, kx=kx, ky=ky)
npg = len(weightpg)

CmplxConj = 0
ncik = 1
nci = ncisp * ncipg * ncik
npoints = ngrid[0] * ngrid[1] * npg
nBroyVec = 30
