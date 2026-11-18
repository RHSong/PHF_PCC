"""Projection settings for the HF / PHF steps (Main_HF.py, Main_PHF.py, PHFFock.py); integrals from FC.p."""
import pickle
from pyscf import gto, scf
from ao2mo import *
from makeH import *
from Spin_Proj import *
from PG_Proj import *
from Cmplx_Proj import *

# frozen-core Hamiltonian of Main_FC.py, in the orthonormal RHF orbital basis
Enuc, h1, eri, ncipg, weightpg, Rpg, NAO, NOccSO, Spin2 = pickle.load(open("FC.p", "rb"))[:9]
mol = gto.Mole()
mol.nelectron = NOccSO
mol.spin = Spin2
ovlp = np.eye(NAO)
OrthAO = np.eye(NAO)
Xinv = np.linalg.inv(OrthAO)
Ovlp = ao2mo(ovlp, OrthAO, 2)
h1 = ao2mo(h1, OrthAO, 2)
eri = ao2mo(eri, OrthAO, 4)
NSO = 2 * NAO
NVrtSO = NSO - NOccSO
OV = NOccSO * NVrtSO
NOccA = (NOccSO + Spin2) // 2
NOccB = (NOccSO - Spin2) // 2
NVrtA = NAO - NOccA
NVrtB = NAO - NOccB
NOccAO = NOccA
is_RHF = False

# spin projection: SP = 0 none, 1 SGHF (grid [Lebedev, gamma]), 2 SUHF (grid [1, n]), 3 Sz
SP = 2
ngrid = [1, 8]
J, M = [Spin2 / 2, Spin2 / 2]
ncisp, roota, rootb, rooty, weightsp, R1, R2 = Spin_Proj(SP, ngrid, NAO, NSO, J)
npg = len(weightpg)          # point-group projector from FC.p (identity unless PG was set in parameter.py)

# complex conjugation: 0 real, 1 complex, 2 complex projection
CmplxConj = 0
ncik = 2 if CmplxConj == 2 else 1
nci = ncisp * ncipg * ncik
npoints = ngrid[0] * ngrid[1] * npg
nBroyVec = 30
