"""PHF_UFC step 4: parameters of the PHF / PCC step in the ACTIVE spin-orbital basis of FC_SO.p (ufc_core.py).
Copy this file (with parameter.py, the drivers of PHF_UFC and FC_SO.p / GHFMO.p) into the case folder and edit the
projection settings below (SP, ngrid, J, CmplxConj) as in the usual parameter_fc.py.
Differences from the standard parameter_fc.py: h1 / eri are NOT spatial here -- the drivers build HOne / HTwo from
FC_SO (ufc_core.hone_so / htwo_mulliken_so); the spin rotations R1, R2 are transformed with ufc_core.spin_rot_so."""
import pickle
import numpy as np
from pyscf import gto, scf
from ao2mo import *
from makeH import *
from Spin_Proj import *
from PG_Proj import *
from Cmplx_Proj import *
from ufc_core import is_fc_so, spin_rot_list_so, hone_so, htwo_mulliken_so, ghf_fock_so, dirac_mo_eri_so
FC_SO = pickle.load(open("FC_SO.p", "rb"))
assert is_fc_so(FC_SO), "FC_SO.p is not a spin-orbital frozen-core file (run ufc_core.py)"
Enuc = FC_SO['Enuc']; NAO = int(FC_SO['NAO']); NOccSO = int(FC_SO['NOccSO']); Spin2 = int(FC_SO.get('Spin2', 0))
ncipg = int(FC_SO['ncipg']); weightpg = np.asarray(FC_SO['weightpg']); Rpg = np.array(FC_SO['Rpg'], dtype=complex)
mol = gto.Mole(); mol.nelectron = NOccSO
h1 = None; eri = None            # no spatial integrals in this basis: use FC_SO (see ufc_core.py)
ovlp = np.eye(NAO); OrthAO = np.eye(NAO); Xinv = np.eye(NAO); Ovlp = np.eye(NAO)
NSO = 2 * NAO
NVrtSO = NSO - NOccSO
NOccA = (NOccSO + Spin2) // 2
NOccB = (NOccSO - Spin2) // 2
NOccAO = NOccSO // 2
OV = NOccSO * NVrtSO
NVrtA = NAO - NOccA
NVrtB = NAO - NOccB
'''
Spin Projection, SP = 0, No projection
				 SP = 1, SGHF
				 SP = 2, SUHF
				 SP = 3, SzPHF
for SUHF & SzPHF, set grid pts as [1,n]
'''
SP = 2
ngrid = [1, 12]
J, M = [Spin2 / 2, 0]
ncisp, roota, rootb, rooty, weightsp, R1, R2 = Spin_Proj(SP, ngrid, NAO, NSO, J)
R1 = spin_rot_list_so(R1, FC_SO['O_ab'], NAO)      # product-form rotations -> UHF spin-orbital basis (truncated)
R2 = spin_rot_list_so(R2, FC_SO['O_ab'], NAO)
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
