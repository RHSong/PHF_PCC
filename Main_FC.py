from parameter import *
import pickle
import time
import os
import pyscf
from scipy.linalg import block_diag
from FrozenCore import *
from MatchOrb import *
from PG_irrep import PG_Proj_SO
	
## RHF
RHF = scf.RHF(mol) if (mol.spin == 0) else scf.ROHF(mol)   # ROHF for open shells (name kept)
#RHF = scf.RHF(mol).sfx2c1e()
#RHF.get_hcore = lambda *args: h1
#RHF.get_ovlp = lambda *args: ovlp
#RHF._eri = pyscf.ao2mo.restore(8, eri, NAO)
#RHF.irrep_nelec = {'Ag': 14}
RHF.diis_space = 20
if (os.path.exists("RHF.p")):
	dm = pickle.load(open( "RHF.p", "rb" ))
	print("load RHF DM")
	RHF.kernel(dm)
else:
	RHF.kernel()
# check RHF stab
if (mol.spin == 0):
	RHFstab = RHF.stability(external=True)
else:   # PySCF has no external (ROHF -> UHF) stability analysis; internal only, UHF guess from the ROHF orbitals
	_mo = RHF.stability(internal=True, external=False)[0]
	RHFstab = (_mo, (_mo, _mo))
# save and print results
print("E(RHF)=",RHF.e_tot)
pickle.dump(RHF.make_rdm1(), open( "RHF.p", "wb" ))
if (os.path.exists("FCMO.p")):
	Ref = pickle.load(open( "FCMO.p", "rb" ))
	RHF.mo_coeff = MatchOrb(RHF.mo_coeff,Ref,NOccAO)
pickle.dump(RHF.mo_coeff, open( "FCMO.p", "wb" ))
# update h1 for x2c
h1 = RHF.get_hcore()

## FC
E0, H1, H2 =  FrozenCore(h1,eri,RHF.mo_coeff,NFC,NFV)
# redefine parameters
NAOnew = NAO - NFC - NFV
NOccSOnew = NOccSO - 2 * NFC
Enuc = Enuc + E0
print("FC ENuc=", Enuc)

## transfer PG symm operators
if (PG in ('C4v', 'D4h') and mol.groupname in ('Dooh', 'Coov')):
	# linear molecule (mol.symmetry = True): operators from the PySCF irrep labels of the (symmetry-adapted) MOs
	ncipg, weightpg, Rpg = PG_Proj_SO(mol, RHF, PG, Irrep)
else:
	ncipg, weightpg, Rpg = PG_Proj(mol, PG, Irrep, ovlp, RHF.mo_coeff, NAO, nx=nx, ny=ny, kx=kx, ky=ky)
Rpgnew = np.zeros([npg,2*NAOnew,2*NAOnew],dtype=complex)
for i in range(npg):
	tmp = Rpg[i,NFC:NAO-NFV,NFC:NAO-NFV]
	Rpgnew[i,:,:] = block_diag(tmp, tmp)
if (PG in ('C4v', 'D4h') and mol.groupname in ('Dooh', 'Coov')):
	# the point-group operators must commute with the Hamiltonian
	err = max(abs(Rpgnew[i,:NAOnew,:NAOnew].real.T @ H1 @ Rpgnew[i,:NAOnew,:NAOnew].real - H1).max() for i in range(npg))
	print("PG operators: max |R^T H1 R - H1| =", err)
	if (err > 1e-6):
		raise RuntimeError("point-group operators do not commute with H1 (%.2e)" % err)
pickle.dump([Enuc,H1,H2,ncipg, weightpg, Rpgnew,NAOnew,NOccSOnew,mol.spin],open("FC.p","wb"),protocol=4)   # 9th entry: 2S
