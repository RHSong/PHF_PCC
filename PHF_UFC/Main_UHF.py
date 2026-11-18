"""PHF_UFC step 2: RHF and (AFM) UHF of the molecule with pyscf directly (not from FC.p).

Reads parameter.py (mol, NFC) and FCMO.p (the RHF orbitals of Main_FC_full.py, so that the RHF basis of FC.p and the
RHF used here are the same object).  Runs the all-electron UHF exactly as Main_HF.py does (RHF -> RHF.stability(external=True) -> UHF with diis_space 10 and
level_shift 0.1 from the instability orbitals; UHF stability only reported), with the guess from
    UHFdm_ae.p   (AO density, if present)                      else
    ATOM_GUESS   (parameter.py: dict {'spin': 2S_atom, 'irrep_nelec': {...}} -> two high-spin atoms, spins opposed)  else
    the RHF -> UHF instability of the RHF.
Saves MOL_MO.p = {'C_rhf', 'C_a', 'C_b', 'e_a', 'e_b', 'S', 'E_RHF', 'E_UHF', 'S2'} (AO coefficients, columns by energy)."""
from parameter import *
import pickle, os
from scipy.linalg import block_diag
X2C = globals().get('X2C', True)          # Cr2 uses sfX2C-1e; set X2C = False in parameter.py otherwise
def _x2c(mf): return mf.sfx2c1e() if X2C else mf
C_rhf = pickle.load(open('FCMO.p', 'rb'))
S = mol.intor('int1e_ovlp')
# point-group symmetry is ignored from the RHF on (a symmetry-adapted RHF only probes symmetry-preserving rotations in
# the stability analysis and misses the AFM instability); the RHF is re-solved from the FCMO.p density
molu = mol.copy(); molu.symmetry = False; molu.build()
RHF = _x2c(scf.RHF(molu))
dm_r = C_rhf[:, :NOccSO // 2] @ C_rhf[:, :NOccSO // 2].T * 2
E_RHF = RHF.kernel(dm0=dm_r); print('E(RHF, from the FCMO.p density, no symmetry)= %.10f  converged %s' % (E_RHF, RHF.converged))
RHFstab = RHF.stability(external=True)                        # as in Main_HF.py
UHF = _x2c(scf.UHF(molu)); UHF.diis_space = 10; UHF.level_shift = 0.1   # Main_HF.py settings
if os.path.exists('UHFdm_ae.p'):
    dm0 = pickle.load(open('UHFdm_ae.p', 'rb')); print('UHF guess: UHFdm_ae.p')
elif 'ATOM_GUESS' in globals():
    ag = ATOM_GUESS; sym = mol.atom_symbol(0)
    atom = gto.M(atom='%s 0 0 0' % sym, basis=mol.basis, spin=ag['spin'], symmetry='D2h', verbose=0)
    mfa = _x2c(scf.UHF(atom)); mfa.irrep_nelec = ag['irrep_nelec']
    ea = mfa.kernel(); da, db = mfa.make_rdm1(); print('UHF guess: two %s atoms (E_atom = %.10f), spins opposed' % (sym, ea))
    dm0 = np.array([block_diag(da, db), block_diag(db, da)])
else:
    mo_u = RHFstab[1]                                            # orbitals rotated along the RHF -> UHF instability (Main_HF.py)
    dm0 = np.array([mo_u[0][:, :mol.nelec[0]] @ mo_u[0][:, :mol.nelec[0]].T, mo_u[1][:, :mol.nelec[1]] @ mo_u[1][:, :mol.nelec[1]].T])
    print('UHF guess: RHF -> UHF external instability')
E_UHF = UHF.kernel(dm0=dm0)
UHFstab = UHF.stability(external=True)   # as in Main_HF.py: the stability is only reported
S2 = UHF.spin_square()[0]
print('E(UHF all-electron)= %.10f  <S^2> %.4f  converged %s' % (E_UHF, S2, UHF.converged))
pickle.dump(UHF.make_rdm1(), open('UHFdm_ae.p', 'wb'))
Ca, Cb = UHF.mo_coeff; ea, eb = UHF.mo_energy
ia, ib = np.argsort(ea), np.argsort(eb)
pickle.dump({'C_rhf': C_rhf, 'C_a': Ca[:, ia], 'C_b': Cb[:, ib], 'e_a': ea[ia], 'e_b': eb[ib], 'S': S, 'E_RHF': E_RHF, 'E_UHF': E_UHF, 'S2': S2,
             'nelec': mol.nelec}, open('MOL_MO.p', 'wb'))
print('MOL_MO.p written')
