"""python semicanon_plusT_ufc.py FC_SO.p MO_used.p SP OUT.p   (PHF_UFC version of SUCCSD_new/semicanon_plusT.py)
Semicanonical orbitals of the PLAIN (bare) Fock inside the CC orbital basis (e.g. the PHF effective-Fock MOs), with the
Fock built exactly as run_pcc_t_ufc.py does (ufc_core.ghf_fock_so on load_fc(FC_SO.p or FC.p)).  oo and vv blocks diagonalized
(SP = 2: alpha and beta blocks separately, OA OB VA VB; otherwise one oo and one vv block); occupied/virtual spaces kept.
Writes OUT.p = [MO_sc, fsp, fpg, fk]; prints the largest off-diagonal oo/vv bare-Fock element before and after."""
import sys, os, pickle, numpy as np
here = os.path.dirname(os.path.abspath(__file__)); sys.path.insert(0, os.path.dirname(here)); sys.path.insert(0, here)   # repo root, then PHF_UFC
from ufc_core import load_fc, ghf_fock_so
from lean_ints import semicanon_from_fock
fcf, mof, SP, out = sys.argv[1], sys.argv[2], int(sys.argv[3]), sys.argv[4]
fc = load_fc(fcf); NAO, NOcc, Spin2 = int(fc['NAO']), int(fc['NOccSO']), int(fc['Spin2']); NOccA = (NOcc + Spin2) // 2
d = pickle.load(open(mof, 'rb')); MO = np.array(d[0], dtype=complex)
F = ghf_fock_so(fc['h_a'], fc['h_b'], fc['eri_aa'], fc['eri_ab'], fc['eri_bb'], MO[:, :NOcc], NAO)
def offd(C):
    Fm = C.conj().T @ F @ C; o = NOcc
    return max(np.abs(Fm[:o, :o] - np.diag(np.diag(Fm[:o, :o]))).max(), np.abs(Fm[o:, o:] - np.diag(np.diag(Fm[o:, o:]))).max())
MOsc = semicanon_from_fock(F, MO, NOcc, SP, NOccA)
P0 = MO[:, :NOcc] @ MO[:, :NOcc].conj().T; P1 = MOsc[:, :NOcc] @ MOsc[:, :NOcc].conj().T
print("%s: max|offdiag oo/vv bare Fock| CC basis %.2e -> semicanonical %.2e   occupied-space change %.1e"
      % (mof, offd(MO), offd(MOsc), np.abs(P1 - P0).max()))
pickle.dump([MOsc] + list(d[1:4]), open(out, 'wb'))
