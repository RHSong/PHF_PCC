"""python semicanon_plusT.py FC.p MO_used.p SP OUT.p   (orbitals for +(T))
Semicanonical orbitals of the PLAIN (bare) Fock inside the orbital basis of MO_used.p (the CC basis, e.g. the PHF
effective-Fock MOs): F = Fock of the reference determinant (first NOcc columns), oo and vv blocks diagonalized
(SP = 2: alpha and beta blocks separately, OA OB VA VB layout; otherwise one oo and one vv block).  The occupied and
virtual spaces are unchanged, so run_pcc_t.py --mo MO_used.p --newmo OUT.p rotates T1/T2 into them and builds the ERIs there.
Writes OUT.p = [MO_sc, fsp, fpg, fk] and prints the largest off-diagonal oo/vv bare-Fock element before and after."""
import sys, os, pickle, numpy as np
root = os.path.dirname(os.path.dirname(os.path.abspath(__file__))); sys.path.insert(0, root)   # repo root (this file: SUCCSD_new/)
from lean_ints import unpack_eri, ghf_fock_ao, semicanon_from_fock
fcf, mof, SP, out = sys.argv[1], sys.argv[2], int(sys.argv[3]), sys.argv[4]
fc = pickle.load(open(fcf, 'rb')); h1, eri, NAO, NOcc = fc[1], fc[2], int(fc[6]), int(fc[7])
Spin2 = int(fc[8]) if len(fc) > 8 else 0; NOccA = (NOcc + Spin2) // 2
d = pickle.load(open(mof, 'rb')); MO = np.array(d[0], dtype=complex)
F = ghf_fock_ao(h1, unpack_eri(eri, NAO), MO[:, :NOcc], NAO)
def offd(C):
    Fm = C.conj().T @ F @ C; o = NOcc
    return max(np.abs(Fm[:o, :o] - np.diag(np.diag(Fm[:o, :o]))).max(), np.abs(Fm[o:, o:] - np.diag(np.diag(Fm[o:, o:]))).max())
MOsc = semicanon_from_fock(F, MO, NOcc, SP, NOccA)
P0 = MO[:, :NOcc] @ MO[:, :NOcc].conj().T; P1 = MOsc[:, :NOcc] @ MOsc[:, :NOcc].conj().T
print("%s: max|offdiag oo/vv bare Fock| CC basis %.2e -> semicanonical %.2e   occupied-space change %.1e"
      % (mof, offd(MO), offd(MOsc), np.abs(P1 - P0).max()))
pickle.dump([MOsc] + list(d[1:4]), open(out, 'wb'))
