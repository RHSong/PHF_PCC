"""Spin-orbital frozen core for PHF / PCC: freeze a spin-polarized (UHF) core.

The standard frozen core (FrozenCore.py) contracts a closed-shell core of spatial orbitals and leaves spatial
integrals (h1, eri) that the drivers expand as 'spatial x spin'.  Here the core may differ for alpha and beta
(e.g. the 1s2s2p of an AFM UHF), so the active one-electron operator becomes spin dependent and the active
two-electron integrals are spin-block dependent:

    basis:   spin orbitals  p^alpha = C_a[:, p] alpha,  p^beta = C_b[:, p] beta   (p over the NAO active orbitals;
             alpha set orthonormal, beta set orthonormal, alpha _|_ beta by spin; C_a and C_b are the active
             columns of the alpha and beta UHF orbitals)
    h_a, h_b:          active one-electron integrals including the core Coulomb/exchange potential
    eri_aa, eri_ab, eri_bb:  Mulliken (pq|rs) with the first pair in the alpha/alpha/beta set and the second pair
                       in the alpha/beta/beta set  ((pq|rs)_ba = (rs|pq)_ab)
    O_ab[p, q] = <p^alpha | q^beta>_spatial = (C_a^T S C_b)[p, q]

Spin rotations R(Omega) = r(Omega) (x) 1 in the 'spatial x spin' basis become, in this basis,
    R'[sigma, tau] = r_{sigma tau} O_{sigma tau},   O_aa = O_bb = 1, O_ab as above, O_ba = O_ab^T,
truncated to the active space (the neglected core-active mixing is of the order 1 - sigma_min(O_core_ab)).
Point-group operators are transformed the same way when the AO representation is available.

FC_SO.p is a dict with 'tag' = 'FC_SO' (see build_fc_so); the drivers (run_pcc.py, run_pcc_t.py, PHFFock.py via
parameter_fc.py) detect it and use the routines below instead of the 'spatial x spin' expansion.
"""
import numpy as np
from pyscf import ao2mo as pyscf_ao2mo
from pyscf import lib

TAG = 'FC_SO'


def _mo_eri(eri_ao, C1, C2, C3, C4):
    """Mulliken (pq|rs) over real MO coefficient sets, full 4-index array."""
    nao = C1.shape[0]
    eri8 = pyscf_ao2mo.restore(8, np.asarray(eri_ao), nao) if np.asarray(eri_ao).ndim == 4 else eri_ao
    out = pyscf_ao2mo.general(eri8, (C1, C2, C3, C4), compact=False)
    return out.reshape(C1.shape[1], C2.shape[1], C3.shape[1], C4.shape[1])


def build_fc_so(h1_ao, eri_ao, S, Ca, Cb, NFC, Enuc, NOccSO, Spin2=0, NFV=0, Rpg_ao=None, ncipg=1, weightpg=None):
    """Spin-orbital frozen-core Hamiltonian from UHF orbitals Ca, Cb (AO x nmo, both orthonormal in the S metric).

    The first NFC columns of Ca (Cb) are the frozen alpha (beta) core; columns NFC:nmo-NFV are active.
    Returns the FC_SO dict:  Enuc (nuclear + core energy), h_a, h_b, eri_aa, eri_ab, eri_bb (active), O_ab,
    NAO (= number of active spatial functions per spin), NOccSO (active electrons), Spin2, ncipg, weightpg,
    Rpg (NSO x NSO, identity if no point group), C_a, C_b (AO coefficients of the active sets).
    Rpg_ao: optional list of point-group operators as AO matrices acting on the ket (R phi = phi R_ao);
    they are transformed as (C_sigma^T S R_ao C_tau) blocks."""
    nmo = Ca.shape[1]; act = slice(NFC, nmo - NFV); core = slice(0, NFC)
    h_a = Ca.T @ h1_ao @ Ca; h_b = Cb.T @ h1_ao @ Cb
    eaa = _mo_eri(eri_ao, Ca, Ca, Ca, Ca); ebb = _mo_eri(eri_ao, Cb, Cb, Cb, Cb); eab = _mo_eri(eri_ao, Ca, Ca, Cb, Cb)
    ia = np.arange(NFC)
    # core energy: sum_i h_ii + 1/2 sum_ij [(ii|jj) - (ij|ji) delta_spin]
    E0 = np.trace(h_a[core, core]) + np.trace(h_b[core, core])
    E0 += 0.5 * (np.einsum('iijj', eaa[core, core, core, core]) - np.einsum('ijji', eaa[core, core, core, core])
                 + np.einsum('iijj', ebb[core, core, core, core]) - np.einsum('ijji', ebb[core, core, core, core])
                 + 2.0 * np.einsum('iijj', eab[core, core, core, core]))
    # active one-electron operator with the core potential
    ha = h_a[act, act] + np.einsum('pqii->pq', eaa[act, act, core, core]) - np.einsum('piiq->pq', eaa[act, core, core, act]) \
        + np.einsum('pqjj->pq', eab[act, act, core, core])
    hb = h_b[act, act] + np.einsum('pqii->pq', ebb[act, act, core, core]) - np.einsum('piiq->pq', ebb[act, core, core, act]) \
        + np.einsum('jjpq->pq', eab[core, core, act, act])
    Ca_act, Cb_act = Ca[:, act], Cb[:, act]
    NAO = Ca_act.shape[1]; NSO = 2 * NAO
    O_ab = Ca_act.T @ S @ Cb_act
    if Rpg_ao is None:
        Rpg = np.zeros([1, NSO, NSO], dtype=complex); Rpg[0] = np.eye(NSO); ncipg = 1; weightpg = np.ones([1, 1, 1])
    else:
        Rpg = np.zeros([len(Rpg_ao), NSO, NSO], dtype=complex)
        for g, R in enumerate(Rpg_ao):
            Rpg[g, :NAO, :NAO] = Ca_act.T @ S @ R @ Ca_act
            Rpg[g, NAO:, NAO:] = Cb_act.T @ S @ R @ Cb_act
    core_ov = np.linalg.svd(Ca[:, core].T @ S @ Cb[:, core], compute_uv=False)
    return {'tag': TAG, 'Enuc': Enuc + E0, 'E0_core': E0, 'h_a': ha, 'h_b': hb,
            'eri_aa': np.ascontiguousarray(eaa[act, act, act, act]), 'eri_ab': np.ascontiguousarray(eab[act, act, act, act]),
            'eri_bb': np.ascontiguousarray(ebb[act, act, act, act]), 'O_ab': O_ab, 'NAO': NAO, 'NOccSO': NOccSO - 2 * NFC,
            'Spin2': Spin2, 'ncipg': ncipg, 'weightpg': weightpg, 'Rpg': Rpg, 'C_a': Ca_act, 'C_b': Cb_act, 'NFC': NFC,
            'core_ab_overlap_sv': core_ov}


def is_fc_so(fc):
    return isinstance(fc, dict) and fc.get('tag') == TAG


def uhf_det_so(NAO, NOccA, NOccB):
    """The UHF determinant in its own spin-orbital basis, OA OB VA VB layout (unit columns)."""
    NSO = 2 * NAO
    cols = list(range(NOccA)) + [NAO + i for i in range(NOccB)] + list(range(NOccA, NAO)) + [NAO + i for i in range(NOccB, NAO)]
    return np.eye(NSO)[:, cols]


def hone_so(h_a, h_b):
    NAO = h_a.shape[0]; NSO = 2 * NAO
    H = np.zeros([NSO, NSO]); H[:NAO, :NAO] = h_a; H[NAO:, NAO:] = h_b
    return H


def htwo_mulliken_so(eri_aa, eri_ab, eri_bb):
    """Dense NSO^4 Mulliken tensor (pq|rs) in the spin-orbital basis (for the Fock-PHF driver)."""
    NAO = eri_aa.shape[0]; NSO = 2 * NAO; a, b = slice(0, NAO), slice(NAO, NSO)
    H = np.zeros([NSO] * 4)
    H[a, a, a, a] = eri_aa; H[b, b, b, b] = eri_bb; H[a, a, b, b] = eri_ab; H[b, b, a, a] = np.transpose(eri_ab, (2, 3, 0, 1))
    return H


def ghf_fock_so(h_a, h_b, eri_aa, eri_ab, eri_bb, MOocc, NAO):
    """Spin-orbital Fock matrix F_ik = h_ik + sum_jl <ij||kl> rho_lj in the spin-orbital basis (cf. lean_ints.ghf_fock_ao)."""
    NSO = 2 * NAO; a, b = slice(0, NAO), slice(NAO, NSO)
    rho = MOocc @ MOocc.T.conj()
    J = lambda e, d: np.einsum('ikjl,lj->ik', e, d)      # sum (ik|jl) d_lj, first pair = row spin, second pair = spin of d
    K = lambda e, d: np.einsum('iljk,lj->ik', e, d)      # sum (il|jk) d_lj
    eba = np.transpose(eri_ab, (2, 3, 0, 1))
    F = np.zeros([NSO, NSO], dtype=complex)
    F[a, a] = h_a + J(eri_aa, rho[a, a]) + J(eri_ab, rho[b, b]) - K(eri_aa, rho[a, a])
    F[b, b] = h_b + J(eba, rho[a, a]) + J(eri_bb, rho[b, b]) - K(eri_bb, rho[b, b])
    F[a, b] = -K(eri_ab, rho[a, b])     # (i^a l^a | j^b k^b) rho[l^a, j^b]... see below
    F[b, a] = -K(eba, rho[b, a])
    return F


def dirac_mo_eri_so(eri_aa, eri_ab, eri_bb, MO, NAO):
    """H2[p,q,r,s] = <pq||rs> in the (complex) MO basis MO (NSO x NSO, rows = alpha set then beta set) from the
    spin-block Mulliken integrals; same slab algorithm and memory layout as lean_ints.dirac_mo_eri_F."""
    NSO = MO.shape[1]; n2 = NAO * NAO
    Ca, Cb = MO[:NAO], MO[NAO:]
    Pa = np.einsum('mp,nq->mnpq', Ca.conj(), Ca).reshape(n2, NSO * NSO, order='F')
    Pb = np.einsum('mp,nq->mnpq', Cb.conj(), Cb).reshape(n2, NSO * NSO, order='F')
    Eaa = np.asarray(eri_aa).reshape(n2, n2, order='F'); Ebb = np.asarray(eri_bb).reshape(n2, n2, order='F')
    Eab = np.asarray(eri_ab).reshape(n2, n2, order='F')          # rows: alpha pair, cols: beta pair
    PTa = np.ascontiguousarray(Pa.T); PTb = np.ascontiguousarray(Pb.T); del Pa, Pb
    # (PQ|RS) = Pa^T Eaa Pa + Pa^T Eab Pb + Pb^T Eab^T Pa + Pb^T Ebb Pb
    #        = PTa . Ya[:, RS] + PTb . Yb[:, RS],   Ya = Eaa Pa + Eab Pb,  Yb = Eab^T Pa + Ebb Pb
    # transposed (rows RS, cols pair):  Yta = PTa Eaa + PTb Eab^T,  Ytb = PTa Eab + PTb Ebb   (Eaa, Ebb symmetric)
    Yta = lib.dot(PTa, Eaa); Yta += lib.dot(PTb, np.ascontiguousarray(Eab.T))
    Ytb = lib.dot(PTa, Eab); Ytb += lib.dot(PTb, Ebb)
    H2 = np.empty([NSO, NSO, NSO, NSO], dtype=complex, order='F')
    for s in range(NSO):
        A = lib.dot(PTa, Yta[s * NSO:(s + 1) * NSO].T) + lib.dot(PTb, Ytb[s * NSO:(s + 1) * NSO].T)
        A = A.reshape(NSO, NSO, NSO, order='F')                   # A[p,q,r] = (pq|rs)
        H2[:, :, :, s] = np.transpose(A, (0, 2, 1)) - np.transpose(A, (2, 0, 1))
    del Yta, Ytb, PTa, PTb
    return H2


def spin_rot_so(R, O_ab, NAO):
    """Transform a product-form spin rotation R = r (x) 1 (NSO x NSO) into the spin-orbital basis with alpha/beta
    overlap O_ab (truncated to the active space)."""
    NSO = 2 * NAO; a, b = slice(0, NAO), slice(NAO, NSO)
    r = np.array([[R[0, 0], R[0, NAO]], [R[NAO, 0], R[NAO, NAO]]])
    for (s, t), blk in (((0, 0), R[a, a]), ((0, 1), R[a, b]), ((1, 0), R[b, a]), ((1, 1), R[b, b])):
        assert np.allclose(blk, r[s, t] * np.eye(NAO), atol=1e-12), "spin_rot_so: R is not of product form r (x) 1"
    Rn = np.zeros([NSO, NSO], dtype=complex)
    Rn[a, a] = r[0, 0] * np.eye(NAO); Rn[b, b] = r[1, 1] * np.eye(NAO)
    Rn[a, b] = r[0, 1] * O_ab; Rn[b, a] = r[1, 0] * O_ab.T
    return Rn


def spin_rot_list_so(Rlist, O_ab, NAO):
    out = np.array(Rlist, dtype=complex)
    for i in range(len(out)): out[i] = spin_rot_so(out[i], O_ab, NAO)
    return out


def spin_so(MO, NOcc, Ca, Cb):
    """<S^2> and <S_x,y,z> of a determinant (columns of MO in the spin-orbital basis, first NOcc occupied) when the
    alpha and beta basis sets are different spatial functions Ca[:, p] and Cb[:, p] (orthonormal metric assumed for the
    columns of Ca, Cb, e.g. the RHF-basis coefficients stored in FC_SO.p).  Replaces Spin.calcS / calcS2, which assume
    p^alpha and p^beta share the spatial function."""
    from pyscf.scf import ghf
    NAO = Ca.shape[1]; occ = np.asarray(MO)[:, :NOcc]
    Pa = Ca @ occ[:NAO]; Pb = Cb @ occ[NAO:]                      # spatial parts of the occupied spin orbitals
    n = Ca.shape[0]
    ss, mult = ghf.spin_square(np.vstack([Pa, Pb]), s=np.eye(n))
    Sz = 0.5 * (np.trace(Pa.conj().T @ Pa) - np.trace(Pb.conj().T @ Pb)).real
    Sp = np.trace(Pa.conj().T @ Pb)                               # <S+> = sum_i <alpha part_i | beta part_i>
    return ss, (Sp.real, Sp.imag, Sz)


def match_orb_blocks(C, Ref, blocks):
    """MatchOrb.py applied block-wise: rotate the columns of C inside each block (list of (start, stop)) by the
    orthogonal Procrustes solution of the SVD of Ref_block^T C_block, so that the new orbitals have maximal overlap
    with the reference ones (orbital following between geometries / runs).  Plain coefficient overlaps, as MatchOrb."""
    C = np.array(C, copy=True)
    for a, b in blocks:
        m = Ref[:, a:b].T.conj() @ C[:, a:b]
        u, sv, vh = np.linalg.svd(m)
        C[:, a:b] = C[:, a:b] @ vh.conj().T @ u.conj().T
        print('  block %3d:%3d  overlap singular values min %.5f max %.5f' % (a, b, sv.min(), sv.max()))
    return C


def load_fc(path):
    """Uniform view of FC.p (standard) or FC_SO.p: returns dict with keys
    Enuc, h_a, h_b, eri_aa, eri_ab, eri_bb (full 4-index, may alias), O_ab, ncipg, weightpg, Rpg, NAO, NOccSO, Spin2, so."""
    import pickle
    fc = pickle.load(open(path, 'rb'))
    if is_fc_so(fc):
        d = dict(fc); d['so'] = True; return d
    Enuc, h1, eri, ncipg, weightpg, Rpg, NAO, NOccSO = fc[:8]
    NAO = int(NAO)
    eri4 = pyscf_ao2mo.restore(1, np.asarray(eri), NAO) if np.asarray(eri).ndim != 4 else np.asarray(eri)
    return {'Enuc': Enuc, 'h_a': h1, 'h_b': h1, 'eri_aa': eri4, 'eri_ab': eri4, 'eri_bb': eri4, 'O_ab': np.eye(NAO),
            'ncipg': ncipg, 'weightpg': weightpg, 'Rpg': Rpg, 'NAO': NAO, 'NOccSO': int(NOccSO),
            'Spin2': int(fc[8]) if len(fc) > 8 else 0, 'so': False}


if __name__ == "__main__":
    """PHF_UFC step 3: FC.p (full RHF basis, Main_FC_full.py) + MOL_MO.p (Main_UHF.py) + NFC (parameter.py) ->
         usage: python ufc_core.py [MOL_MO.p | -] [REF_MOL_MO.p]   (REF: previous geometry; its UHF sets define the orbital order/phases)
         U_SO.p   : U = C_rhf^T S C_uhf per spin, as one NSO_full x NSO_full matrix, columns OA OB VA VB
                    (core alpha first in OA, core beta first in OB)
         FC_SO.p  : spin-orbital frozen-core Hamiltonian in the active UHF spin-orbital basis (dict, tag FC_SO)
         GHFMO.p  : the UHF determinant in that active basis (unit columns, OA OB VA VB) = start for PHF."""
    import pickle, sys
    from parameter import NFC, NFV
    fc = pickle.load(open('FC.p', 'rb'))
    Enuc, h1, eri, ncipg, weightpg, Rpg, NAOf, NOccSOf = fc[:8]
    NAOf = int(NAOf); NOccSOf = int(NOccSOf); Spin2 = int(fc[8]) if len(fc) > 8 else 0
    m = pickle.load(open(sys.argv[1] if len(sys.argv) > 1 and sys.argv[1] != '-' else 'MOL_MO.p', 'rb'))   # argv[1]: MOL_MO file ('-' = default); argv[2]: reference MOL_MO.p to match (orbital following)
    S = m['S']; C_rhf = m['C_rhf']; na_tot, nb_tot = m['nelec']
    import os
    if os.path.exists('FCMO.p'):
        # FC.p holds h1/eri in the RHF basis of the Main_FC_full.py run that wrote FCMO.p; U must use exactly that basis.
        # (MOL_MO.p's C_rhf can differ by signs / degenerate rotations when FC.p is regenerated in another run.)
        C_fc = np.asarray(pickle.load(open('FCMO.p', 'rb')))
        print('RHF basis of FC.p (FCMO.p) vs MOL_MO.p C_rhf: max|diff| %.2e  -> using FCMO.p' % np.max(np.abs(C_fc - C_rhf)))
        C_rhf = C_fc
    if len(sys.argv) > 2:
        # orbital following (MatchOrb.py): align the UHF alpha / beta sets to those of a reference run (MOL_MO.p of the
        # previous geometry), block-wise so that the frozen core remains the span of the NFC lowest orbitals
        ref = pickle.load(open(sys.argv[2], 'rb')); nmo = m['C_a'].shape[1]
        print('matching the UHF orbitals to', sys.argv[2])
        m['C_a'] = match_orb_blocks(m['C_a'], ref['C_a'], [(0, NFC), (NFC, na_tot), (na_tot, nmo)])
        m['C_b'] = match_orb_blocks(m['C_b'], ref['C_b'], [(0, NFC), (NFC, nb_tot), (nb_tot, nmo)])
        pickle.dump(m, open('MOL_MO.p', 'wb'))
    Ua = C_rhf.T @ S @ m['C_a']; Ub = C_rhf.T @ S @ m['C_b']          # UHF orbitals in the (orthonormal) RHF MO basis
    print('U unitarity: max|U^T U - 1| = %.1e / %.1e' % (abs(Ua.T @ Ua - np.eye(NAOf)).max(), abs(Ub.T @ Ub - np.eye(NAOf)).max()))
    # sorted full transformation (for reference): OA = core-a + active occ a, OB = core-b + active occ b, VA, VB
    NSOf = 2 * NAOf
    Ufull = np.zeros([NSOf, NSOf]); Ufull[:NAOf, :na_tot] = Ua[:, :na_tot]; Ufull[NAOf:, na_tot:na_tot + nb_tot] = Ub[:, :nb_tot]
    Ufull[:NAOf, na_tot + nb_tot:na_tot + nb_tot + NAOf - na_tot] = Ua[:, na_tot:]; Ufull[NAOf:, na_tot + nb_tot + NAOf - na_tot:] = Ub[:, nb_tot:]
    pickle.dump(Ufull, open('U_SO.p', 'wb'))
    # spin-orbital frozen core: first NFC columns (lowest energies) of each spin; integrals of FC.p are in the RHF basis (S = 1)
    Rpg_sp = [np.asarray(R)[:NAOf, :NAOf] for R in Rpg] if (int(ncipg) > 1 or len(Rpg) > 1) else None
    d = build_fc_so(h1, eri, np.eye(NAOf), Ua, Ub, NFC, Enuc, NOccSOf, Spin2, NFV, Rpg_ao=Rpg_sp, ncipg=ncipg, weightpg=weightpg)
    pickle.dump(d, open('FC_SO.p', 'wb'), protocol=4)
    NAO = d['NAO']; NOccA, NOccB = na_tot - NFC, nb_tot - NFC
    pickle.dump(uhf_det_so(NAO, NOccA, NOccB), open('GHFMO.p', 'wb'))
    print('FC_SO.p: NFC %d  active NAO %d  NOccSO %d (%d a, %d b)  core energy %.10f  Enuc+E0 %.10f' % (NFC, NAO, d['NOccSO'], NOccA, NOccB, d['E0_core'], d['Enuc']))
    print('core alpha/beta overlap singular values: min %.6f ; active O_ab singular values: min %.6f' % (d['core_ab_overlap_sv'].min(), np.linalg.svd(d['O_ab'], compute_uv=False).min()))
    # check: the UHF determinant in this model must reproduce the all-electron UHF energy
    from makeH import Mulliken2Dirac, MakeEne
    HOne = hone_so(d['h_a'], d['h_b']); HTwo = Mulliken2Dirac(htwo_mulliken_so(d['eri_aa'], d['eri_ab'], d['eri_bb']))
    E = MakeEne(HOne, HTwo, uhf_det_so(NAO, NOccA, NOccB)[:, :NOccA + NOccB]).real + d['Enuc']
    print('check: E(UHF det in FC_SO model) = %.10f  vs all-electron UHF %.10f  diff %.2e' % (E, m['E_UHF'], abs(E - m['E_UHF'])))
    if abs(E - m['E_UHF']) > 1e-8:
        os.remove('FC_SO.p'); os.remove('GHFMO.p')
        sys.exit('ufc_core: FC_SO.p inconsistent with MOL_MO.p (UHF energy check %.2e > 1e-8): FC_SO.p removed' % abs(E - m['E_UHF']))
