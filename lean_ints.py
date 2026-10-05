"""Memory-lean spin-orbital integrals for the PCC driver (Main_PCC.py).

The old path built the real spin-orbital Mulliken tensor (NSO^4), a Dirac copy,
then ao2mo by einsum with complex NSO^4 intermediates and a C-ordered result that
f2py copied once more to Fortran order (peak ~4 NSO^4 complex arrays).

Here only ONE NSO^4 complex array is ever allocated: the final antisymmetrised
Dirac tensor <pq||rs>, in Fortran order, which f2py passes to PGCC without a copy.

  * the spatial AO integrals are unpacked with pyscf.ao2mo.restore (NAO^4 real);
  * pyscf's ao2mo transforms take real coefficients only, so the complex GHF
    transform is done as two GEMMs over pair densities
        P[(mu nu),(p q)] = sum_sigma conj(C_sigma[mu,p]) C_sigma[nu,q]
        (pq|rs)          = P^T . ERI . P
    with pyscf.lib.dot (threaded BLAS);
  * the second GEMM is done one MO index s at a time, and every slab
    (pq|r s) is antisymmetrised straight into H2[:,:,:,s]:
        <pq||rs> = (pr|qs) - (ps|qr) = A[p,r,q] - A[q,r,p],  A[a,b,c] = (ab|c s)
    using (ps|qr) = (qr|ps) (real AO integrals).
Peak memory (inside the slab loop): H2 (NSO^4 complex) + 2 * NAO^2 NSO^2 complex + NAO^4 real.
For Cr2 / NFC=10 (NAO=76): 8.5 + 2*2.1 + 0.3 = 13 GB.
"""
import numpy as np
from pyscf import ao2mo as pyscf_ao2mo
from pyscf import lib


def unpack_eri(eri, NAO):
    """Full 4-index real (ij|kl) array from any pyscf storage (s1, s4, s8)."""
    eri = np.asarray(eri)
    if eri.ndim == 4:
        return eri
    return pyscf_ao2mo.restore(1, eri, NAO)


def ghf_fock_ao(h1, eri4, MOocc, NAO):
    """Spin-orbital AO Fock, identical to makeH.MakeFock(HOne, Mulliken2Dirac(HTwo), MOocc)
    but built from the spatial integrals:  F_ik = h_ik + sum_jl <ij||kl> rho_lj."""
    NSO = 2 * NAO
    rho = MOocc @ MOocc.T.conj()
    a, b = slice(0, NAO), slice(NAO, NSO)
    J = lambda d: np.einsum('ikjl,lj->ik', eri4, d)      # sum (ik|jl) d_lj
    K = lambda d: np.einsum('iljk,lj->ik', eri4, d)      # sum (il|jk) d_lj
    F = np.zeros([NSO, NSO], dtype=complex)
    Jt = J(rho[a, a] + rho[b, b])
    F[a, a] = h1 + Jt - K(rho[a, a])
    F[b, b] = h1 + Jt - K(rho[b, b])
    F[a, b] = -K(rho[a, b])
    F[b, a] = -K(rho[b, a])
    return F


def semicanon_from_fock(F, MO, NOcc, SP):
    """makeH.SemiCanon with the AO Fock supplied (no NSO^4 tensor needed)."""
    NOccA = NOcc // 2
    NVrtA = (len(F) - NOcc) // 2
    F = MO.T.conj() @ F @ MO
    Foo = F[:NOcc, :NOcc]
    Fvv = F[NOcc:, NOcc:]
    if SP == 2:
        ea, Ua = np.linalg.eigh(Foo[:NOccA, :NOccA])
        eb, Ub = np.linalg.eigh(Foo[NOccA:, NOccA:])
        mo_o = np.hstack((MO[:, :NOccA] @ Ua, MO[:, NOccA:NOcc] @ Ub))
        ea, Ua = np.linalg.eigh(Fvv[:NVrtA, :NVrtA])
        eb, Ub = np.linalg.eigh(Fvv[NVrtA:, NVrtA:])
        mo_v = np.hstack((MO[:, NOcc:NOcc + NVrtA] @ Ua, MO[:, NOcc + NVrtA:] @ Ub))
    else:
        e, U = np.linalg.eigh(Foo)
        mo_o = MO[:, :NOcc] @ U
        e, U = np.linalg.eigh(Fvv)
        mo_v = MO[:, NOcc:] @ U
    return np.hstack((mo_o, mo_v))


def dirac_mo_eri_F(eri4, MO, NAO):
    """H2[p,q,r,s] = <pq||rs> in the (complex) spin-orbital MO basis MO (NSO x NSO,
    rows = alpha AOs then beta AOs).  Returns complex128, Fortran-ordered."""
    NSO = MO.shape[1]
    n2 = NAO * NAO
    Ca, Cb = MO[:NAO], MO[NAO:]
    # pair density P[mu,nu,p,q], flattened F-order: row mu+NAO*nu, col p+NSO*q
    P = np.einsum('mp,nq->mnpq', Ca.conj(), Ca)
    P += np.einsum('mp,nq->mnpq', Cb.conj(), Cb)
    P = P.reshape(n2, NSO * NSO, order='F')
    E = np.asarray(eri4).reshape(n2, n2, order='F')      # row mu+NAO*nu, col lam+NAO*kap (symmetric)
    PT = np.ascontiguousarray(P.T)                       # (NSO^2, NAO^2): row p+NSO*q
    del P
    # Yt[(r s),(mu nu)] = sum_(lam kap) P[(lam kap),(r s)] (lam kap|mu nu)  (= (E P)^T, E symmetric);
    # rows r+NSO*s with fixed s are contiguous, so each slab is a plain GEMM operand
    Yt = lib.dot(PT, E)
    H2 = np.empty([NSO, NSO, NSO, NSO], dtype=complex, order='F')
    for s in range(NSO):
        # A[p,q,r] = (pq|rs) for this s
        A = lib.dot(PT, Yt[s * NSO:(s + 1) * NSO].T).reshape(NSO, NSO, NSO, order='F')
        H2[:, :, :, s] = np.transpose(A, (0, 2, 1)) - np.transpose(A, (2, 0, 1))
    del Yt, PT
    return H2


def check_against_dense(h1, eri, MO, NAO, NOcc, SP):
    """Small-system check vs. the original dense path (makeH + ao2mo). Not for Cr2."""
    from makeH import Mulliken2Dirac, SemiCanon
    from ao2mo import ao2mo
    NSO = 2 * NAO
    eri4 = unpack_eri(eri, NAO)
    HOne = np.zeros([NSO, NSO]); HOne[:NAO, :NAO] = h1; HOne[NAO:, NAO:] = h1
    HTwo = np.zeros([NSO] * 4)
    HTwo[:NAO, :NAO, :NAO, :NAO] = eri4; HTwo[NAO:, NAO:, NAO:, NAO:] = eri4
    HTwo[NAO:, NAO:, :NAO, :NAO] = eri4; HTwo[:NAO, :NAO, NAO:, NAO:] = eri4
    HTwo = Mulliken2Dirac(HTwo)
    MO = np.array(MO, dtype=complex)
    mo_ref = SemiCanon(HOne, HTwo, MO, NOcc, SP)
    mo_new = semicanon_from_fock(ghf_fock_ao(h1, eri4, MO[:, :NOcc], NAO), MO, NOcc, SP)
    H2_ref = ao2mo(HTwo, mo_ref, 4)
    H2_new = dirac_mo_eri_F(eri4, mo_ref, NAO)
    return np.max(np.abs(mo_ref - mo_new)), np.max(np.abs(H2_ref - H2_new))
