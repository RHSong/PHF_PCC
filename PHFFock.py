"""
Fock-matrix (diagonalisation based) PHF solver.

Uses FockTools_F (PHF_fort/FockTools.f90), which implements Eqs. 29-33 of
PHFB_formula.pdf for the PHF special case.  The one- and two-body integrals,
the projection operators R1/R2/Rpg and the orbitals MO are all kept in the
frozen-core (FC) basis; nothing is re-transformed between iterations.

Iteration:
    1. buildhsg  -> H, S, Fock integrand, rho integrand over the CI index
    2. EigenSolver(H, S) -> E0 and CI vector -> fsp, fpg, fk
    3. buildfock -> F = F0 + F+ + F-   (Eq. 29)
    4. SymmFilter(F)  (real / Sz-block-diagonal / RHF, mirroring PHF.py)
    5. transform F to the current MO basis, decouple the NFO frozen orbitals
       per spin (FrozenProj analogue), convergence test on dE and the active
       occupied-virtual block of F
    6. diagonalise the active block(s), occupy the lowest eigenvectors (aufbau)

Run as a script:  python PHFFock.py [guess.p] [eps]
    guess.p defaults to SGHF.p if present, else GHFMO.p.
    eps (optional) applies a random orthogonal rotation exp(eps*K) mixing
    occupied-alpha/virtual-beta and occupied-beta/virtual-alpha orbitals of the
    guess (SpinFlipPerturb), which lets the iteration leave the UHF manifold.
    The result is written to SGHFFock.p in the same layout as SGHF.p, and the
    eigenvalues of the final PHF Fock matrix to SGHFFock_moe.p.
"""
import numpy as np
import os
import sys
import time
import pickle
from parameter_fc import *
from PHFTools import EigenSolver, VecDecomp
from makeH import Mulliken2Dirac
from Spin import calcS, calcS2
from FockTools_F import focktools

NOccB = NOccSO - NOccA
NVrtB = NAO - NOccB


def BuildFockMatrix(HOne, HTwo, MO):
    """One evaluation of the PHF energy and Fock matrix at orbitals MO."""
    Hmat, Smat, Fockmat, Rdmmat = focktools.buildhsg(HOne, HTwo, MO, ngrid, nci, CmplxConj, ncisp, NOccSO,
                                                     R1, R2, Rpg, weightsp, weightpg, roota, rootb, rooty, J, SP,
                                                     npg, ncipg, NSO)
    evals, evecs = EigenSolver(Hmat, Smat)
    E0 = evals[0].real
    fsp, fpg, fk = VecDecomp(evecs[:, 0])
    F = focktools.buildfock(Fockmat, Rdmmat, Smat, E0, fsp, fpg, fk, CmplxConj, NSO, nci, ncisp, ncipg)
    return E0, F, fsp, fpg, fk, evals


def SymmFilter(F):
    """Impose on F the same constraints PHF.py imposes on the Thouless matrix Z."""
    F = np.array(F, dtype=complex)
    if (CmplxConj == 0):
        F = F.real.astype(complex)
    if (SP == 2 or is_RHF):
        # no alpha-beta mixing (SzProj)
        F[:NAO, NAO:] = 0
        F[NAO:, :NAO] = 0
    if (is_RHF):
        Favg = 0.5 * (F[:NAO, :NAO] + F[NAO:, NAO:])
        F[:NAO, :NAO] = Favg
        F[NAO:, NAO:] = Favg
    return F


NFO = 0
"""Number of frozen occupied orbitals per spin, as in PHF.py: the first NFO alpha
and the first NFO beta occupied columns of the guess (OA OB VA VB layout) are kept
fixed and never mix with the rest.  The guess should be RHF/UHF type."""


def FrozenIndex(nfo):
    """Column indices of the frozen orbitals in the OA OB VA VB layout."""
    return list(range(nfo)) + list(range(NOccA, NOccA + nfo))


def FockMO(F, MO, nfo, level_shift=0.0):
    """F in the current MO basis, with the frozen orbitals decoupled from everything
    else (their rows and columns are zeroed except the diagonal; this is the Fock
    analogue of FrozenProj in PHF.py) and an optional level shift on the virtuals."""
    Fmo = MO.T.conj() @ F @ MO
    frz = FrozenIndex(nfo)
    if (len(frz) > 0):
        d = np.diag(Fmo).copy()
        Fmo[frz, :] = 0
        Fmo[:, frz] = 0
        for i in frz:
            Fmo[i, i] = d[i]
    if (level_shift != 0.0):
        for v in range(NOccSO, NSO):
            Fmo[v, v] += level_shift
    return Fmo


def Diagonalize(Fmo, MO, nfo):
    """Diagonalise the active block(s) of Fmo (F in the current MO basis) and
    occupy the lowest eigenvectors (aufbau).  Frozen columns are copied over
    unchanged.  For SP=2 / is_RHF the alpha and beta blocks are diagonalised
    separately.  Returns (moe, newMO) in the OA OB VA VB column layout, where
    moe holds the eigenvalues and, for frozen orbitals, the diagonal element."""
    frz = FrozenIndex(nfo)
    newMO = np.zeros([NSO, NSO], dtype=complex)
    moe = np.zeros(NSO)
    for i in frz:
        newMO[:, i] = MO[:, i]
        moe[i] = Fmo[i, i].real
    if (SP == 2 or is_RHF):
        blocks = [(range(0, NOccA), range(NOccSO, NOccSO + NVrtA)),
                  (range(NOccA, NOccSO), range(NOccSO + NVrtA, NSO))]
    else:
        blocks = [(range(0, NOccSO), range(NOccSO, NSO))]
    for occ_cols, vir_cols in blocks:
        occ_act = [i for i in occ_cols if i not in frz]
        vir_cols = list(vir_cols)
        act = occ_act + vir_cols
        e, U = np.linalg.eigh(Fmo[np.ix_(act, act)])
        C = MO[:, act] @ U
        n = len(occ_act)
        newMO[:, occ_act] = C[:, :n]
        moe[occ_act] = e[:n]
        newMO[:, vir_cols] = C[:, n:]
        moe[vir_cols] = e[n:]
    return moe, newMO


def Residual(Fmo):
    """[F, rho] in the MO basis: the active occupied-virtual block of Fmo."""
    P = np.zeros([NSO, NSO])
    P[:NOccSO, :NOccSO] = np.eye(NOccSO)
    return Fmo @ P - P @ Fmo


def SpinFlipPerturb(MO, eps, seed=0):
    """Random orthogonal rotation restricted to the spin-flip occupied-virtual
    block: occupied-alpha with virtual-beta and occupied-beta with virtual-alpha.
    Assumes the OA OB VA VB column layout (UHF-type guess).  A UHF determinant
    is a stationary point of the SGHF energy along every other direction, so
    this is the block that has to be perturbed to leave the UHF manifold."""
    from scipy.linalg import expm
    rng = np.random.default_rng(seed)
    K = np.zeros([NSO, NSO])
    oa = slice(0, NOccA)
    ob = slice(NOccA, NOccSO)
    va = slice(NOccSO, NOccSO + NVrtA)
    vb = slice(NOccSO + NVrtA, NSO)
    K[oa, vb] = rng.standard_normal([NOccA, NVrtB])
    K[ob, va] = rng.standard_normal([NOccB, NVrtA])
    K = K - K.T
    K = K / np.max(np.abs(K))
    return MO @ expm(eps * K)


def Density(MO):
    occ = MO[:, :NOccSO]
    return occ @ occ.T.conj()


def DIIS(Flist, Errlist):
    """Pulay DIIS extrapolation of the Fock matrix.
    Flist: previous Fock matrices, Errlist: matching error vectors [F, rho]."""
    n = len(Flist)
    B = -np.ones([n + 1, n + 1], dtype=complex)
    B[n, n] = 0
    for i in range(n):
        for j in range(n):
            B[i, j] = np.vdot(Errlist[i].ravel(), Errlist[j].ravel())
    rhs = np.zeros(n + 1, dtype=complex)
    rhs[n] = -1
    c = np.linalg.solve(B, rhs)
    F = np.zeros_like(Flist[0])
    for i in range(n):
        F += c[i] * Flist[i]
    return F


def optPHFFock(HOne, HTwo, MO0, maxiter=100, tol_e=1e-6, tol_g=5e-5, level_shift=0.0,
               use_diis=True, ndiis=8, diis_start=1, nfo=None, log='PHFFockout'):
    """Self-consistent PHF via repeated diagonalisation of the PHF Fock matrix.

    DIIS (Pulay, error = [F, rho]) is applied from iteration diis_start on,
    keeping the last ndiis Fock matrices.  nfo frozen occupied orbitals per spin
    (default: module variable NFO) are kept fixed, see FrozenIndex/FockMO.
    After every diagonalisation the new orbitals are written to newmo.p
    ([MO, fsp, fpg, fk], same layout as SGHFFock.p), analogous to newsol.p; the orbitals with the
    lowest energy evaluated so far are written to newmo_best.p (with DIIS the energy need not
    decrease monotonically).
    Returns E0 (electronic, without Enuc), MO, fsp, fpg, fk, converged, moe,
    where moe are the eigenvalues ("MO energies") of the final PHF Fock matrix
    in the OA OB VA VB order used for MO."""
    if (CmplxConj == 2):
        sys.exit("PHFFock: the complex-conjugation projection (CmplxConj=2) has not been derived for the Fock solver")
    if (nfo is None):
        nfo = NFO
    MO = np.array(MO0, dtype=complex)
    if (CmplxConj == 0):
        MO = MO.real.astype(complex)
    Eold = 0.0
    converged = False
    Flist, Errlist = [], []
    output = open(log, 'a')
    if (nfo > 0):
        print("frozen orbitals (per spin):", nfo, " columns", FrozenIndex(nfo), file=output)
    moe = None
    Ebest = None
    for it in range(maxiter):
        t1 = time.time()
        E0, F, fsp, fpg, fk, evals = BuildFockMatrix(HOne, HTwo, MO)
        F = SymmFilter(F)
        Fmo = FockMO(F, MO, nfo)
        err_mo = Residual(Fmo)
        gmax = np.max(np.abs(err_mo))
        dE = E0 - Eold
        t2 = time.time()
        msg = "iter %3d  E= %.10f  dE= % .3e  max|[F,rho]|= %.3e  time= %.1f" % (it, E0 + Enuc, dE, gmax, t2 - t1)
        print(msg, flush=True)
        print(msg, file=output, flush=True)
        if (Ebest is None or E0 < Ebest):
            # orbitals whose energy was just evaluated (newmo.p below holds the NEXT, not yet evaluated, iterate)
            Ebest = E0
            pickle.dump([MO, fsp, fpg, fk], open("newmo_best.p", "wb"))
        moe = Diagonalize(Fmo, MO, nfo)[0]
        if (it > 0 and abs(dE) < tol_e and gmax < tol_g):
            converged = True
            break
        if (use_diis and it >= diis_start):
            Flist.append(F)
            Errlist.append(MO @ err_mo @ MO.T.conj())
            if (len(Flist) > ndiis):
                Flist.pop(0)
                Errlist.pop(0)
            if (len(Flist) > 1):
                F = DIIS(Flist, Errlist)
                Fmo = FockMO(F, MO, nfo, level_shift)
        if (level_shift != 0.0 and not (use_diis and it >= diis_start and len(Flist) > 1)):
            Fmo = FockMO(F, MO, nfo, level_shift)
        e, MO = Diagonalize(Fmo, MO, nfo)
        pickle.dump([MO, fsp, fpg, fk], open("newmo.p", "wb"))   # like newsol.p in PHF.py: latest orbitals, overwritten every iteration
        Eold = E0
    print("MO energies of the final PHF Fock matrix (occupied | virtual):", file=output)
    print(np.array2string(moe[:NOccSO].real, precision=6, max_line_width=120), file=output)
    print(np.array2string(moe[NOccSO:].real, precision=6, max_line_width=120), file=output)
    output.close()
    if (not converged):
        print("PHFFock: not converged after", maxiter, "iterations")
    return E0, MO, fsp, fpg, fk, converged, moe


if __name__ == "__main__":
    t1 = time.time()
    if (len(sys.argv) > 1):
        guess = sys.argv[1]
    elif (os.path.exists("SGHF.p")):
        guess = "SGHF.p"
    else:
        guess = "GHFMO.p"
    data = pickle.load(open(guess, "rb"))
    if (isinstance(data, (list, tuple))):
        MO0 = data[0]
    else:
        MO0 = data
    print("initial guess from", guess, " NFO =", NFO)
    if (len(sys.argv) > 2):
        eps = float(sys.argv[2])
        MO0 = SpinFlipPerturb(MO0, eps)
        print("guess perturbed in the spin-flip ov block with eps =", eps)

    HOne = np.zeros([NSO, NSO])
    HTwo = np.zeros([NSO, NSO, NSO, NSO])
    HOne[:NAO, :NAO] = h1[:, :]
    HOne[NAO:, NAO:] = h1[:, :]
    HTwo[:NAO, :NAO, :NAO, :NAO] = eri[:, :, :, :]
    HTwo[NAO:, NAO:, NAO:, NAO:] = eri[:, :, :, :]
    HTwo[NAO:, NAO:, :NAO, :NAO] = eri[:, :, :, :]
    HTwo[:NAO, :NAO, NAO:, NAO:] = eri[:, :, :, :]
    HTwo = Mulliken2Dirac(HTwo)

    EPHF, MO, fsp, fpg, fk, conv, moe = optPHFFock(HOne, HTwo, MO0)
    S = calcS(MO, NOccSO, NAO)
    SS = calcS2(MO, NOccSO, NAO)
    t2 = time.time()
    print("E(PHF-Fock)=", EPHF + Enuc, " converged:", conv)
    print("det S^2, Sxyz=", SS, S)
    print("PHFFock time=", t2 - t1)
    print("MO energies (occupied):", np.array2string(moe[:NOccSO].real, precision=6, max_line_width=120))
    print("MO energies (virtual): ", np.array2string(moe[NOccSO:].real, precision=6, max_line_width=120))
    # [MO, fsp, fpg, fk, moe]: moe = eigenvalues of the PHF effective Fock (oo/vv blocks) in the
    # order of the MO columns; readers that expect 4 entries should unpack data[:4]
    pickle.dump([MO, fsp, fpg, fk, np.real(moe)], open("SGHFFock.p", "wb"))
    pickle.dump(np.real(moe), open("SGHFFock_moe.p", "wb"))
