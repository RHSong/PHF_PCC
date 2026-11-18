"""(T) correction for PCC from saved PCCSD amplitudes (no PCCSD rerun).

Coherent (T): the projected (T) of the PCCSD state, in the CC orbital basis (usually with the PHF effective-Fock
orbital energies, --moe eff). +(T): single-determinant (T) (--sp 0 --grid 1 1) on the amplitudes rotated into the
semicanonical orbitals of the bare Fock (semicanon_plusT.py, --newmo) with --moe fock.

Usage
-----
    coherent (T):
    python run_pcc_t.py --fc FC.p --mo pcc/MO_used.p --t1t2 pcc/T1T2.p --newmo SGHFFock.p --moe eff --sp 2 --grid 1 8
    +(T):
    python run_pcc_t.py --fc FC.p --mo pcc/MO_used.p --t1t2 pcc/T1T2.p --newmo MO_sc.p --moe fock --sp 0 --grid 1 1

Inputs
------
--mo      the exact orbital basis in which T1/T2 were converged (MO_used.p of run_pcc.py)
--newmo   optional second basis spanning the same occupied / virtual spaces; T1/T2 are rotated into it
--moe     'fock': denominators from the diagonal of the bare Fock matrix in the final basis
          'eff' : the PHF effective-Fock energies, 5th entry of the --newmo (or --mo) pickle [MO, fsp, fpg, fk, moe]
          file.p: pickled array of NSO orbital energies in the column order of the final basis
--pg 1    also apply the point-group projector stored in FC.p
--inc     also evaluate the old per-grid (incoherent) (T) sum
--epcc    PCCSD energy, to add the total energies to the result
--zero-overlap TAU
          zero-overlap-safe coherent (T) (pgcc_t_safe), deformed Thouless step below TAU

Output
------
result_T<tag>.json in the folder of --t1t2.
"""
import argparse
import json
import os
import pickle
import sys

import numpy as np

SUCCSD_DIR = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, SUCCSD_DIR)
from run_pcc import FrozenCoreData, load_phf_solution, get_spin_projection, get_pg_projection
from ao2mo import ao2mo
from lean_ints import unpack_eri, ghf_fock_ao, dirac_mo_eri_F


def parse_args():
    """Command-line options (see the module docstring); file paths are made absolute."""
    parser = argparse.ArgumentParser(description="(T) correction for PCC from saved amplitudes")
    parser.add_argument("--module", default="PUCC_X2", help="f2py module name")
    parser.add_argument("--libdir", default=SUCCSD_DIR, help="folder that holds the module")
    parser.add_argument("--fc", required=True, help="FC.p")
    parser.add_argument("--mo", required=True, help="orbital basis of T1/T2 (MO_used.p)")
    parser.add_argument("--t1t2", required=True, help="T1T2.p")
    parser.add_argument("--newmo", default=None, help="basis to rotate T1/T2 into")
    parser.add_argument("--moe", default="fock", help="fock | eff | pickle of NSO orbital energies")
    parser.add_argument("--sp", type=int, default=2, help="spin projection: 0 none, 1 SGHF, 2 SUHF, 3 Sz")
    parser.add_argument("--grid", type=int, nargs=2, default=[1, 8], help="spin-projection grid")
    parser.add_argument("--J", type=float, default=None, help="total spin (default from FC.p)")
    parser.add_argument("--pg", type=int, default=0, help="1: point-group projection from FC.p")
    parser.add_argument("--inc", action="store_true", help="also the incoherent per-grid (T)")
    parser.add_argument("--epcc", type=float, default=None, help="PCCSD energy")
    parser.add_argument("--tag", default="", help="suffix of the result file")
    parser.add_argument("--zero-overlap", type=float, default=None, metavar="TAU",
                        help="zero-overlap-safe coherent (T) (pgcc_t_safe) with threshold TAU")
    args = parser.parse_args()
    args.fc = os.path.abspath(args.fc)
    args.mo = os.path.abspath(args.mo)
    args.t1t2 = os.path.abspath(args.t1t2)
    args.libdir = os.path.abspath(args.libdir)
    if args.newmo:
        args.newmo = os.path.abspath(args.newmo)
    if args.moe not in ("fock", "eff"):
        args.moe = os.path.abspath(args.moe)
    return args


def load_effective_fock_energies(path):
    """PHF effective-Fock orbital energies: the 5th entry of [MO, fsp, fpg, fk, moe] (PHFFock.py); None if absent."""
    data = pickle.load(open(path, "rb"))
    if isinstance(data, (list, tuple)) and len(data) >= 5:
        return np.real(data[4])
    return None


def rotate_amplitudes(t1, t2, mo_old, mo_new, nocc):
    """Rotate T1/T2 from the orbitals mo_old into mo_new (same occupied and virtual spaces).

    With rot_occ = C_old,occ^H C_new,occ (new i' = sum_i old i rot_occ[i, i']) and rot_vir likewise,
        T1'[A, I] = sum_ai conj(rot_vir[a, A]) T1[a, i] rot_occ[i, I],   T2' likewise on each index.

    Returns:
        t1, t2 in the new basis; stops if the two bases do not span the same spaces.
    """
    nvir = mo_old.shape[1] - nocc
    rot_occ = mo_old[:, :nocc].conj().T @ mo_new[:, :nocc]
    rot_vir = mo_old[:, nocc:].conj().T @ mo_new[:, nocc:]
    err_occ = np.max(np.abs(rot_occ.conj().T @ rot_occ - np.eye(nocc)))
    err_vir = np.max(np.abs(rot_vir.conj().T @ rot_vir - np.eye(nvir)))
    err = max(err_occ, err_vir)
    print("basis rotation: max deviation from unitarity %.2e (occupied/virtual spaces must coincide)" % err)
    if err > 1e-6:
        sys.exit("--newmo does not span the same occupied/virtual spaces as --mo")
    t1_new = rot_vir.conj().T @ t1 @ rot_occ
    t2_new = np.einsum("aA,bB,abij,iI,jJ->ABIJ", rot_vir.conj(), rot_vir.conj(), t2, rot_occ, rot_occ, optimize=True)
    return t1_new, t2_new


def get_orbital_energies(moe_option, fock_mo, nocc, moe_eff):
    """Orbital energies for the (T) denominators.

    Args:
        moe_option: 'fock', 'eff' or the path of a pickled NSO array
        fock_mo: bare Fock matrix in the final basis
        nocc: number of occupied spin orbitals
        moe_eff: effective-Fock energies of the final basis (None if not available)
    """
    if moe_option == "fock":
        fock_oo = fock_mo[:nocc, :nocc]
        fock_vv = fock_mo[nocc:, nocc:]
        offdiag_oo = np.max(np.abs(fock_oo - np.diag(np.diag(fock_oo))))
        offdiag_vv = np.max(np.abs(fock_vv - np.diag(np.diag(fock_vv))))
        print("max |off-diagonal oo/vv Fock| in this basis: %.2e" % max(offdiag_oo, offdiag_vv))
        return np.real(np.diag(fock_mo)).copy()
    if moe_option == "eff":
        if moe_eff is None:
            sys.exit("--moe eff needs a 5-entry pickle [MO, fsp, fpg, fk, moe] (PHFFock.py output)")
        return moe_eff.copy()
    moe = np.real(np.asarray(pickle.load(open(moe_option, "rb")))).copy()
    assert moe.shape == (fock_mo.shape[0],)
    return moe


def get_t_solver(pcc_module, tau):
    """(T) entry of the Fortran module: pgcc_t, or pgcc_t_safe (zero-overlap-safe) when tau is given (PCC_ZO_TAU)."""
    if tau is None:
        return pcc_module.pgcc_t
    os.environ["PCC_ZO_TAU"] = repr(tau)
    print("zero-overlap-safe coherent (T) (pgcc_t_safe), tau =", tau, flush=True)
    return pcc_module.pgcc_t_safe


def main():
    args = parse_args()
    sys.path.insert(0, args.libdir)
    pcc_module = __import__(args.module)

    fc = FrozenCoreData(args.fc)
    nao = fc.nao
    nso = 2 * nao
    nocc = fc.nocc
    nvir = nso - nocc
    if args.J is None:
        args.J = fc.spin2 / 2

    # amplitudes and their orbital basis
    mo, fsp, fpg_cc, fk = load_phf_solution(args.mo)
    t1, t2 = pickle.load(open(args.t1t2, "rb"))
    t1 = np.array(t1, dtype=complex)
    t2 = np.array(t2, dtype=complex)
    assert t1.shape == (nvir, nocc) and t2.shape == (nvir, nvir, nocc, nocc), (t1.shape, t2.shape)
    if args.newmo:
        mo_new = load_phf_solution(args.newmo)[0]
        t1, t2 = rotate_amplitudes(t1, t2, mo, mo_new, nocc)
        mo = mo_new
        moe_eff = load_effective_fock_energies(args.newmo)
    else:
        moe_eff = load_effective_fock_energies(args.mo)

    # projectors (point-group coefficients set to one, as in the original driver)
    ncisp, roota, rootb, rooty, weightsp, rot_alpha, rot_beta = get_spin_projection(args.sp, args.grid, nao, args.J)
    if fsp is None:
        fsp = np.ones(ncisp, dtype=complex)
    ncipg, weightpg, rot_pg, fpg = get_pg_projection(fc, args.pg, None, nso)

    # integrals, Fock matrix and orbital energies in the final basis
    eri_ao = unpack_eri(fc.eri, nao)
    hone_so = np.zeros([nso, nso])
    hone_so[:nao, :nao] = fc.h1
    hone_so[nao:, nao:] = fc.h1
    h1_mo = ao2mo(hone_so, mo, 2)
    h2_mo = dirac_mo_eri_F(eri_ao, mo, nao)
    fock_mo = mo.conj().T @ ghf_fock_ao(fc.h1, eri_ao, mo[:, :nocc], nao) @ mo
    moe = get_orbital_energies(args.moe, fock_mo, nocc, moe_eff)
    del eri_ao
    for i in range(args.grid[0]):
        rot_alpha[i] = ao2mo(rot_alpha[i], mo, 2)
    for i in range(args.grid[1]):
        rot_beta[i] = ao2mo(rot_beta[i], mo, 2)
    for i in range(len(rot_pg)):
        rot_pg[i] = ao2mo(rot_pg[i], mo, 2)

    from mpi4py import MPI
    comm = MPI.COMM_WORLD
    mpirank = comm.Get_rank()
    mpisize = comm.Get_size()
    t_solver = get_t_solver(pcc_module, args.zero_overlap)
    et_coh, et_inc, s00, max_zbar_dev = t_solver(h1_mo, h2_mo, nao, nocc, args.J, args.sp, args.grid,
                                                          rot_alpha, rot_beta, rot_pg, roota, rootb, rooty,
                                                          weightsp, weightpg, fsp, fpg, t1, t2, moe,
                                                          int(args.inc), comm.py2f())
    et_coh = complex(et_coh)
    s00 = complex(s00)

    if args.moe in ("fock", "eff"):
        moe_label = args.moe
    else:
        moe_label = os.path.basename(args.moe)
    result = {"mpi_ranks": mpisize, "ET_coherent": et_coh.real, "ET_coherent_imag": et_coh.imag,
              "ET_incoherent": complex(et_inc).real if args.inc else None, "S00": s00.real, "S00_imag": s00.imag,
              "max_Zbar_minus_T": float(max_zbar_dev), "moe": moe_label,
              "newmo": os.path.basename(args.newmo) if args.newmo else None,
              "sp": args.sp, "grid": args.grid, "pg": args.pg, "npg": len(rot_pg), "zero_overlap": args.zero_overlap}
    if args.epcc is not None:
        result["E_PCCSD"] = args.epcc
        result["E_PCCSD_T_coherent"] = args.epcc + result["ET_coherent"]
        if args.inc:
            result["E_PCCSD_T_incoherent"] = args.epcc + result["ET_incoherent"]
    result_path = os.path.join(os.path.dirname(args.t1t2), "result_T%s.json" % args.tag)
    if mpirank == 0:
        json.dump(result, open(result_path, "w"), indent=1)
        print("RESULT_T", json.dumps(result), flush=True)


if __name__ == "__main__":
    main()
