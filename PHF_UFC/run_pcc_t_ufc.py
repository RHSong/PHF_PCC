"""(T) correction for PCC from saved PCCSD amplitudes with the UHF frozen core (FC_SO.p) or a standard FC.p.

Same options and output as SUCCSD_new/run_pcc_t.py (which provides the amplitude rotation and the choice of the
orbital energies), without --pg (no point-group projection here). The integrals come from the spin-block FC_SO
Hamiltonian and the spin rotations are transformed to the alpha / beta active orbital sets (ufc_core.py).

Usage
-----
    coherent (T):
    python PHF_UFC/run_pcc_t_ufc.py --fc FC_SO.p --mo pcc/MO_used.p --t1t2 pcc/T1T2.p --newmo SGHFFock.p \
        --moe eff --sp 2 --grid 1 8
    +(T) (orbitals from semicanon_plusT_ufc.py):
    python PHF_UFC/run_pcc_t_ufc.py --fc FC_SO.p --mo pcc/MO_used.p --t1t2 pcc/T1T2.p --newmo MO_sc.p \
        --moe fock --sp 0 --grid 1 1
"""
import argparse
import json
import os
import pickle
import sys

import numpy as np

from run_pcc_ufc import REPO_ROOT, UfcData
from run_pcc import load_phf_solution, get_spin_projection, get_pg_projection
from run_pcc_t import load_effective_fock_energies, rotate_amplitudes, get_orbital_energies, get_t_solver
from ao2mo import ao2mo
from ufc_core import hone_so, ghf_fock_so, dirac_mo_eri_so, spin_rot_list_so


def parse_args():
    """Command-line options (as run_pcc_t.py, without --pg); file paths are made absolute."""
    parser = argparse.ArgumentParser(description="(T) correction for PCC from saved amplitudes, UHF frozen core")
    parser.add_argument("--module", default="PUCC_X2", help="f2py module name")
    parser.add_argument("--libdir", default=os.path.join(REPO_ROOT, "SUCCSD_new"), help="folder of the module")
    parser.add_argument("--fc", required=True, help="FC_SO.p (or FC.p)")
    parser.add_argument("--mo", required=True, help="orbital basis of T1/T2 (MO_used.p)")
    parser.add_argument("--t1t2", required=True, help="T1T2.p")
    parser.add_argument("--newmo", default=None, help="basis to rotate T1/T2 into")
    parser.add_argument("--moe", default="fock", help="fock | eff | pickle of NSO orbital energies")
    parser.add_argument("--sp", type=int, default=2, help="spin projection: 0 none, 1 SGHF, 2 SUHF, 3 Sz")
    parser.add_argument("--grid", type=int, nargs=2, default=[1, 8], help="spin-projection grid")
    parser.add_argument("--J", type=float, default=None, help="total spin (default from the FC file)")
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


def main():
    args = parse_args()
    sys.path.insert(0, args.libdir)
    pcc_module = __import__(args.module)

    fc = UfcData(args.fc)
    if fc.is_so:
        print("integrals: spin-orbital frozen core (FC_SO)", flush=True)
    else:
        print("integrals: spatial (FC.p)", flush=True)
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

    # projectors: spin rotations transformed to the alpha / beta orbital sets, no point group
    ncisp, roota, rootb, rooty, weightsp, rot_alpha, rot_beta = get_spin_projection(args.sp, args.grid, nao, args.J)
    rot_alpha = spin_rot_list_so(rot_alpha, fc.ovlp_ab, nao)
    rot_beta = spin_rot_list_so(rot_beta, fc.ovlp_ab, nao)
    if fsp is None:
        fsp = np.ones(ncisp, dtype=complex)
    ncipg, weightpg, rot_pg, fpg = get_pg_projection(fc, 0, None, nso)

    # integrals, Fock matrix and orbital energies in the final basis
    h1_mo = ao2mo(hone_so(fc.h_a, fc.h_b), mo, 2)
    h2_mo = dirac_mo_eri_so(fc.eri_aa, fc.eri_ab, fc.eri_bb, mo, nao)
    fock_so = ghf_fock_so(fc.h_a, fc.h_b, fc.eri_aa, fc.eri_ab, fc.eri_bb, mo[:, :nocc], nao)
    fock_mo = mo.conj().T @ fock_so @ mo
    moe = get_orbital_energies(args.moe, fock_mo, nocc, moe_eff)
    for i in range(args.grid[0]):
        rot_alpha[i] = ao2mo(rot_alpha[i], mo, 2)
    for i in range(args.grid[1]):
        rot_beta[i] = ao2mo(rot_beta[i], mo, 2)
    rot_pg[0] = ao2mo(rot_pg[0], mo, 2)

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
              "newmo": os.path.basename(args.newmo) if args.newmo else None, "sp": args.sp, "grid": args.grid,
              "zero_overlap": args.zero_overlap}
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
