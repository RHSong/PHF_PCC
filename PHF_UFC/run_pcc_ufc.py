"""PCCSD on a saved PHF solution with the UHF frozen core (FC_SO.p of ufc_core.py) or a standard FC.p.

Same options and outputs as SUCCSD_new/run_pcc.py, which provides the shared steps (PHF solution, projectors,
work folder, restart guess, iteration count). Only the integrals differ: with FC_SO.p the one-electron operator is
spin dependent (h_alpha, h_beta), the two-electron integrals are kept per spin block, and the spin rotations are
transformed to the alpha / beta active orbital sets (ufc_core.py). A standard FC.p gives the same result as
run_pcc.py.

Usage
-----
    PCC_MOE_FILE=$PWD/moe.txt python PHF_UFC/run_pcc_ufc.py --module PUCC_X2 --libdir SUCCSD_new \
        --fc FC_SO.p --mo SGHFFock.p --sp 2 --grid 1 8 --semicanon 0 --work pcc

Run it from the repository (PHF_UFC/run_pcc_ufc.py) or from a case folder that holds a copy of the repository
layout; it needs SUCCSD_new/run_pcc.py and PHF_UFC/ufc_core.py.
"""
import argparse
import json
import os
import pickle
import sys
import time

import numpy as np


def find_repo_root():
    """Repository root: the folder with SUCCSD_new/run_pcc.py, searched next to and above this file and in the cwd."""
    here = os.path.dirname(os.path.abspath(__file__))
    for candidate in (os.path.dirname(here), here, os.getcwd()):
        if os.path.exists(os.path.join(candidate, "SUCCSD_new", "run_pcc.py")):
            return candidate
    sys.exit("run_pcc_ufc.py: cannot find the repository (SUCCSD_new/run_pcc.py); run it from PHF_UFC/")


REPO_ROOT = find_repo_root()
sys.path.insert(0, os.path.join(REPO_ROOT, "SUCCSD_new"))
sys.path.insert(0, os.path.join(REPO_ROOT, "PHF_UFC"))
from run_pcc import (load_phf_solution, get_spin_projection, get_pg_projection, prepare_work_dir,
                     write_guess_from_neighbour, count_iterations)
from ao2mo import ao2mo
from lean_ints import semicanon_from_fock
from ufc_core import load_fc, hone_so, ghf_fock_so, dirac_mo_eri_so, spin_rot_list_so


def parse_args():
    """Command-line options (as run_pcc.py, without --lean); file and folder paths are made absolute."""
    parser = argparse.ArgumentParser(description="PCCSD on a saved PHF solution, UHF frozen core")
    parser.add_argument("--module", required=True, help="f2py module name, e.g. PUCC_X2")
    parser.add_argument("--libdir", required=True, help="folder that holds the module")
    parser.add_argument("--fc", required=True, help="FC_SO.p (or FC.p)")
    parser.add_argument("--mo", required=True, help="PHF solution (SGHF.p, SGHFFock.p)")
    parser.add_argument("--sp", type=int, default=2, help="spin projection: 0 none, 1 SGHF, 2 SUHF, 3 Sz")
    parser.add_argument("--grid", type=int, nargs=2, default=[1, 8], help="spin-projection grid")
    parser.add_argument("--J", type=float, default=None, help="total spin (default from the FC file)")
    parser.add_argument("--cmplx", type=int, default=0, help="complex-conjugation projection (CmplxConj)")
    parser.add_argument("--pg", type=int, default=0, help="1: point-group projection from the FC file")
    parser.add_argument("--work", required=True, help="work folder")
    parser.add_argument("--nbroy", type=int, default=30, help="Broyden history length")
    parser.add_argument("--semicanon", type=int, default=1, help="1: semicanonicalize the bare Fock first")
    parser.add_argument("--keep-chk", action="store_true", help="restart from Chk_PUCC_* in --work")
    parser.add_argument("--guess-from", default=None, help="work folder of a neighbouring geometry")
    args = parser.parse_args()
    args.fc = os.path.abspath(args.fc)
    args.mo = os.path.abspath(args.mo)
    args.libdir = os.path.abspath(args.libdir)
    args.work = os.path.abspath(args.work)
    if args.guess_from:
        args.guess_from = os.path.abspath(args.guess_from)
    return args


class UfcData:
    """Frozen-core Hamiltonian from FC_SO.p (or FC.p) via ufc_core.load_fc.

    Attributes:
        enuc: nuclear repulsion plus frozen-core energy
        h_a, h_b: active one-electron integrals of the alpha and beta orbital sets
        eri_aa, eri_ab, eri_bb: Mulliken spin-block integrals (pq|rs)
        ovlp_ab: spatial overlap <p^alpha|q^beta> of the two active sets (identity for FC.p)
        ncipg, weightpg, rot_pg: point-group projector
        nao, nocc, spin2: active orbitals per spin, occupied spin orbitals, 2S
        is_so: True for FC_SO.p
    """

    def __init__(self, path):
        fc = load_fc(path)
        self.enuc = fc["Enuc"]
        self.h_a = fc["h_a"]
        self.h_b = fc["h_b"]
        self.eri_aa = fc["eri_aa"]
        self.eri_ab = fc["eri_ab"]
        self.eri_bb = fc["eri_bb"]
        self.ovlp_ab = fc["O_ab"]
        self.ncipg = fc["ncipg"]
        self.weightpg = fc["weightpg"]
        self.rot_pg = fc["Rpg"]
        self.nao = int(fc["NAO"])
        self.nocc = int(fc["NOccSO"])
        self.spin2 = int(fc["Spin2"])
        self.is_so = fc["so"]


def make_mo_integrals_so(fc, mo, nocc_alpha, sp, semicanon):
    """One- and two-electron integrals in the CC orbital basis from the spin-block integrals.

    Returns:
        mo: the CC orbitals (semicanonical in the bare Fock if semicanon), h1_mo, h2_mo (<pq||rs>, Fortran order)
    """
    nao = fc.nao
    nocc = fc.nocc
    if semicanon:
        fock_so = ghf_fock_so(fc.h_a, fc.h_b, fc.eri_aa, fc.eri_ab, fc.eri_bb, mo[:, :nocc], nao)
        mo = semicanon_from_fock(fock_so, mo, nocc, sp, nocc_alpha)
    h1_mo = ao2mo(hone_so(fc.h_a, fc.h_b), mo, 2)
    h2_mo = dirac_mo_eri_so(fc.eri_aa, fc.eri_ab, fc.eri_bb, mo, nao)
    return mo, h1_mo, h2_mo


def main():
    args = parse_args()
    sys.path.insert(0, args.libdir)
    pcc_module = __import__(args.module)

    # frozen-core integrals and the PHF reference
    fc = UfcData(args.fc)
    if fc.is_so:
        print("integrals: spin-orbital frozen core (FC_SO)", flush=True)
    else:
        print("integrals: spatial (FC.p)", flush=True)
    nao = fc.nao
    nso = 2 * nao
    nocc = fc.nocc
    nocc_alpha = (nocc + fc.spin2) // 2
    if args.J is None:
        args.J = fc.spin2 / 2
    if args.sp == 2 and abs(args.J - fc.spin2 / 2) > 1e-8:
        sys.exit("SUHF-based PCC needs J = (n_alpha - n_beta)/2")
    mo, fsp, fpg_phf, fk = load_phf_solution(args.mo)

    # projectors: spin rotations transformed to the alpha / beta orbital sets, point group, Rk = 1
    ncisp, roota, rootb, rooty, weightsp, rot_alpha, rot_beta = get_spin_projection(args.sp, args.grid, nao, args.J)
    rot_alpha = spin_rot_list_so(rot_alpha, fc.ovlp_ab, nao)
    rot_beta = spin_rot_list_so(rot_beta, fc.ovlp_ab, nao)
    if fsp is None:
        fsp = np.ones(ncisp, dtype=complex)
    ncipg, weightpg, rot_pg, fpg = get_pg_projection(fc, args.pg, fpg_phf, nso)
    npg = len(weightpg)
    ncik = 1
    fk = np.ones(1, dtype=complex)
    rot_k = np.zeros([1, nso, nso], dtype=complex)
    rot_k[0] = np.eye(nso)

    # integrals and projection operators in the CC orbital basis
    mo, h1_mo, h2_mo = make_mo_integrals_so(fc, mo, nocc_alpha, args.sp, args.semicanon)
    for i in range(args.grid[0]):
        rot_alpha[i] = ao2mo(rot_alpha[i], mo, 2)
    for i in range(args.grid[1]):
        rot_beta[i] = ao2mo(rot_beta[i], mo, 2)
    for i in range(npg):
        rot_pg[i] = ao2mo(rot_pg[i], mo, 2)
    x_orth = np.eye(nso, dtype=complex)
    x_orth_inv = np.eye(nso, dtype=complex)

    from mpi4py import MPI
    comm = MPI.COMM_WORLD
    mpirank = comm.Get_rank()
    mpisize = comm.Get_size()
    prepare_work_dir(args.work, mpirank, args.keep_chk or args.guess_from)
    if args.guess_from:
        write_guess_from_neighbour(args.guess_from, mo, nocc)

    print("module %s  NSO %d NOcc %d  SP %d grid %s  ncisp %d  npg %d"
          % (args.module, nso, nocc, args.sp, args.grid, ncisp, npg), flush=True)
    time_start = time.time()
    pcc_out = pcc_module.pgcc(mo, h1_mo, h2_mo, nao, nocc, args.J, args.cmplx, args.sp, args.grid,
                              rot_alpha, rot_beta, rot_pg, rot_k, roota, rootb, rooty, weightsp, weightpg,
                              fsp, fpg, fk, args.nbroy, x_orth, x_orth_inv, fc.enuc, comm.py2f(),
                              nso, npg, ncisp, ncipg, ncik)
    wall_time = time.time() - time_start
    e_pcc = complex(pcc_out[0])
    t1 = pcc_out[2]
    t2 = pcc_out[3]

    result = {"mpi_ranks": mpisize, "module": args.module, "E_PCC": e_pcc.real, "E_imag": e_pcc.imag,
              "wall_s": wall_time, "iterations": count_iterations("Output"),
              "maxT1": float(np.max(np.abs(t1))), "maxT2": float(np.max(np.abs(t2))),
              "NSO": nso, "NOcc": nocc, "grid": args.grid}
    if mpirank == 0:
        json.dump(result, open("result.json", "w"), indent=1)
        pickle.dump([t1, t2], open("T1T2.p", "wb"), protocol=4)
        pickle.dump([mo, fsp, fpg, fk], open("MO_used.p", "wb"))
    print("RESULT", json.dumps(result), flush=True)


if __name__ == "__main__":
    main()
