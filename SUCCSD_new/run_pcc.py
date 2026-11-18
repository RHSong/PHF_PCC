"""PCCSD on a saved PHF solution (Fortran module PUCC_X2 of SUCCSD_new).

Usage
-----
    python run_pcc.py --module PUCC_X2 --libdir SUCCSD_new --fc FC.p --mo SGHFFock.p \
                      --sp 2 --grid 1 8 --semicanon 0 --work pcc

Inputs
------
--fc          FC.p of Main_FC.py: [Enuc, h1, eri, ncipg, weightpg, Rpg, NAO, NOccSO(, 2S)]
--mo          PHF solution [MO, fsp, fpg, fk(, moe)] (SGHF.p, SGHFFock.p) or a bare MO matrix
--sp, --grid  spin projection (0 none, 1 SGHF, 2 SUHF, 3 Sz) and its grid
--J           total spin of the projected state (default: S = (n_alpha - n_beta)/2 from FC.p)
--pg 1        also apply the point-group projector stored in FC.p
--semicanon   1: semicanonicalize the bare Fock in the PHF orbital spaces first;
              0: keep the orbitals of --mo (use 0 for the PHF effective-Fock orbitals)
--lean        1: lean_ints integral path (one complex NSO^4 array); 0: old dense path
--guess-from  work folder of a converged run at a neighbouring geometry; its amplitudes, rotated into the
              present orbitals, start the solver
--keep-chk    keep the Fortran restart files Chk_PUCC_* of an earlier run in --work (restart from them)
--zero-overlap TAU
              zero-overlap-safe kernels (pgcc_safe): on grid points where the overlap matrix has singular values
              below TAU the Thouless step uses the deformed determinant (FPUCC::EvalKernelSafe)

Environment (read by the Fortran code)
--------------------------------------
PCC_MOE_FILE  text file with NSO orbital energies in the CC orbital order (e.g. the PHF effective-Fock
              eigenvalues) used in the Broyden preconditioner instead of the bare-Fock diagonal
PCC_CYCMAX    maximum number of Broyden cycles
PCC_MIX       Broyden mixing parameter (default 0.5)
PCC_ZO_TAU    zero-overlap threshold (set by --zero-overlap)

Outputs (in --work)
-------------------
result.json   energy, wall time, number of iterations, largest amplitudes
T1T2.p        [T1, T2], T1[a, i], T2[a, b, i, j]
MO_used.p     [MO, fsp, fpg, fk]: the exact orbital basis of T1/T2 (input of run_pcc_t.py)
Output, Kern.log, RKstat, Chk_PUCC_*: Fortran output and restart files
With several MPI ranks, rank r > 0 works in --work/rank<r>.
"""
import argparse
import json
import os
import pickle
import sys
import time

import numpy as np

REPO_ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))   # this file: <repo>/SUCCSD_new/
sys.path.insert(0, REPO_ROOT)
from ao2mo import ao2mo
from makeH import Mulliken2Dirac, SemiCanon
from Spin_Proj import Spin_Proj
from lean_ints import unpack_eri, ghf_fock_ao, semicanon_from_fock, dirac_mo_eri_F

RESTART_FILES = ["Chk_PUCC_Ene", "Chk_PUCC_T1", "Chk_PUCC_T2"]
OUTPUT_FILES = ["Output", "RKstat", "Kern.log"]


def parse_args():
    """Command-line options (see the module docstring); file and folder paths are made absolute."""
    parser = argparse.ArgumentParser(description="PCCSD on a saved PHF solution")
    parser.add_argument("--module", required=True, help="f2py module name, e.g. PUCC_X2")
    parser.add_argument("--libdir", required=True, help="folder that holds the module")
    parser.add_argument("--fc", required=True, help="FC.p")
    parser.add_argument("--mo", required=True, help="PHF solution (SGHF.p, SGHFFock.p)")
    parser.add_argument("--sp", type=int, default=2, help="spin projection: 0 none, 1 SGHF, 2 SUHF, 3 Sz")
    parser.add_argument("--grid", type=int, nargs=2, default=[1, 8], help="spin-projection grid")
    parser.add_argument("--J", type=float, default=None, help="total spin (default from FC.p)")
    parser.add_argument("--cmplx", type=int, default=0, help="complex-conjugation projection (CmplxConj)")
    parser.add_argument("--pg", type=int, default=0, help="1: point-group projection from FC.p")
    parser.add_argument("--work", required=True, help="work folder")
    parser.add_argument("--nbroy", type=int, default=30, help="Broyden history length")
    parser.add_argument("--semicanon", type=int, default=1, help="1: semicanonicalize the bare Fock first")
    parser.add_argument("--lean", type=int, default=1, help="1: lean_ints integral path, 0: dense path")
    parser.add_argument("--keep-chk", action="store_true", help="restart from Chk_PUCC_* in --work")
    parser.add_argument("--guess-from", default=None, help="work folder of a neighbouring geometry")
    parser.add_argument("--zero-overlap", type=float, default=None, metavar="TAU",
                        help="zero-overlap-safe kernels (pgcc_safe) with singular-value threshold TAU")
    args = parser.parse_args()
    args.fc = os.path.abspath(args.fc)
    args.mo = os.path.abspath(args.mo)
    args.libdir = os.path.abspath(args.libdir)
    args.work = os.path.abspath(args.work)
    if args.guess_from:
        args.guess_from = os.path.abspath(args.guess_from)
    return args


class FrozenCoreData:
    """Contents of FC.p (Main_FC.py).

    Attributes:
        enuc: nuclear repulsion plus frozen-core energy
        h1, eri: one- and two-electron integrals over the NAO active spatial orbitals (eri packed or 4-index)
        ncipg, weightpg, rot_pg: point-group projector (number of irrep components, weights, operators)
        nao, nocc: number of active spatial orbitals and of occupied spin orbitals
        spin2: 2S = n_alpha - n_beta
    """

    def __init__(self, path):
        fc = pickle.load(open(path, "rb"))
        self.enuc, self.h1, self.eri, self.ncipg, self.weightpg, self.rot_pg = fc[:6]
        self.nao = int(fc[6])
        self.nocc = int(fc[7])
        if len(fc) > 8:
            self.spin2 = int(fc[8])
        else:
            self.spin2 = 0


def load_phf_solution(path):
    """Read the PHF orbitals and projection coefficients.

    Args:
        path: [MO, fsp, fpg, fk(, moe)] pickle, or a pickled MO matrix (or [MO])

    Returns:
        mo (complex NSO x NSO), fsp, fpg, fk (None when the file holds only orbitals)
    """
    data = pickle.load(open(path, "rb"))
    fsp = fpg = fk = None
    if isinstance(data, (list, tuple)) and len(data) >= 4:
        mo, fsp, fpg, fk = data[:4]
    elif isinstance(data, (list, tuple)):
        mo = data[0]
    else:
        mo = data
    return np.array(mo, dtype=complex), fsp, fpg, fk


def get_spin_projection(sp, grid, nao, spin_j):
    """Spin-projection grid and rotation matrices (Spin_Proj.py) in the spin-orbital AO-like basis.

    Spin_Proj reads lebedev_*.txt relative to the repository root, so it runs from there.

    Returns:
        ncisp, roota, rootb, rooty, weightsp, rot_alpha (R1), rot_beta (R2)
    """
    cwd = os.getcwd()
    os.chdir(REPO_ROOT)
    projection = Spin_Proj(sp, grid, nao, 2 * nao, spin_j)
    os.chdir(cwd)
    return projection


def get_pg_projection(fc, use_pg, fpg_phf, nso):
    """Point-group projector: the one in FC.p (use_pg = 1) or the identity.

    Returns:
        ncipg, weightpg, rot_pg (npg x NSO x NSO, complex), fpg (PHF coefficients, or ones)
    """
    if use_pg:
        rot_pg = np.array(fc.rot_pg, dtype=complex)
        weightpg = np.array(fc.weightpg)
        ncipg = fc.ncipg
        if fpg_phf is None:
            fpg = np.ones(ncipg, dtype=complex)
        else:
            fpg = np.array(fpg_phf, dtype=complex).reshape(ncipg)
        return ncipg, weightpg, rot_pg, fpg
    rot_pg = np.zeros([1, nso, nso], dtype=complex)
    rot_pg[0] = np.eye(nso)
    weightpg = np.ones([1, 1, 1])
    fpg = np.ones(1, dtype=complex)
    return 1, weightpg, rot_pg, fpg


def make_mo_integrals(fc, mo, nocc_alpha, sp, semicanon, lean):
    """One- and two-electron integrals in the CC orbital basis.

    Args:
        fc: FrozenCoreData
        mo: PHF orbitals (NSO x NSO, alpha block first in the rows)
        nocc_alpha: number of occupied alpha orbitals (SUHF semicanonicalization)
        sp: spin projection type
        semicanon: 1 to semicanonicalize the bare Fock in the occupied / virtual spaces of mo first
        lean: 1 for the lean_ints path, 0 for the dense NSO^4 real + complex path

    Returns:
        mo: the CC orbitals, h1_mo (NSO x NSO), h2_mo (<pq||rs>, Fortran order)
    """
    nao = fc.nao
    nso = 2 * nao
    nocc = fc.nocc
    hone_so = np.zeros([nso, nso])
    hone_so[:nao, :nao] = fc.h1
    hone_so[nao:, nao:] = fc.h1
    if lean:
        eri_ao = unpack_eri(fc.eri, nao)
        if semicanon:
            fock_ao = ghf_fock_ao(fc.h1, eri_ao, mo[:, :nocc], nao)
            mo = semicanon_from_fock(fock_ao, mo, nocc, sp, nocc_alpha)
        h1_mo = ao2mo(hone_so, mo, 2)
        h2_mo = dirac_mo_eri_F(eri_ao, mo, nao)
        return mo, h1_mo, h2_mo
    # dense path: Mulliken (pq|rs) in the 'spatial x spin' basis, blocks alpha-alpha|alpha-alpha, beta-beta|beta-beta,
    # beta-beta|alpha-alpha, alpha-alpha|beta-beta
    htwo_so = np.zeros([nso] * 4)
    for spin_pq, spin_rs in [(0, 0), (1, 1), (1, 0), (0, 1)]:
        pq = slice(spin_pq * nao, (spin_pq + 1) * nao)
        rs = slice(spin_rs * nao, (spin_rs + 1) * nao)
        htwo_so[pq, pq, rs, rs] = fc.eri
    htwo_so = Mulliken2Dirac(htwo_so)
    if semicanon:
        mo = SemiCanon(hone_so, htwo_so, mo, nocc, sp, nocc_alpha)
    h1_mo = ao2mo(hone_so, mo, 2)
    h2_mo = np.asfortranarray(ao2mo(htwo_so, mo, 4))
    return mo, h1_mo, h2_mo


def prepare_work_dir(work, mpirank, keep_restart):
    """Create and enter the work folder of this rank (rank 0: work, rank r: work/rank<r>).

    Old Fortran output files are removed; the restart files Chk_PUCC_* too, unless keep_restart.
    """
    if mpirank == 0:
        work_dir = work
    else:
        work_dir = os.path.join(work, "rank%d" % mpirank)
    os.makedirs(work_dir, exist_ok=True)
    os.chdir(work_dir)
    stale_files = list(OUTPUT_FILES)
    if not keep_restart:
        stale_files += RESTART_FILES
    for name in stale_files:
        if os.path.exists(name):
            os.remove(name)


def write_guess_from_neighbour(guess_dir, mo_new, nocc):
    """Write the Fortran restart files from the amplitudes of a neighbouring geometry.

    The old amplitudes live in the old CC orbitals C_old (guess_dir/MO_used.p), the new run uses C_new; both are
    expanded in the frozen-core orbital basis, which orbital following keeps aligned between geometries. With
    ovlp = C_new^H C_old (old orbital q = sum_p new p ovlp[p, q]), and the occupied-virtual blocks dropped,
        T1'[b, j] = sum_ai ovlp_vir[b, a] T1[a, i] conj(ovlp_occ[j, i])
        T2'[c, d, k, l] = sum_abij ovlp_vir[c, a] ovlp_vir[d, b] T2[a, b, i, j] conj(ovlp_occ[k, i]) conj(ovlp_occ[l, j])

    Args:
        guess_dir: work folder of the converged neighbouring run (MO_used.p, T1T2.p, result.json)
        mo_new: CC orbitals of this run
        nocc: number of occupied spin orbitals

    Writes Chk_PUCC_Ene, Chk_PUCC_T1, Chk_PUCC_T2 into the current folder.
    """
    from scipy.io import FortranFile
    mo_old = np.array(pickle.load(open(os.path.join(guess_dir, "MO_used.p"), "rb"))[0], dtype=complex)
    t1_old, t2_old = pickle.load(open(os.path.join(guess_dir, "T1T2.p"), "rb"))
    t1_old = np.array(t1_old, dtype=complex)
    t2_old = np.array(t2_old, dtype=complex)
    e_old = json.load(open(os.path.join(guess_dir, "result.json")))["E_PCC"]

    ovlp = mo_new.conj().T @ mo_old
    ovlp_occ = ovlp[:nocc, :nocc]
    ovlp_vir = ovlp[nocc:, nocc:]
    t1_new = ovlp_vir @ t1_old @ ovlp_occ.conj().T
    t2_new = np.einsum("ca,db,abij,ki,lj->cdkl", ovlp_vir, ovlp_vir, t2_old, ovlp_occ.conj(), ovlp_occ.conj(),
                       optimize=True)

    restart_data = [("Chk_PUCC_Ene", np.array([e_old + 0j])), ("Chk_PUCC_T1", t1_new), ("Chk_PUCC_T2", t2_new)]
    for name, array in restart_data:
        restart_file = FortranFile(name, "w")
        restart_file.write_record(np.asfortranarray(array).ravel(order="F"))
        restart_file.close()

    sv_occ = np.linalg.svd(ovlp_occ, compute_uv=False)
    sv_vir = np.linalg.svd(ovlp_vir, compute_uv=False)
    print("guess from %s: E_old %.8f" % (guess_dir, e_old))
    print("  occupied overlap singular values %.5f .. %.5f" % (sv_occ.min(), sv_occ.max()))
    print("  virtual overlap singular values  %.5f .. %.5f" % (sv_vir.min(), sv_vir.max()))
    print("  max|T1| %.4f -> %.4f   |T2| %.5f -> %.5f" % (np.abs(t1_old).max(), np.abs(t1_new).max(),
                                                         np.linalg.norm(t2_old), np.linalg.norm(t2_new)), flush=True)


def get_cc_solver(pcc_module, tau):
    """PCCSD entry of the Fortran module: pgcc, or pgcc_safe (zero-overlap-safe kernels) when tau is given.

    The Fortran code reads tau from the environment variable PCC_ZO_TAU.
    """
    if tau is None:
        return pcc_module.pgcc
    os.environ["PCC_ZO_TAU"] = repr(tau)
    print("zero-overlap-safe kernels (pgcc_safe), tau =", tau, flush=True)
    return pcc_module.pgcc_safe


def count_iterations(output_path):
    """Number of Broyden iterations in the Fortran Output file (lines '* <iter> <energy> ...'); -1 if absent."""
    if not os.path.exists(output_path):
        return -1
    niter = 0
    for line in open(output_path):
        words = line.split()
        if not line.strip().startswith("*") or "." not in line or len(words) < 2:
            continue
        if words[1].replace("-", "").replace(".", "").isdigit():
            niter += 1
    return niter


def main():
    args = parse_args()
    sys.path.insert(0, args.libdir)
    pcc_module = __import__(args.module)

    # frozen-core integrals and the PHF reference
    fc = FrozenCoreData(args.fc)
    nao = fc.nao
    nso = 2 * nao
    nocc = fc.nocc
    nocc_alpha = (nocc + fc.spin2) // 2
    if args.J is None:
        args.J = fc.spin2 / 2
    if args.sp == 2 and abs(args.J - fc.spin2 / 2) > 1e-8:
        sys.exit("SUHF-based PCC needs J = (n_alpha - n_beta)/2")
    mo, fsp, fpg_phf, fk = load_phf_solution(args.mo)

    # projectors: spin (R1, R2 over the grid), point group, no complex conjugation (Rk = 1)
    ncisp, roota, rootb, rooty, weightsp, rot_alpha, rot_beta = get_spin_projection(args.sp, args.grid, nao, args.J)
    if fsp is None:
        fsp = np.ones(ncisp, dtype=complex)
    ncipg, weightpg, rot_pg, fpg = get_pg_projection(fc, args.pg, fpg_phf, nso)
    npg = len(weightpg)
    ncik = 1
    fk = np.ones(1, dtype=complex)
    rot_k = np.zeros([1, nso, nso], dtype=complex)
    rot_k[0] = np.eye(nso)

    # integrals and projection operators in the CC orbital basis
    mo, h1_mo, h2_mo = make_mo_integrals(fc, mo, nocc_alpha, args.sp, args.semicanon, args.lean)
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
    cc_solver = get_cc_solver(pcc_module, args.zero_overlap)
    pcc_out = cc_solver(mo, h1_mo, h2_mo, nao, nocc, args.J, args.cmplx, args.sp, args.grid,
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
              "NSO": nso, "NOcc": nocc, "grid": args.grid, "zero_overlap": args.zero_overlap}
    if mpirank == 0:
        json.dump(result, open("result.json", "w"), indent=1)
        pickle.dump([t1, t2], open("T1T2.p", "wb"), protocol=4)
        pickle.dump([mo, fsp, fpg, fk], open("MO_used.p", "wb"))
    print("RESULT", json.dumps(result), flush=True)


if __name__ == "__main__":
    main()
