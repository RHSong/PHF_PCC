"""PCCSD driver: run the PCC module (PUCC_X2 of SUCCSD_new) on a saved PHF solution.
python run_pcc.py --module PUCC_X2 --libdir SUCCSD_new --fc FC.p --mo SGHF.p --sp 2 --grid 1 8 --work run_dir
Writes run_dir/result.json with energy, wall time, iterations, max residual; Fortran output files land in run_dir."""
import sys, os, time, json, argparse, pickle
import numpy as np
ap = argparse.ArgumentParser()
ap.add_argument("--module", required=True); ap.add_argument("--libdir", required=True)
ap.add_argument("--fc", required=True); ap.add_argument("--mo", required=True)
ap.add_argument("--sp", type=int, default=2); ap.add_argument("--grid", type=int, nargs=2, default=[1, 8])
ap.add_argument("--J", type=float, default=None); ap.add_argument("--cmplx", type=int, default=0)   # J: total spin, default S from FC.p
ap.add_argument("--work", required=True); ap.add_argument("--nbroy", type=int, default=30)
ap.add_argument("--semicanon", type=int, default=1)
ap.add_argument("--pg", type=int, default=0)
ap.add_argument("--nullvec", default=None)   # nullvec.txt (nullvec.py): run the null-space-projected solver (pgcc_null) with these directions
ap.add_argument("--lean", type=int, default=1)   # 1: lean_ints integral path (needed for Cr2), 0: old dense path
ap.add_argument("--keep-chk", action="store_true")   # keep Chk_PUCC_* in the work dir: the Fortran code restarts from them
ap.add_argument("--guess-from", default=None)   # work dir of a converged run at a neighbouring geometry (MO_used.p, T1T2.p,
                                                # result.json): its amplitudes, mapped to this run's orbitals, start the solver
args = ap.parse_args()
if args.guess_from: args.guess_from = os.path.abspath(args.guess_from)
if args.nullvec: args.nullvec = os.path.abspath(args.nullvec)
for k in ("fc", "mo", "libdir", "work"): setattr(args, k, os.path.abspath(getattr(args, k)))
root = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))   # PHF_PCC repo root (this file: SUCCSD_new/)
sys.path.insert(0, root); sys.path.insert(0, os.path.abspath(args.libdir))
from ao2mo import ao2mo
from makeH import Mulliken2Dirac, SemiCanon
from Spin_Proj import Spin_Proj
mod = __import__(args.module)
_fc = pickle.load(open(args.fc, "rb"))
Enuc, h1, eri, ncipg, weightpg, Rpg, NAO, NOccSO = _fc[:8]
Spin2 = int(_fc[8]) if len(_fc) > 8 else 0   # 2S = n_alpha - n_beta
NOccA = (NOccSO + Spin2) // 2
if args.J is None: args.J = Spin2 / 2
if args.sp == 2 and abs(args.J - Spin2 / 2) > 1e-8: sys.exit("SUHF-based PCC needs J = (n_alpha - n_beta)/2")
d = pickle.load(open(args.mo, "rb"))
if isinstance(d, (list, tuple)) and len(d) >= 4: MO, fsp, fpg, fk = d[:4]
else: MO = d[0] if isinstance(d, (list, tuple)) else d; fsp = fpg = fk = None
NSO = int(2 * NAO)
cwd0 = os.getcwd(); os.chdir(root)   # SpinSymm reads lebedev_*.txt relative to the repo root
ncisp, roota, rootb, rooty, weightsp, R1, R2 = Spin_Proj(args.sp, args.grid, NAO, NSO, args.J)
os.chdir(cwd0)
if args.pg:
    Rpg = np.array(Rpg, dtype=complex); weightpg = np.array(weightpg); npg = len(weightpg)
    fpg = np.ones(ncipg, dtype=complex) if fpg is None else np.array(fpg, dtype=complex).reshape(ncipg)
else:
    ncipg = 1; weightpg = np.ones([1, 1, 1]); Rpg = np.zeros([1, NSO, NSO], dtype=complex); Rpg[0] = np.eye(NSO)
    npg = 1; fpg = np.ones(1, dtype=complex)
ncik = 1
if fsp is None: fsp = np.ones(ncisp, dtype=complex)
fk = np.ones(1, dtype=complex)
HOne = np.zeros([NSO, NSO]); HOne[:NAO, :NAO] = h1; HOne[NAO:, NAO:] = h1
MO = np.array(MO, dtype=complex)
if args.lean:
    # one complex NSO^4 array only (the final <pq||rs>, Fortran order): see lean_ints.py
    from lean_ints import unpack_eri, ghf_fock_ao, semicanon_from_fock, dirac_mo_eri_F
    eri4 = unpack_eri(eri, NAO); del eri
    if args.semicanon: MO = semicanon_from_fock(ghf_fock_ao(h1, eri4, MO[:, :NOccSO], NAO), MO, NOccSO, args.sp, NOccA)
    H1 = ao2mo(HOne, MO, 2); H2 = dirac_mo_eri_F(eri4, MO, NAO); del eri4
else:
    HTwo = np.zeros([NSO] * 4)
    for a, b in [(0, 0), (1, 1), (1, 0), (0, 1)]:
        HTwo[a*NAO:(a+1)*NAO, a*NAO:(a+1)*NAO, b*NAO:(b+1)*NAO, b*NAO:(b+1)*NAO] = eri
    HTwo = Mulliken2Dirac(HTwo)
    if args.semicanon: MO = SemiCanon(HOne, HTwo, MO, NOccSO, args.sp, NOccA)
    H1 = ao2mo(HOne, MO, 2); H2 = np.asfortranarray(ao2mo(HTwo, MO, 4)); del HTwo
for i in range(args.grid[0]): R1[i] = ao2mo(R1[i], MO, 2)
for i in range(args.grid[1]): R2[i] = ao2mo(R2[i], MO, 2)
for i in range(npg): Rpg[i] = ao2mo(Rpg[i], MO, 2)
Rk = np.zeros([1, NSO, NSO], dtype=complex); Rk[0] = np.eye(NSO)
X = np.eye(NSO, dtype=complex); Xinv = np.eye(NSO, dtype=complex)
from mpi4py import MPI
fcomm = MPI.COMM_WORLD.py2f()
comm = MPI.COMM_WORLD; mpirank = comm.Get_rank(); mpisize = comm.Get_size()
# each rank runs the Fortran driver in its own folder (it writes Output/Chk files on every rank);
# rank 0 keeps the top-level work folder, which holds result.json, T1T2.p and MO_used.p
wd = args.work if mpirank == 0 else os.path.join(args.work, "rank%d" % mpirank)
os.makedirs(wd, exist_ok=True); os.chdir(wd)
for f in (["Output", "RKstat", "Kern.log"] if (args.keep_chk or args.guess_from) else ["Chk_PUCC_Ene", "Chk_PUCC_T1", "Chk_PUCC_T2", "Output", "RKstat", "Kern.log"]):
    if os.path.exists(f): os.remove(f)
if args.guess_from:
    # Amplitudes of a neighbouring geometry in its own (semicanonical) orbitals C_old; this run uses C_new, both expanded
    # in the frozen-core orbital basis, which orbital following keeps aligned between geometries.  With
    # S = C_new^H C_old (old orbital q = sum_p new p S_pq) the excitation operator becomes
    #   t'(b,j) = sum_ai S_vv(b,a) t(a,i) conj(S_oo(j,i)),  t'(cd,kl) = sum S_vv S_vv t(ab,ij) conj(S_oo) conj(S_oo);
    # the occupied-virtual blocks of S are dropped.  Written as the Fortran restart files Chk_PUCC_*.
    from scipy.io import FortranFile
    Cold = np.array(pickle.load(open(os.path.join(args.guess_from, "MO_used.p"), "rb"))[0], dtype=complex)
    T1o, T2o = [np.array(x, dtype=complex) for x in pickle.load(open(os.path.join(args.guess_from, "T1T2.p"), "rb"))]
    Eold = json.load(open(os.path.join(args.guess_from, "result.json")))["E_PCC"]
    no = NOccSO; S = MO.conj().T @ Cold; Soo = S[:no, :no]; Svv = S[no:, no:]
    T1g = Svv @ T1o @ Soo.conj().T
    T2g = np.tensordot(Svv, T2o, axes=(1, 0)); T2g = np.tensordot(Svv, T2g, axes=(1, 1)).transpose(1, 0, 2, 3)
    T2g = np.tensordot(T2g, Soo.conj(), axes=(2, 1)); T2g = np.tensordot(T2g, Soo.conj(), axes=(2, 1))
    for name, arr in (("Chk_PUCC_Ene", np.array([Eold + 0j])), ("Chk_PUCC_T1", T1g), ("Chk_PUCC_T2", T2g)):
        f = FortranFile(name, "w"); f.write_record(np.asfortranarray(arr).ravel(order="F")); f.close()
    sv = lambda A: np.linalg.svd(A, compute_uv=False)
    print("guess from %s: E_old %.8f  occ overlap sv %.5f..%.5f  vir overlap sv %.5f..%.5f  max|T1| %.4f -> %.4f  |T2| %.5f -> %.5f"
          % (args.guess_from, Eold, sv(Soo).min(), sv(Soo).max(), sv(Svv).min(), sv(Svv).max(), np.abs(T1o).max(),
             np.abs(T1g).max(), np.linalg.norm(T2o), np.linalg.norm(T2g)), flush=True)
    del T2o, T2g
print("module %s  NSO %d NOcc %d  SP %d grid %s  ncisp %d  npg %d" % (args.module, NSO, NOccSO, args.sp, args.grid, ncisp, npg), flush=True)
t0 = time.time()
if os.path.exists("nullvec.txt"): os.remove("nullvec.txt")
if args.nullvec:
    import shutil; shutil.copy(args.nullvec, "nullvec.txt"); driver = mod.pgcc_null
    print("null-space-projected solver (pgcc_null) with", args.nullvec, flush=True)
else:
    driver = mod.pgcc
out = driver(MO, H1, H2, NAO, NOccSO, args.J, args.cmplx, args.sp, args.grid, R1, R2, Rpg, Rk, roota, rootb, rooty,
             weightsp, weightpg, fsp, fpg, fk, args.nbroy, X, Xinv, Enuc, fcomm, NSO, npg, ncisp, ncipg, ncik)
dt = time.time() - t0
E = complex(out[0]); T1 = out[2]; T2 = out[3]
niter = sum(1 for l in open("Output") if l.strip().startswith("*") and "." in l and l.split()[1].replace("-", "").replace(".", "").isdigit()) if os.path.exists("Output") else -1
res = {"mpi_ranks": mpisize, "module": args.module, "E_PCC": E.real, "E_imag": E.imag, "wall_s": dt, "iterations": niter,
       "maxT1": float(np.max(np.abs(T1))), "maxT2": float(np.max(np.abs(T2))), "NSO": int(NSO), "NOcc": int(NOccSO), "grid": args.grid}
if mpirank == 0: json.dump(res, open("result.json", "w"), indent=1)
if mpirank == 0: pickle.dump([T1, T2], open("T1T2.p", "wb"), protocol=4)
if mpirank == 0: pickle.dump([MO, fsp, fpg, fk], open("MO_used.p", "wb"))   # exact basis of T1/T2 (after SemiCanon)
print("RESULT", json.dumps(res), flush=True)
