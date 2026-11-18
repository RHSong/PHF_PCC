"""(T) driver for PCC (coherent (T) and, with --sp 0 --grid 1 1 --newmo, +(T)); --pg 1: point-group projector from FC.p.
Coherent (T) for PCC from saved, converged PCCSD amplitudes (no PCCSD rerun).

python run_pcc_t.py --fc FC.p --mo MO_used.p --t1t2 T1T2.p --sp 2 --grid 1 8 \
                [--moe fock | --moe moe.p] [--newmo MO_eff.p] [--inc] [--epcc -109.0]

--mo      the exact orbital basis in which T1/T2 were converged (test_pcc.py writes MO_used.p)
--newmo   optional second basis spanning the same occupied / virtual spaces (e.g. the
          semicanonical orbitals of the PHF effective Fock); T is rotated into it
--moe     'fock' : denominators from the diagonal of the bare Fock matrix in the final basis
          file.p : pickled array of NSO orbital energies in the column order of the final basis
--inc     also evaluate the old per-grid (incoherent) (T) sum
Writes result_T.json into the directory of --t1t2.
"""
import sys, os, json, argparse, pickle
import numpy as np

ap = argparse.ArgumentParser()
ap.add_argument("--module", default="PUCC_X2"); ap.add_argument("--libdir", default="../../SUCCSD_new")
ap.add_argument("--fc", required=True); ap.add_argument("--mo", required=True); ap.add_argument("--t1t2", required=True)
ap.add_argument("--newmo", default=None)
ap.add_argument("--moe", default="fock", help="fock | pickle of NSO energies | 'eff' (5th entry of the --newmo/--mo pickle)")
ap.add_argument("--sp", type=int, default=2); ap.add_argument("--grid", type=int, nargs=2, default=[1, 8])
ap.add_argument("--J", type=float, default=None); ap.add_argument("--inc", action="store_true")
ap.add_argument("--epcc", type=float, default=None); ap.add_argument("--tag", default="")
ap.add_argument("--pg", type=int, default=0)
args = ap.parse_args()
for k in ("fc", "mo", "t1t2", "libdir"):
    setattr(args, k, os.path.abspath(getattr(args, k)))
if args.newmo: args.newmo = os.path.abspath(args.newmo)
if args.moe not in ("fock", "eff"): args.moe = os.path.abspath(args.moe)
root = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))   # PHF_PCC repo root (this file: SUCCSD_new/)
sys.path.insert(0, root); sys.path.insert(0, args.libdir)
from ao2mo import ao2mo
from Spin_Proj import Spin_Proj
from lean_ints import unpack_eri, ghf_fock_ao, dirac_mo_eri_F
mod = __import__(args.module)

_fc = pickle.load(open(args.fc, "rb"))
Enuc, h1, eri, ncipg, weightpg, Rpg, NAO, NOcc = _fc[:8]
Spin2 = int(_fc[8]) if len(_fc) > 8 else 0
if args.J is None: args.J = Spin2 / 2
NAO = int(NAO); NOcc = int(NOcc); NSO = 2 * NAO; o = NOcc; v = NSO - NOcc
d = pickle.load(open(args.mo, "rb"))
if isinstance(d, (list, tuple)) and len(d) >= 4:
    MO, fsp, fpg, fk = d[:4]
else:
    MO = d[0] if isinstance(d, (list, tuple)) else d; fsp = None
MO = np.array(MO, dtype=complex)
T1, T2 = pickle.load(open(args.t1t2, "rb"))
T1 = np.array(T1, dtype=complex); T2 = np.array(T2, dtype=complex)
assert T1.shape == (v, o) and T2.shape == (v, v, o, o), (T1.shape, T2.shape)

# optional rotation of T into a second basis (same occupied and virtual spaces)
if args.newmo:
    dn = pickle.load(open(args.newmo, "rb"))
    MOn = np.array(dn[0] if isinstance(dn, (list, tuple)) else dn, dtype=complex)
    moe_eff = np.real(dn[4]) if (isinstance(dn, (list, tuple)) and len(dn) >= 5) else None
    Uo = MO[:, :o].conj().T @ MOn[:, :o]          # new_i' = sum_i old_i Uo(i,i')
    Uv = MO[:, o:].conj().T @ MOn[:, o:]
    err = max(np.max(np.abs(Uo.conj().T @ Uo - np.eye(o))), np.max(np.abs(Uv.conj().T @ Uv - np.eye(v))))
    print("basis rotation: max deviation from unitarity %.2e (occupied/virtual spaces must coincide)" % err)
    if err > 1e-6: sys.exit("--newmo does not span the same occupied/virtual spaces as --mo")
    T1 = Uv.conj().T @ T1 @ Uo
    T2 = np.einsum("aA,bB,abij,iI,jJ->ABIJ", Uv.conj(), Uv.conj(), T2, Uo, Uo, optimize=True)
    MO = MOn

cwd0 = os.getcwd(); os.chdir(root)
ncisp, roota, rootb, rooty, weightsp, R1, R2 = Spin_Proj(args.sp, args.grid, NAO, NSO, args.J)
os.chdir(cwd0)
if fsp is None: fsp = np.ones(ncisp, dtype=complex)
if args.pg:
    Rpg = np.array(Rpg, dtype=complex); weightpg = np.array(weightpg); fpg = np.ones(int(ncipg), dtype=complex)
else:
    fpg = np.ones(1, dtype=complex)
    weightpg = np.ones([1, 1, 1]); Rpg = np.zeros([1, NSO, NSO], dtype=complex); Rpg[0] = np.eye(NSO)

eri4 = unpack_eri(eri, NAO)
HOne = np.zeros([NSO, NSO]); HOne[:NAO, :NAO] = h1; HOne[NAO:, NAO:] = h1
H1 = ao2mo(HOne, MO, 2)
H2 = dirac_mo_eri_F(eri4, MO, NAO)
Fmo = MO.conj().T @ ghf_fock_ao(h1, eri4, MO[:, :o], NAO) @ MO
if args.moe == "fock":
    moe = np.real(np.diag(Fmo)).copy()
    print("max |off-diagonal oo/vv Fock| in this basis: %.2e" % max(np.max(np.abs(Fmo[:o, :o] - np.diag(np.diag(Fmo[:o, :o])))),
                                                                     np.max(np.abs(Fmo[o:, o:] - np.diag(np.diag(Fmo[o:, o:]))))))
elif args.moe == "eff":
    if not args.newmo: moe_eff = np.real(d[4]) if (isinstance(d, (list, tuple)) and len(d) >= 5) else None
    if moe_eff is None: sys.exit("--moe eff needs a 5-entry pickle [MO,fsp,fpg,fk,moe] (PHFFock.py output or moe_feff.py)")
    moe = moe_eff.copy()
else:
    moe = np.real(np.asarray(pickle.load(open(args.moe, "rb")))).copy()
    assert moe.shape == (NSO,)
del eri4
for i in range(args.grid[0]): R1[i] = ao2mo(R1[i], MO, 2)
for i in range(args.grid[1]): R2[i] = ao2mo(R2[i], MO, 2)
for i in range(len(Rpg)): Rpg[i] = ao2mo(Rpg[i], MO, 2)

from mpi4py import MPI
comm = MPI.COMM_WORLD; mpirank = comm.Get_rank(); mpisize = comm.Get_size()
tfun = mod.pgcc_t
et_coh, et_inc, s00, devz = tfun(H1, H2, NAO, NOcc, args.J, args.sp, args.grid, R1, R2, Rpg,
                                       roota, rootb, rooty, weightsp, weightpg, fsp, fpg, T1, T2, moe, int(args.inc),
                                       comm.py2f())
res = {"mpi_ranks": mpisize, "ET_coherent": complex(et_coh).real, "ET_coherent_imag": complex(et_coh).imag,
       "ET_incoherent": complex(et_inc).real if args.inc else None, "S00": complex(s00).real, "S00_imag": complex(s00).imag,
       "max_Zbar_minus_T": float(devz), "moe": args.moe if args.moe in ("fock","eff") else os.path.basename(args.moe),
       "newmo": os.path.basename(args.newmo) if args.newmo else None, "sp": args.sp, "grid": args.grid, "pg": args.pg, "npg": len(Rpg)}
if args.epcc is not None:
    res["E_PCCSD"] = args.epcc; res["E_PCCSD_T_coherent"] = args.epcc + res["ET_coherent"]
    if args.inc: res["E_PCCSD_T_incoherent"] = args.epcc + res["ET_incoherent"]
out = os.path.join(os.path.dirname(args.t1t2), "result_T%s.json" % args.tag)
if mpirank == 0:
    json.dump(res, open(out, "w"), indent=1)
    print("RESULT_T", json.dumps(res), flush=True)
