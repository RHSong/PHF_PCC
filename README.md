Projected HF & CC

## Workflow

| step | script | output |
|---|---|---|
| integrals, frozen core | `Main_FC.py` (settings in `parameter.py`, `NFC`) | `FC.p`, `FCMO.p` |
| RHF / UHF / GHF | `Main_HF.py` (settings in `parameter_fc.py`) | `GHFMO.p`, density guesses |
| PHF, gradient-based (BFGS, basin hopping) | `Main_PHF.py` -> `PHF.py` + `PHF_fort/PHFTools.f90` | `SGHF.p`, `newsol.p` |
| PHF, Fock-based (DIIS) | `PHFFock.py [guess.p]` -> `PHF_fort/FockTools.f90` | `SGHFFock.p`, `newmo.p`, `newmo_best.p` |
| projected CC | `Main_PCC.py` -> `SUCCSD_new/` (module `PUCC_X2`) | energy, T1/T2 |

Spin projection `SP`: 0 none, 1 SGHF (full SU(2) grid, `ngrid = [Lebedev, gamma]`), 2 SUHF (`[1, n]`), 3 Sz.
Run the scripts from the calculation folder so that its `parameter.py` / `parameter_fc.py` are imported.

## Update: Fock-based PHF, optimized PCC

### Fock-based PHF (`PHFFock.py`, `PHF_fort/FockTools.f90`)
* PHF from the PHF effective Fock matrix F = F0 + F+ + F- (PHFB formulation); residual [F, rho];
  Roothaan steps with DIIS (`use_diis`, `ndiis`, `diis_start`), optional level shift,
  frozen occupied orbitals (`NFO`), symmetry filtering, aufbau occupation.
* Reproduces the gradient-based PHF energies (N2 SUHF and SGHF) to within 5e-7 Ha.
* Every iteration writes `newmo.p` (latest orbitals, like `newsol.p`) and `newmo_best.p`
  (lowest energy evaluated so far; with DIIS the energy need not decrease monotonically).
  Output: `SGHFFock.p = [MO, fsp, fpg, fk, moe]`, with `moe` the eigenvalues of the PHF Fock matrix.
  Readers of `SGHF.p`-type files take the first four entries.
* Grid loop parallelised with OpenMP over all (Lebedev, gamma) points (`collapse(2)`).
* Build: `PHF_fort/makefile` (FockTools line), then copy `FockTools_F*.so` to the repo root.

### PHF utilities
* `FixGauge.py`: `FixSz(MOs)` applies the Sz boost exp(a Sz) with a chosen so that <Sz> = 0 of the
  determinant (the J = 0 projected energy is unchanged up to grid error). `FixGauge` minimises the
  determinant <S^2> over complexified spin rotations, which also sets <S> = 0.
* `Spin.py`: `calcS` used `BuildSz` for the y component; now `BuildSy`.
* `PHF.py`: `nhop` = number of basin-hopping iterations (0 = plain BFGS).

### Projected CC (`SUCCSD_new/`, f2py module `PUCC_X2`)
Restructured version of `SUCCSD_Boson_cmplx` with the same equations and results
(energies and T1/T2 agree with the original to 1e-13 / 1e-15 on N2, including (T) and CmplxConj = 1):
* no NSO^4 or v^4 intermediates in the kernels: o/v integral blocks extracted once (`ERIBlocks.f90`),
  vvvv ladder read directly from H2, build-W (`TransT2`) and the overlap residual (`EvalOvlp`)
  formed block by block, `IntTran4` and the CC residual as ZGEMM;
* MPI over the beta grid (kernels and overlap residual); OpenMP inside each rank;
* test knobs: `PCC_CYCMAX=n` (max Broyden cycles), `PCC_TRUNODE=3` (in-loop per-grid (T)).
* Build: `SUCCSD_new/build.sh [mkl]` locally, `build_hpc.sh` on UCI HPC (Intel 2021.4, Intel MPI,
  MKL, f2py meson backend, `-heap-arrays`). Details in `SUCCSD_new/README.md`.

**Timing: use the new code linked against MKL.** N2 / cc-pVDZ (NSO = 56, NOcc = 14), one Broyden
iteration, 4 OpenMP threads, same run in all three columns (energies and amplitudes identical):

| run | original code (bundled BLAS) | new code (bundled BLAS) | new code + MKL |
|---|---|---|---|
| SUCCSD, SP=2, 1x8 grid (16 kernel evaluations) | 117.5 s | 28.1 s | 3.0 s |
| SUCCSD with in-loop (T), SP=2, 1x8 grid | 298.7 s | 105.8 s | 9.0 s |
| SGCCSD with in-loop (T), SP=1, 14x8 grid (224 kernel evaluations) | 75.7 min | 27.8 min | 2.1 min |
| peak memory (Python + Fortran) | 1.17 GB | 0.72 GB | 0.74 GB |

The bundled `libAllLinAlg.so` is a single-threaded reference BLAS; the new kernels are written in
ZGEMM/ZGERU, so linking MKL (`sh build.sh mkl`, or `build_hpc.sh` on the cluster) gives a further
factor of about 10 (about 35 times faster than the original code). Larger case: Cr2 / cc-pVDZ-DK, NFC = 10 (NSO = 152, NOcc = 28), SGCCSD on the 14x6
grid converged in 34 iterations in 26.5 h with 2 MPI ranks x 8 threads (MKL, Haswell node),
19 GB per rank; the original code would need about 60 GB.

### Lean integrals (`lean_ints.py`)
Spin-orbital MO integrals <pq||rs> built as one complex Fortran-ordered NSO^4 array (pair-density
GEMMs, slab-wise antisymmetrisation), semicanonicalisation from the spatial integrals; used by
`Main_PCC.py` (also sets the Broyden history to 20). Cr2 NFC = 10 (NSO = 152): about 13 GB instead of
about 30 GB for the integral step.
