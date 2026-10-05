# SUCCSD_new — OpenMP-only restructuring of SUCCSD_Boson_cmplx

Same algorithm and same numbers as `SUCCSD_base` (PGCC interface unchanged:
`from PUCC_X2 import pgcc`), but the per-grid-point kernel no longer builds
any NSO^4 or v^4 intermediate.  Verified on N2 (NSO=56, NOcc=14, SP=2,
grid 1x8, one Broyden iteration): energies agree to 1e-13, T1/T2 to 1e-15,
for CCSD, CCSD(T) (TrunODE=3) and the CmplxConj=1 branch.  See `tests/`.

## Build
    sh build.sh          # bundled libAllLinAlg.so BLAS (reference, single-threaded)
    sh build.sh mkl      # MKL ZGEMM (needs LD_LIBRARY_PATH to include the MKL lib dir)
At run time `LD_LIBRARY_PATH` must contain this directory (pucc_X2.so, libAllLinAlg.so).

## Test knobs (environment variables, read once at start of PGCC / BroydenIter)
    PCC_CYCMAX=n     stop after n+1 Broyden iterations (default 300)
    PCC_TRUNODE=3    switch on the (T) kernel (Main.f90 default is 1)

## What changed (file by file)
* `ERIBlocks.f90` (new): the o/v blocks of <pq||rs> that the CC kernel needs
  (oooo, ooov, oovv, ovvo, ovoo, vovv, vvvo, vvoo and three permuted copies;
  plus two (T) blocks when TrunODE=3) are extracted ONCE per PbarHbarOlap call
  instead of being sliced from the NSO^4 array at every grid point.  The vvvv
  block is never stored (largest block is o v^3).  `Conj=.True.` builds the
  blocks of conjg(H2) element-wise - no conjg(HTwo) NSO^4 temporary.
* `CCRes.f90`: `CCRes12` is a Stanton-Gauss CCSD residual on blocks.  The
  particle-particle ladder reads <ab||ef> straight from H2 with a strided
  leading dimension (one ZGEMM per b, threaded over b), so no v^4 Wvvvv/itm3
  copies.  Ring terms via ZGEMM on permuted (v o) x (v o) copies instead of
  scalar loops.  `CCEnergyBlocks` uses the oovv block.
* `IntTrans.f90`: `IntTran4` is four ZGEMM passes; new `TransBlockT2` forms a
  single block of U T2 U^-1 directly (bra indices first).
* `PCCCI.f90`: `TransT2` (build W) never forms the NSO^4 transformed tensor;
  all blocks of the second-order term come from `TransBlockT2`, and the vvvv
  piece is done "contract first" (O(o^2 v^3) instead of O(o^2 v^4 + v^5)).
* `FPUCC_Tools.f90` / `FPUCC.f90`: `BuildKernSD`, `EvalKernel`, `EvalOvlp`,
  `BuildSC_test` carry doubles as (v,v,o,o) only (old code used NSO^4 arrays
  eps2t/lambda2t/a2/C2V1).  `MKFock_MO` and ESCF are computed once per
  PbarHbarOlap instead of once per grid point.
* `CCSD_T.f90`: same W/V formula, but with pre-permuted integral blocks
  (`Tvvvo`, `Tovoo`, `Vvvoo`) and pre-permuted T2, so every operand is a
  contiguous/strided ZGEMM or ZGERU argument; the a>=b>=c triangle is
  flattened and scheduled dynamically; no nested OpenMP.
* `Wrappers.f90`: `PZGEMM` (column-split threaded ZGEMM for the bundled BLAS).
* `Broyden.f90`, `Main.f90`: `ShutDownERIBlocks`; PCC_CYCMAX / PCC_TRUNODE.
* `Constants.f90`: `ZZero`.

## Measured on N2 (4 threads, one Broyden iteration = 16 kernel evaluations)
| run     | base (bundled BLAS) | new (bundled BLAS) | new (MKL) | per grid TransT2 (base/new/mkl) | per grid CC (base/new/mkl) |
|---------|---------------------|--------------------|-----------|---------------------------------|----------------------------|
| CCSD    | 117.5 s             | 28.1 s             | 3.0 s     | 2.62 / 0.43 / 0.07 s            | 1.98 / 1.27 / 0.09 s       |
| CCSD(T) | 298.7 s             | 105.8 s            | 9.0 s     | 2.63 / 0.42 / 0.07 s            | 13.3 / 6.1 / 0.46 s        |
Peak RSS (Python + Fortran): base 1.18 GB, new 0.71 GB (H2 itself is 157 MB for N2).
The bundled libAllLinAlg.so BLAS is a single-threaded reference implementation
(~2 GFlop/s); with the kernel now expressed entirely in ZGEMM/ZGERU, linking MKL
(`sh build.sh mkl`) gives the remaining factor of ~10.


### SGCCSD(T) on N2 (SP=1, grid 14x8 = 224 kernel evaluations, one Broyden iteration, 4 threads)
| code               | wall      | peak RSS | per grid: TransT2 | CC+(T)  |
|--------------------|-----------|----------|-------------------|---------|
| old, bundled BLAS  | 1 h 15 min| 1.17 GB  | 2.61 s            | 14.9 s  |
| new, bundled BLAS  | 27.8 min  | 0.72 GB  | 0.42 s            | 7.0 s   |
| new, MKL           | 2.1 min   | 0.74 GB  | 0.07 s            | 0.47 s  |
E_PCC = -109.1337024057189 for all three; T1/T2 agree with the old code to 2e-14.

## Coherent (T) for PCC (new, 2026-10-01)
Theory: first-order Löwdin/Lagrangian triples correction with the projector metric,
closed with the Hermitian (Λ ≈ T†) approximation and lowest order in H_N.  Because
W and V are linear in the amplitudes the grid sum acts on the amplitudes:
    Zbar_n = Σ_g W_g e^{z0(g)} Z_n(g) / S00,   S00 = Σ_g W_g e^{z0(g)}
    E(T)   = standard (T) formula on (Zbar1, Zbar2), denominators from MOE
(`FPUCC::CoherentT`, f2py entry `pgcc_t` in MainT.f90, driver tests/pcc_t.py).
MOE can be the bare semicanonical Fock diagonal (`--moe fock`) or the PHF
effective-Fock eigenvalues (`--moe eff`, 5th entry of the PHFFock.py pickle, with
`--newmo` the orbitals that diagonalize its oo/vv blocks; T is rotated into them).
`--inc` also returns the old per-grid sum Σ_g W_g e^{z0} E_T(Z(g))/S00.
Cost: one (T) per PCC solution instead of one per grid point.  SP=0 reduces to CCSD(T)
(checked: coherent = per-grid = pyscf UCCSD(T) to 1e-5, limited by PCC convergence).
PCCSD+T: the same driver with the projector off (`--sp 0 --grid 1 1`): standard (T) on the
PHF determinant with the converged PCC amplitudes.  Note S00 includes |fsp|^2 (the CI
coefficient carried in the orbital pickle); it cancels in every reported quantity.

N2 / cc-pVDZ / 3.0 bohr, all-electron, SUHF reference (tests_N2/R3.0_bohr_ae), FCI -109.088876:
| method                              | E(T)       | total          | vs FCI   |
|-------------------------------------|------------|----------------|----------|
| UCCSD / UCCSD(T) (pyscf, UHF ref)   | -0.015984  | -109.068864    | +20.0 mHa|
| SUCCSD (grid 1x8 = 1x16)            |            | -109.077871    | +11.0 mHa|
| SUCCSD(T) coherent, bare Fock denom | -0.015527  | -109.093397    | -4.5 mHa |
| SUCCSD(T) coherent, PHF-Fock denom  | -0.015536  | -109.093407    | -4.5 mHa |
| SUCCSD(T) old per-grid sum          | -0.026075  | -109.103946    | -15.1 mHa|
| PCCSD+T, bare Fock denom            | -0.008134  | -109.086004    | +2.9 mHa |
| PCCSD+T, PHF-Fock denom             | -0.008112  | -109.085983    | +2.9 mHa |

## Memory estimate for Cr2 (NSO=136, NOcc=12, v=124; H2 = 5.5 GB complex)
Fortran side of the new code: ERIBlocks ~0.8 GB (1.1 GB with (T)), TransT2
scratch ~1.5 GB (a few o v^3 blocks), CC residual < 0.2 GB; no NSO^4 or v^4
temporaries, so roughly H2 + 3 GB.  The old code holds several NSO^4
complex temporaries (conjg(HTwo), IntTran4Sim scratch, U/C2, C2V1 ...),
i.e. 5.5 GB each, which is where the ~60 GB came from.
NOTE the Python driver still dominates for Cr2: `Main_PCC.py` builds the
real NSO^4 HTwo (2.7 GB), `ao2mo(...,4)` via einsum makes complex NSO^4
intermediates, and f2py copies H2 to Fortran order unless it is passed as
`np.asfortranarray(H2)`.  Use the lean integral path of `Cr2_new/Cr2_int.py`
and pass a Fortran-ordered H2 to keep the Python peak near one copy of H2.
