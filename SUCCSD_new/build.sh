#!/bin/sh
# Build the OpenMP-only PCC library (pucc_X2.so) and the f2py module PUCC_X2.
#   sh build.sh          -> BLAS/LAPACK from the bundled libAllLinAlg.so (reference, single-threaded)
#   sh build.sh mkl      -> ZGEMM etc. from Intel MKL (gfortran interface, GNU OpenMP threading);
#                           libAllLinAlg.so is still linked for the routines MKL does not provide.
# Run inside SUCCSD_new; needs mpif90, f2py (SETUPTOOLS_USE_DISTUTILS=stdlib for new setuptools).
set -e
export SETUPTOOLS_USE_DISTUTILS=stdlib
MKLROOT=${MKLROOT:-/home/song/intel/oneapi/mkl/2023.0.0}
if [ "$1" = "mkl" ]; then
  BLASLIBS="-L$MKLROOT/lib/intel64 -lmkl_gf_lp64 -lmkl_gnu_thread -lmkl_core -L. -lAllLinAlg -lgomp -lpthread -lm -ldl"
else
  BLASLIBS="-L. -lAllLinAlg -lgomp"
fi
mkdir -p build && cd build && rm -f *.o *.mod
mpif90 -fPIC -c -fopenmp -O3 ../Precision.f90 ../Constants.f90 ../Wrappers.f90 ../DIIS.f90 ../UHF.f90 \
   ../IntTrans.f90 ../ERIBlocks.f90 ../CCRes.f90 ../PCCCI.f90 ../CCSD_T.f90 ../FPUCC_Tools.f90 \
   ../FPUCC.f90 ../Broyden.f90 ../Spin.f90
cd ..
mpif90 -shared build/*.o -o pucc_X2.so $BLASLIBS
f2py --f90exec=mpif90 --opt='-O3' --f90flags="-fopenmp -Ibuild" -lgomp -m PUCC_X2 -c Precision.f90 Constants.f90 Main.f90 MainT.f90 pucc_X2.so > build/f2py.log 2>&1
echo "built PUCC_X2 ($([ "$1" = mkl ] && echo MKL || echo bundled BLAS)); at run time LD_LIBRARY_PATH must contain $(pwd)$([ "$1" = mkl ] && echo " and $MKLROOT/lib/intel64")"
