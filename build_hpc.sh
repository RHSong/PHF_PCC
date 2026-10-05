#!/bin/bash
# Build PHFTools_F, FockTools_F (PHF_fort) and PUCC_X2 (SUCCSD_new) on UCI HPC.
# Intel 2021.4 ifort/mpiifort + Intel MPI + MKL (replaces the bundled libAllLinAlg.so),
# numpy 2.x f2py with the meson backend.  Run from the PHF_PCC root, ideally in a Slurm job.
set -euo pipefail
module load python/3.10.18
module load intel/2021.4 intelmpi/2021.4 mkl
export PATH=$HOME/.local/bin:$PATH                     # meson, ninja (pip --user)
ROOT=$(pwd)
ILIB=/sopt/IntelOneAPI/2021.4/compiler/latest/linux/compiler/lib/intel64_lin
export CC=/sopt/GCC/9.4.0/bin/gcc                      # icc cannot compile on this system
export FC=mpiifort
export LDFLAGS="-L$ILIB -L$MKLROOT/lib/intel64 -L$I_MPI_ROOT/lib/release -L$I_MPI_ROOT/lib"
LIBS="-lmkl_rt -liomp5 -lmpifort -lmpi"
F2PY="python3 -m numpy.f2py --backend meson"

echo "=== PHF_fort: PHFTools_F"
cd $ROOT/PHF_fort && rm -rf bld_* *.so
FFLAGS="-qopenmp -O3 -fPIC" $F2PY -c Precision.f90 Constants.f90 BuildSC_dr.f90 Wrappers.f90 PHFTools.f90 \
    -m PHFTools_F --build-dir bld_phf $LIBS > bld_phf.log 2>&1 || { tail -30 bld_phf.log; exit 1; }
echo "=== PHF_fort: FockTools_F"
FFLAGS="-qopenmp -O3 -fPIC" $F2PY -c Precision.f90 Constants.f90 Wrappers.f90 FockTools.f90 \
    -m FockTools_F --build-dir bld_fock $LIBS > bld_fock.log 2>&1 || { tail -30 bld_fock.log; exit 1; }
cp PHFTools_F*.so FockTools_F*.so $ROOT/

echo "=== SUCCSD_new: libpucc_X2.so + PUCC_X2"
cd $ROOT/SUCCSD_new && mkdir -p build && rm -f build/*.o build/*.mod libpucc_X2.so PUCC_X2*.so
( cd build && mpiifort -fPIC -c -qopenmp -O3 -heap-arrays ../Precision.f90 ../Constants.f90 ../Wrappers.f90 ../DIIS.f90 ../UHF.f90 \
    ../IntTrans.f90 ../ERIBlocks.f90 ../CCRes.f90 ../PCCCI.f90 ../CCSD_T.f90 ../FPUCC_Tools.f90 ../FPUCC.f90 \
    ../Broyden.f90 ../Spin.f90 )
mpiifort -shared -qopenmp build/*.o -o libpucc_X2.so -qmkl
FFLAGS="-qopenmp -O3 -fPIC -heap-arrays -I$ROOT/SUCCSD_new/build" $F2PY -c Precision.f90 Constants.f90 Main.f90 MainT.f90 \
    -m PUCC_X2 --build-dir bld_f2py -L$ROOT/SUCCSD_new -lpucc_X2 $LIBS > bld_f2py.log 2>&1 || { tail -30 bld_f2py.log; exit 1; }
echo "=== import check"
cd $ROOT && export LD_LIBRARY_PATH=$ROOT/SUCCSD_new:$LD_LIBRARY_PATH
env -u I_MPI_PMI_LIBRARY python3 -c "
import sys; sys.path.insert(0,'SUCCSD_new')
import PHFTools_F, FockTools_F, PUCC_X2
print('PHFTools_F :', [n for n in dir(PHFTools_F) if not n.startswith('_')])
print('FockTools_F:', [n for n in dir(FockTools_F) if not n.startswith('_')])
print('PUCC_X2    :', [n for n in dir(PUCC_X2) if not n.startswith('_')])"
echo "=== build done"
