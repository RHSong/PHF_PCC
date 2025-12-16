#!/bin/sh
#f2py --fcompiler=pg --f90flags="-mp" -lgomp -m PUCC_X2 -c Precision.f90 Constants.f90 Main.f90 pucc_X2.so
#f2py --fcompiler=gfortran --f90flags="-fopenmp" -lgomp -m PHFTools -c Precision.f90 Constants.f90 BuildSC_dr.f90 Wrappers.f90 PHFTools.f90 -L . -lAllLinAlg
#f2py --f90exec=mpif90 --fcompiler=gfortran --f90flags="-fopenmp" --opt='-O1' -lgomp -m PHFTools_F -c Precision.f90 Constants.f90 BuildSC_dr.f90 Wrappers.f90 PHFTools.f90 -L. -lAllLinAlg

#f2py --f90exec=mpif90 --fcompiler=intelem --f90flags="-qopenmp" --opt='-O1' -lgomp -m PHFTools_F -c Precision.f90 Constants.f90 BuildSC_dr.f90 Wrappers.f90 PHFTools.f90 -L. -lAllLinAlg
f2py --f90exec=mpif90 --opt='-O3' --f90flags="-fopenmp" -lgomp -m PHFTools_F -c Precision.f90 Constants.f90 BuildSC_dr.f90 Wrappers.f90 PHFTools.f90 -L. -lAllLinAlg

#rm *mod
#rm *.o
