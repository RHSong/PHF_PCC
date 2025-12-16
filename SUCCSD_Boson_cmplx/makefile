#!/bin/sh
#mpif90 -fPIC -c -g -fopenmp -O3 Precision.f90 Constants.f90 CCRes_dr.f90 UtildeAlg_dr.f90 IntTrans.f90 Wrappers.f90
#mpif90 -fPIC -c -g -fopenmp -O3 outmat.f90 AllCase.f90 Functions.f90 HamiltonianData.f90 IO.f90 ObaraSaika.f90 DoIntegrals.f90 DIIS.f90
#mpif90 -fPIC -c -g -fopenmp -O3 UHF.f90 UHFWrap.f90 GHF.f90 MP2.f90 DIISCC.f90 DIISPCC.f90 CCScr.f90 CCSD_DIIS.f90 RootsWeights.f90 RotationMatrices.f90
#mpif90 -fPIC -c -g -fopenmp -O3 CCRes.f90 PCCRes.f90 BuildSC_dr.f90
#mpif90 -fPIC -c -g -fopenmp -O3 PCCCI.f90 IndexCI.f90 Hamiltonian.f90 PHFTools.f90 CI.f90 CITools.f90 UtildeAlg.f90 FPUCC_Tools.f90 FPUCC_FCI_Tools.f90 InvSy.f90 FPUCC.f90 Broyden.f90 Spin.f90
#mpif90 -shared *.o -o pucc_X2.so -L . -l AllLinAlg
#f2py --f90exec=mpif90 --opt='-O3' --f90flags="-fopenmp" -lgomp -m PUCC_X2 -c Precision.f90 Constants.f90 Main.f90 pucc_X2.so
#
#rm *mod
#rm *.o

#!/bin/sh
mpiifort -fPIC -c -g -qopenmp -heap-arrays -O3 Precision.f90 Constants.f90
mpiifort -fPIC -c -g -qopenmp -heap-arrays -O3 Wrappers.f90 DIIS.f90 UHF.f90 IntTrans.f90
mpiifort -fPIC -c -g -qopenmp -heap-arrays -O3 CCRes.f90 PCCCI.f90 CCSD_T.f90
mpiifort -fPIC -c -g -qopenmp -heap-arrays -O3 FPUCC_Tools.f90 FPUCC.f90 Broyden.f90 Spin.f90
mpiifort -shared *.o -o pucc_X2.so -qmkl
f2py --f90exec=mpiifort --f77exec=mpiifort --opt='-O3' --f90flags="-qopenmp" -lgomp -m PUCC_X2 -c Precision.f90 Constants.f90 Main.f90 pucc_X2.so

#rm *mod
#rm *.o
