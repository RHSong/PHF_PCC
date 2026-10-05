      Subroutine PGCC_T(H1,H2,NAO,NSO,NOcc,QJ,SP,NPoints,Npg,ncisp,ncipg,R1,R2,Rpg, &
                        Roota,Rootb,Rooty,Weightsp,Weightpg,fsp,fpg,T1,T2,MOE,DoInc, &
                        ETcoh,ETinc,S00,DevZ,comm)
!  f2py entry: coherent (T) from saved PCC amplitudes (see FPUCC::CoherentT).
!  DoInc /= 0 also evaluates the old per-grid (T) sum for comparison.
      Use Precision
      Use FPUCC
      Implicit None
      Integer, Intent(In) :: NAO, NSO, NOcc, QJ, SP, NPoints(2), Npg, ncisp, ncipg, DoInc, comm
      Complex (Kind=pr), Intent(In) :: H1(NSO,NSO), H2(NSO,NSO,NSO,NSO)
      Complex (Kind=pr), Intent(In) :: R1(:,:,:), R2(:,:,:), Rpg(Npg,NSO,NSO)
      Real (Kind=pr), Intent(In) :: Roota(:), Rootb(:), Rooty(:)
      Real (Kind=pr), Intent(In) :: Weightsp(:,:), Weightpg(Npg,ncipg,ncipg)
      Complex (Kind=pr), Intent(In) :: fsp(ncisp), fpg(ncipg)
      Complex (Kind=pr), Intent(In) :: T1(NOcc+1:NSO,NOcc), T2(NOcc+1:NSO,NOcc+1:NSO,NOcc,NOcc)
      Real (Kind=pr), Intent(In) :: MOE(NSO)
      Complex (Kind=pr), Intent(Out) :: ETcoh, ETinc, S00
      Real (Kind=pr), Intent(Out) :: DevZ
      Call CoherentT(H1,H2,T1,T2,MOE,NOcc,NAO,NSO,QJ,SP,NPoints,Npg,ncisp,ncipg,R1,R2,Rpg, &
                     Roota,Rootb,Rooty,Weightsp,Weightpg,fsp,fpg,(DoInc /= 0),ETcoh,ETinc,S00,DevZ,comm)
      End Subroutine PGCC_T
