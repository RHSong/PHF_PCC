      Subroutine PGCC_T(H1,H2,NAO,NSO,NOcc,QJ,SP,NPoints,Npg,ncisp,ncipg,R1,R2,Rpg, &
                        Roota,Rootb,Rooty,Weightsp,Weightpg,fsp,fpg,T1,T2,MOE,DoInc, &
                        ETcoh,ETinc,S00,DevZ,comm)
!  f2py entry: coherent (T) from saved PCC amplitudes (see FPUCC::CoherentT).
!  DoInc /= 0 also evaluates the old per-grid (T) sum for comparison.
      Use Precision
      Use FPUCC
      Use FPUCC_Tools, Only: SetSzAlpha
      Implicit None
      Integer, Intent(In) :: NAO, NSO, NOcc, SP, NPoints(2), Npg, ncisp, ncipg, DoInc, comm
      Real(kind=pr), Intent(In) :: QJ   ! total spin J, integer or half-integer
      Complex (Kind=pr), Intent(In) :: H1(NSO,NSO), H2(NSO,NSO,NSO,NSO)
      Complex (Kind=pr), Intent(In) :: R1(:,:,:), R2(:,:,:), Rpg(Npg,NSO,NSO)
      Real (Kind=pr), Intent(In) :: Roota(:), Rootb(:), Rooty(:)
      Real (Kind=pr), Intent(In) :: Weightsp(:,:), Weightpg(Npg,ncipg,ncipg)
      Complex (Kind=pr), Intent(In) :: fsp(ncisp), fpg(ncipg)
      Complex (Kind=pr), Intent(In) :: T1(NOcc+1:NSO,NOcc), T2(NOcc+1:NSO,NOcc+1:NSO,NOcc,NOcc)
      Real (Kind=pr), Intent(In) :: MOE(NSO)
      Complex (Kind=pr), Intent(Out) :: ETcoh, ETinc, S00
      Real (Kind=pr), Intent(Out) :: DevZ
! n_alpha for the SUHF Sz-block projections: SUHF needs J = Sz = (n_alpha - n_beta)/2
       If (SP == 2) Then
         Call SetSzAlpha((NOcc + nint(2*QJ)) / 2)
       Else
         Call SetSzAlpha(NOcc / 2)
       End If
      Call CoherentT(H1,H2,T1,T2,MOE,NOcc,NAO,NSO,QJ,SP,NPoints,Npg,ncisp,ncipg,R1,R2,Rpg, &
                     Roota,Rootb,Rooty,Weightsp,Weightpg,fsp,fpg,(DoInc /= 0),ETcoh,ETinc,S00,DevZ,comm)
      End Subroutine PGCC_T


      Subroutine PGCC_T_SAFE(H1,H2,NAO,NSO,NOcc,QJ,SP,NPoints,Npg,ncisp,ncipg,R1,R2,Rpg, &
                        Roota,Rootb,Rooty,Weightsp,Weightpg,fsp,fpg,T1,T2,MOE,DoInc, &
                        ETcoh,ETinc,S00,DevZ,comm)
!  f2py entry: coherent (T) with the zero-overlap-safe Thouless step (FPUCC::CoherentTSafe; tau from PCC_ZO_TAU).
!  DoInc /= 0 also evaluates the old per-grid (T) sum for comparison.
      Use Precision
      Use FPUCC
      Use FPUCC_Tools, Only: SetSzAlpha
      Implicit None
      Integer, Intent(In) :: NAO, NSO, NOcc, SP, NPoints(2), Npg, ncisp, ncipg, DoInc, comm
      Real(kind=pr), Intent(In) :: QJ   ! total spin J, integer or half-integer
      Complex (Kind=pr), Intent(In) :: H1(NSO,NSO), H2(NSO,NSO,NSO,NSO)
      Complex (Kind=pr), Intent(In) :: R1(:,:,:), R2(:,:,:), Rpg(Npg,NSO,NSO)
      Real (Kind=pr), Intent(In) :: Roota(:), Rootb(:), Rooty(:)
      Real (Kind=pr), Intent(In) :: Weightsp(:,:), Weightpg(Npg,ncipg,ncipg)
      Complex (Kind=pr), Intent(In) :: fsp(ncisp), fpg(ncipg)
      Complex (Kind=pr), Intent(In) :: T1(NOcc+1:NSO,NOcc), T2(NOcc+1:NSO,NOcc+1:NSO,NOcc,NOcc)
      Real (Kind=pr), Intent(In) :: MOE(NSO)
      Complex (Kind=pr), Intent(Out) :: ETcoh, ETinc, S00
      Real (Kind=pr), Intent(Out) :: DevZ
! n_alpha for the SUHF Sz-block projections: SUHF needs J = Sz = (n_alpha - n_beta)/2
       If (SP == 2) Then
         Call SetSzAlpha((NOcc + nint(2*QJ)) / 2)
       Else
         Call SetSzAlpha(NOcc / 2)
       End If
      Call CoherentTSafe(H1,H2,T1,T2,MOE,NOcc,NAO,NSO,QJ,SP,NPoints,Npg,ncisp,ncipg,R1,R2,Rpg, &
                     Roota,Rootb,Rooty,Weightsp,Weightpg,fsp,fpg,(DoInc /= 0),ETcoh,ETinc,S00,DevZ,comm)
      End Subroutine PGCC_T_SAFE
