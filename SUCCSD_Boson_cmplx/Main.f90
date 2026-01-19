
      Subroutine PGCC(Ref,H1,H2,NAO,NSO,NOcc,QJ,CC,SP, &
	  			      NPoints, Npg, ncisp, ncipg, ncik, R1, R2, Rpg, Rk, &
                      Roota,Rootb,Rooty,Weightsp,Weightpg,fsp,fpg,fk, &
					  pBrdIn,X,Xinv,ENuc,Ene,S2,T1,T2,comm)
      Use Precision
      Use DIIS
      Use IntTrans
      Use FPUCC
      Use FPUCC_Tools
      Use Broyden
	  Use CCSD_T
      Use Spin
      Implicit None
! Dimensioning variables
      Integer, Intent(In) :: NOcc, NAO, NSO
      Integer, Intent(In) :: NPoints(2),pBrdIn, Npg
      Integer, Intent(In) :: QJ, CC,SP,comm
      Integer :: TrunODE = 1, NPtsPUCC(2)
	  Integer, Intent(In) :: ncisp, ncipg, ncik
      Real (Kind=pr), Intent(In) :: Roota(:), Rootb(:), Rooty(:)
      Real (Kind=pr), Intent(In) :: Weightsp(:,:), Weightpg(npg,ncipg,ncipg)
	  Complex (Kind=pr), Intent(In) :: fsp(ncisp),fpg(ncipg),fk(ncik)
	  Complex (Kind=pr), Intent(In) :: R1(:,:,:), R2(:,:,:)
	  Complex (Kind=pr), Intent(In) :: Rpg(npg,NSO,NSO)
	  Complex (Kind=pr), Intent(In) :: Rk(ncik,NSO,NSO)
! Integrals
      Complex (Kind=pr), Intent(In) :: H1(NSO,NSO), H2(NSO,NSO,NSO,NSO)
      Real (Kind=pr), Intent(In) :: ENuc
      Complex (Kind=pr), Allocatable :: OneH(:,:), ERI(:,:,:,:)
! SCF output
      Complex (Kind=pr), Intent(In) :: Ref(NSO,NSO)
      Complex (Kind=pr) :: Evecs(NSO,NSO)
      Complex (Kind=pr), Allocatable :: Fock(:,:)
      Real    (Kind=pr) :: ESCF, ECorr
      Complex (Kind=pr) :: ESCFc, Sz, ET
! Correlated stuff
      Complex (Kind=pr), Intent(Out) :: Ene, S2
      Complex (Kind=pr), Intent(Out) :: T1(NOcc+1:NSO,NOcc)
      Complex (Kind=pr), Intent(Out) :: T2(NOcc+1:NSO,NOcc+1:NSO,NOcc,NOcc)
      Complex (Kind=pr), Allocatable :: Z1(:,:), Z2(:,:,:,:)
      Real (Kind=pr) :: EneChk
      Real (Kind=pr), Allocatable :: T1Chk(:,:), T2Chk(:,:,:,:)
! Error checking variables
      Integer :: IAlloc
      Integer :: I, J, A, B
      Integer :: p, q, r, s, pA, qA, rA, sA, pB, qB, rB, sB 
      Logical :: Exists, ExistsA, ExistsB, ExistsC, ExistsD
      Logical :: DoCCDIn = .False., Restart = .True.
! Timers
      Integer         :: cnt1, cnt2, clock_rate, clock_max
      Real (Kind=pr)  :: cput1, cput2
! restore T1, T2
      Complex (Kind=pr), Intent(In) :: X(:,:), Xinv(:,:)
      Integer :: npgtmp
      Complex(kind=pr) :: Rpgtmp(1,NSO,NSO)
      Real (kind=pr) :: Weightpgtmp(1,ncipg,ncipg)
! check Fock
	  Complex(kind=pr) :: Foo(NOcc,NOcc), Fvv(NSO-NOcc,NSO-NOcc)

!==========================================!
!  basic information for the calculation.  !
!==========================================!
      Allocate(Fock(NSO,NSO), &
!               T1(NOcc+1:NSO,NOcc),            &
!               T2(NOcc+1:NSO,NOcc+1:NSO,NOcc,NOcc),                &
!               Z2(NOcc,NOcc,NOcc+1:NSO,NOcc+1:NSO),                &
               Stat=IAlloc)
      NPtsPUCC = NPoints
      Evecs = Ref
!      OneH(:,:) = H1(:,:)
!      ERI(:,:,:,:) = H2(:,:,:,:)

!      If (SP == 2) then
!        Call ReOrderMO(Evecs,OneH,ERI,NOcc,NAO,NSO)
!      EndIf
      npgtmp = 1
      Rpgtmp = Zero
      Do I = 1, NSO
        Rpgtmp(1,I,I) = One
      End Do
      Weightpgtmp = One


!=====================================!
!  Now the integral transformations.  !
!=====================================!
!      Call IntTran2(OneH,Evecs,NSO)
!      Call IntTran4(ERI,Evecs,NSO)
! Set up tools for spin projection
!      Call SetUpPUCCTools(Evecs,NAO,NSO)
!      Call SetUpFciPuccTools(NOcc,NAO,NSO)

      Call MKFock_MO(H1,H2,Fock,ESCFc,NOcc,NSO)
      ESCF = Real(ESCFc) + ENuc
      Write(*,*) "ESCF from MkFock", ESCFc

	  Foo = Zero
	  Do I = 1, NOcc
	  	Foo(I,I) = Fock(I,I)
	  End Do

	  Fvv = Zero
	  Do I = NOcc+1, NSO
	  	Fvv(I-NOcc,I-NOcc) = Fock(I,I)
	  End Do

	  Print *, "Test Foo", maxval(abs(Foo-Fock(1:NOcc,1:NOcc)))
	  Print *, "Test Fvv", maxval(abs(Fvv-Fock(NOcc+1:NSO,NOcc+1:NSO)))

!      Inquire(File="Chk_T1",Exist=ExistsA)
!      Inquire(File="Chk_T2",Exist=ExistsB)
!      If (ExistsA.and.ExistsB) then
!        Open(70,File="Chk_T1",form="unformatted")
!        Read(70) T1Chk
!        Close(70)
!        Open(70,File="Chk_T2",form="unformatted")
!        Read(70) T2Chk
!        Close(70)
!        Call ReadT1(T1,T1Chk,Evecs,X,NOcc,NSO)
!        Call ReadT2(T2,T2Chk,Evecs,X,NOcc,NSO)
!        Call DrvCCSD(Fock,ERI,T1,T2,NOcc,NSO,ESCF,ECorr,DoCCDIn)
!        Write(*,*) "Max T1", MaxVal(Real(T1,kind=pr))
!        Call SaveT1(T1,T1Chk,Evecs,Xinv,NOcc,NSO)
!        Call SaveT2(T2,T2Chk,Evecs,Xinv,NOcc,NSO)
!        Open(70,File="Chk_Ene",status="Replace",form="unformatted")
!        Write(70) Cmplx(ESCF+ECorr,Zero,Kind=pr)
!        Close(70)
!        Open(70,File="Chk_T1",status="Replace",form="unformatted")
!        Write(70) T1Chk
!        Close(70)
!        Open(70,File="Chk_T2",status="Replace",form="unformatted")
!        Write(70) T2Chk
!        Close(70)
!      Else
!        T1 = Zero
!        T2 = Zero
!        Call DrvCCSD(Fock,ERI,T1,T2,NOcc,NSO,ESCF,ECorr,DoCCDIn)
!        Write(*,*) "Max T1", MaxVal(Real(T1,kind=pr))
!        Call SaveT1(T1,T1Chk,Evecs,Xinv,NOcc,NSO)
!        Call SaveT2(T2,T2Chk,Evecs,Xinv,NOcc,NSO)
!        Open(70,File="Chk_Ene",status="Replace",form="unformatted")
!        Write(70) Cmplx(ESCF+ECorr,Zero,Kind=pr)
!        Close(70)
!        Open(70,File="Chk_T1",status="Replace",form="unformatted")
!        Write(70) T1Chk
!        Close(70)
!        Open(70,File="Chk_T2",status="Replace",form="unformatted")
!        Write(70) T2Chk
!        Close(70)
!      End If
!	  T1CC = T1
!	  T2CC = T2


!============================!
!  Do the Projected CCSD.    !
!============================!

      Call cpu_time(cput1)
      Call system_clock(cnt1,clock_rate,clock_max)

      Open(9,File="RKstat",Status="Replace")
      Write(9,*) ''
      Close(9)
!      Z1 = Zero
!      Z2 = Zero
!      Call PbarHbarEneLR(OneH,ERI,Evecs,Ene,T1,T2,Z1,Z2,NOcc,NAO,NSO,&
!                        ENuc,QJ,CIf,NPtsPUCC,Roota,Rootb,Rooty,Weights,TrunODE)
!      Print *, "PAV Ene=", Ene

      ExistsA = .True.
      Inquire(File="Chk_PUCC_Ene",Exist=Exists)
      ExistsA = Exists .and. ExistsA
      Inquire(File="Chk_PUCC_T1",Exist=Exists)
      ExistsA = Exists .and. ExistsA
      Inquire(File="Chk_PUCC_T2",Exist=Exists)
      ExistsA = Exists .and. ExistsA
      ExistsB = .True.
      If (ExistsA) then
!        Allocate(T1Chk(NOcc+1:NSO,NOcc), &
!        T2Chk(NOcc+1:NSO,NOcc+1:NSO,NOcc,NOcc))
        Open(70,File="Chk_PUCC_Ene",form="unformatted")
        Read(70) Ene
        Close(70)
        Open(70,File="Chk_PUCC_T1",form="unformatted")
        Read(70) T1
        Close(70)
        Open(70,File="Chk_PUCC_T2",form="unformatted")
        Read(70) T2
        Close(70)
!        Ene = Cmplx(EneChk,Zero)
!        T1 = Cmplx(T1Chk,Zero)
!        T2 = Cmplx(T2Chk,Zero)
!        Call ReadT1(T1,T1Chk,Evecs,X,NOcc,NSO)
!        Call ReadT2(T2,T2Chk,Evecs,X,NOcc,NSO)
!        Deallocate(T1Chk,T2Chk)
      Else
!       Ene = Zero 
        T1 = Zero
        T2 = Zero
!        Z1 = Zero
!        Z2 = Zero
        Ene = EScf
!        Call PbarHbarEneLR(OneH,ERI,Evecs,Ene,T1,T2,Z1,Z2,NOcc,NAO,NSO,ENuc,QJ,&
!                           SP,NPoints, Npg, ncisp, ncipg, ncik, R1, R2, Rpg, Rk, &
!						   Roota,Rootb,Rooty,Weightsp,Weightpg,fsp,fpg,fk,TrunODE,comm)
      EndIf


      Call BroydenIter(H1,H2,Evecs,Xinv,Ene,T1,T2,NOcc,NAO,NSO,ENuc,QJ, &
                       CC,SP,NPoints,Npg, ncisp, ncipg, ncik, R1, R2, Rpg, Rk, &
					   Roota,Rootb,Rooty,Weightsp,Weightpg, fsp,fpg,fk, &
					   TrunODE,pBrdIn,DoCCDIn,comm)


      Call cpu_time(cput2)
      Call system_clock(cnt2,clock_rate,clock_max)
      Write(*,*) 'VAP SUCCSD cpu time =  ', real(cput2-cput1)
      write(*,*) 'VAP SUCCSD real time = ', real(cnt2-cnt1)/real(clock_rate)
! calculate S^2
      S2 = zero
!!      Call CalcS2(Evecs,S2,Sz,T1,T2,NOcc,NAO,NSO,QJ,SP, &
!!                  NPoints,Npg, ncisp, ncipg, ncik, R1, R2, Rpg, Rk, &
!!                  Roota,Rootb,Rooty,Weightsp,Weightpg,fsp,fpg,fk,TrunODE,comm)
!!      Print *, "VAP Ene=", Ene
!!      Print *, "S^2(PCC)= ", real(S2,kind=pr)
!!      Print *, "Sz(PCC)= ", real(Sz,kind=pr)
!	  TrunODE = 3
!      Call PbarHbarEneLR(H1,H2,Evecs,Ene,T1,T2,Z1,Z2,NOcc,NAO,NSO,ENuc,QJ,&
!                         SP,NPoints, Npg, ncisp, ncipg, ncik, R1, R2, Rpg, Rk, &
!                         Roota,Rootb,Rooty,Weightsp,Weightpg,fsp,fpg,fk,TrunODE,comm)
!	  Call CCSD_T_Ene(Fock, H2, T1, T2, NOcc, NSO-NOcc, NSO, ET)
!	  Print *, "(T) Ene=", ET
!      Print *, "CCSD(T) Ene=", Ene


!==============================================!
!  Lastly, deallocate memory and exit safely.  !
!==============================================!
      DeAllocate(Fock,      &
!                 T1, T2,                   &
                 Stat=IAlloc)
      If(IAlloc /= 0) Stop "Could not deallocate in main"
!      Call ShutDownPUCCTools

      End Subroutine PGCC

