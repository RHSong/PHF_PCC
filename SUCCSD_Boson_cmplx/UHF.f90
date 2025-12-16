
   Module UHF

   Use Precision
   Use Constants
   Use DIIS

   Private
   Public  :: DrvUHF, UHFMix, MKFock_MO, MKFock_AO, SemiCanon

   Contains

      Subroutine DrvUHF(Olap,OneH,ERI,EvecsA,EvecsB,ENuc,NOccA,NOccB,NAO,UseDIIS)
      Implicit None
      Integer,           Intent(In)    :: NOccA, NOccB, NAO
      Logical,           Intent(In)    :: UseDIIS
      Real (Kind=pr),    Intent(In)    :: Olap(NAO,NAO), OneH(NAO,NAO)
      Real (Kind=pr),    Intent(In)    :: ERI(NAO,NAO,NAO,NAO)
      Real (Kind=pr),    Intent(In)    :: ENuc
      Real (Kind=pr),    Intent(InOut) :: EvecsA(NAO,NAO), EvecsB(NAO,NAO)
      Real (Kind=pr),    Parameter     :: TolMax1 = 1.0E-07_pr
      Real (Kind=pr),    Parameter     :: TolMax2 = 1.0E-08_pr
      Real (Kind=pr),    Allocatable   :: Trans(:,:)
      Real (Kind=pr),    Allocatable   :: Scr1(:,:), Scr2(:,:)
      Real (Kind=pr),    Allocatable   :: FockA(:,:), DenMatA(:,:), EvalsA(:)
      Real (Kind=pr),    Allocatable   :: FockB(:,:), DenMatB(:,:), EvalsB(:)
      Real (Kind=pr),    Allocatable   :: FockADIIS(:,:), FockAVec(:,:,:)
      Real (Kind=pr),    Allocatable   :: FockBDIIS(:,:), FockBVec(:,:,:)
      Integer           :: IAlloc(13)
      Real (Kind=pr)    :: RMS, ESCF
      Integer           :: NIter

!===================================================!
!  This routine does open-shell SCF calculations.   !
!  At the end of NewDen, DenHF is the new density.  !
!===================================================!

      Write(6,*) "Doing UHF..."
      Open(7,File='Output',Position='Append')
      Allocate(Trans(NAO,NAO),            Stat=IAlloc(1))
      Allocate(Scr1(NAO,NAO),             Stat=IAlloc(2))
      Allocate(Scr2(NAO,NAO),             Stat=IAlloc(3))
      Allocate(FockA(NAO,NAO),            Stat=IAlloc(4))
      Allocate(DenMatA(NAO,NAO),          Stat=IAlloc(5))
      Allocate(EvalsA(NAO),               Stat=IAlloc(6))
      Allocate(FockB(NAO,NAO),            Stat=IAlloc(7))
      Allocate(DenMatB(NAO,NAO),          Stat=IAlloc(8))
      Allocate(EvalsB(NAO),               Stat=IAlloc(9))
      Allocate(FockADIIS(NAO,NAO),        Stat=IAlloc(10))
      Allocate(FockBDIIS(NAO,NAO),        Stat=IAlloc(11))
      Allocate(FockAVec(NAO,NAO,NDIIS),   Stat=IAlloc(12))
      Allocate(FockBVec(NAO,NAO,NDIIS),   Stat=IAlloc(13))
      If(Any(IAlloc /= 0)) Stop 'Could not Allocate in UHF'


!=========================================================!
!  Get the transformation matrix and prepare for output.  !
!=========================================================!

      Write(7,*)
      Write(7,1000)
      Write(7,1010)
      Write(7,1020)
      Write(7,1030)
      Call GetTrans(Trans,Olap,NAO)


!===============!
!  Do the SCF.  !
!===============!

! Initial Guess
      NIter = 0
      DenMatA = Zero
      DenMatB = Zero
      Call NewDen(DenMatA,DenMatB,EvecsA,EvecsB,NOccA,NOccB,   &
                  RMS,NIter,NAO,Scr1,Scr2)
      Call MkFock(OneH,ERI,DenMatA,DenMatB,FockA,NAO)
      Call MkFock(OneH,ERI,DenMatB,DenMatA,FockB,NAO)
      Call SCFEnergy(OneH,FockA,FockB,DenMatA,DenMatB,ENuc,ESCF,NAO)
      Write(7,1040) ESCF,NIter,RMS

! Initialize DIIS
      FockADIIS = FockA
      FockBDIIS = FockB
      FockAVec  = Zero
      FockBVec  = Zero
      ErrVecs   = Zero

! Iterate to convergence
      Do
        Call FockDiag(FockADIIS,EvecsA,EvalsA,Trans,NAO,Scr1,Scr2)
        Call FockDiag(FockBDIIS,EvecsB,EvalsB,Trans,NAO,Scr1,Scr2)
        Call NewDen(DenMatA,DenMatB,EvecsA,EvecsB,NOccA,NOccB,   &
                    RMS,NIter,NAO,Scr1,Scr2)
        Call MkFock(OneH,ERI,DenMatA,DenMatB,FockA,NAO)
        Call MkFock(OneH,ERI,DenMatB,DenMatA,FockB,NAO)
        Call SCFEnergy(OneH,FockA,FockB,DenMatA,DenMatB,ENuc,ESCF,NAO)
        Write(7,1040) ESCF,NIter,RMS
        If(RMS < TolMax1) Exit
        Call DrvDIIS(Olap,Trans,FockA,FockB,FockADIIS,FockBDIIS,FockAVec,FockBVec,  &
                     DenMatA,DenMatB,Scr1,Scr2,RMS,NAO,NIter)
        If(.not. UseDIIS) Then
          FockADIIS = FockA
          FockBDIIS = FockB
        End If
      End Do

! TMH, experimental - turn off DIIS once we hit "convergence"
      Do
        Call FockDiag(FockA,EvecsA,EvalsA,Trans,NAO,Scr1,Scr2)
        Call FockDiag(FockB,EvecsB,EvalsB,Trans,NAO,Scr1,Scr2)
        Call NewDen(DenMatA,DenMatB,EvecsA,EvecsB,NOccA,NOccB,   &
                    RMS,NIter,NAO,Scr1,Scr2)
        Call MkFock(OneH,ERI,DenMatA,DenMatB,FockA,NAO)
        Call MkFock(OneH,ERI,DenMatB,DenMatA,FockB,NAO)
        Call SCFEnergy(OneH,FockA,FockB,DenMatA,DenMatB,ENuc,ESCF,NAO)
        Write(7,1040) ESCF,NIter,RMS
        If(RMS < TolMax2) Exit
      End Do


      Write(7,1020)
      Write(7,1050) NIter
      Write(7,1060) ESCF
      Write(7,1000)
      Close(7)


!===============================!
!  Deallocate and exit safely.  !
!===============================!

      DeAllocate(Trans,       Stat=IAlloc(1))
      DeAllocate(Scr1,        Stat=IAlloc(2))
      DeAllocate(Scr2,        Stat=IAlloc(3))
      DeAllocate(FockA,       Stat=IAlloc(4))
      DeAllocate(DenMatA,     Stat=IAlloc(5))
      DeAllocate(EvalsA,      Stat=IAlloc(6))
      DeAllocate(FockB,       Stat=IAlloc(7))
      DeAllocate(DenMatB,     Stat=IAlloc(8))
      DeAllocate(EvalsB,      Stat=IAlloc(9))
      DeAllocate(FockADIIS,   Stat=IAlloc(10))
      DeAllocate(FockBDIIS,   Stat=IAlloc(11))
      DeAllocate(FockAVec,    Stat=IAlloc(12))
      DeAllocate(FockBVec,    Stat=IAlloc(13))
      If(Any(IAlloc /= 0)) Stop 'Could not DeAllocate in UHF'

1000  Format(14x,'**************************************************')
1010  Format(14X,'*               UHF summary follows              *')
1020  Format(14x,'*------------------------------------------------*')
1030  Format(14X,'*   UHF Energy     Iteration    Biggest change   *')
1040  Format(14X,'* ',F15.10,2x,I5,7X,F14.10,4x,'*')
1050  Format(14x,'* SCF has converged in ',I4,' iterations',11x,'*')
1060  Format(14X,'* Final UHF Energy is ',2x,F18.12,' a.u.',2x,'*')

      Return
      End Subroutine DrvUHF






      Subroutine GetTrans(Trans,Olap,NAO)
      Implicit None
      Integer,        Intent(In)  :: NAO
      Real (Kind=pr), Intent(In)  :: Olap(NAO,NAO)
      Real (Kind=pr), Intent(Out) :: Trans(NAO,NAO)
      Real (Kind=pr), Allocatable :: UMat(:,:), SVals(:), Scr(:,:)
      Real (Kind=pr), Parameter   :: TolMax = 1.0E-9_pr
      Integer :: I
     
!=================================================!
!  This generates the canonically orthogonalized  !
!  transformation matrix, X = U*S^(-1/2).         !
!=================================================!

! Get U, s
      Allocate(UMat(NAO,NAO), SVals(NAO), Scr(NAO,NAO))
      Call DiagR(Olap,SVals,UMat,NAO)

! Form s^(-1/2)
      Scr = Zero
      Do I=1,NAO
       If(SVals(I) <= TolMax) Stop 'Linear dependence in basis set'
       Scr(I,I) = One/Sqrt(SVals(I))
      End Do

! Build X, which is block-diagonal
      Trans = MatMul(UMat,Scr)
      DeAllocate(UMat,SVals,Scr)

      Return
      End Subroutine GetTrans






      Subroutine NewDen(DenA,DenB,EvecsA,EvecsB,NOccA,NOccB,     &
                        RMS,NIter,NAO,DenHF2,DiffM)
      Implicit None
      Integer,        Intent(In)    :: NOccA, NoccB, NAO
      Integer,        Intent(InOut) :: NIter
      Real (Kind=pr), Intent(In)    :: EvecsA(NAO,NAO), EvecsB(NAO,NAO)
      Real (Kind=pr), Intent(InOut) :: DenA(NAO,NAO),   DenB(NAO,NAO)
      Real (Kind=pr), Intent(Out)   :: RMS
      Real (Kind=pr) :: DenHF2(NAO,NAO), DiffM(NAO,NAO)
      Integer, Parameter :: SCFMaxCyc = 1500

!=======================================================================!
!  This generates the density matrix from the SCF eigenvectors, forms   !
!  the convergence criterion (I'm using largest absolute change in the  !
!  density matrix elements) and updates the density matrix.             !
!=======================================================================!

! Do the alpha MOs
      DiffM  = Transpose(EvecsA)
      DenHF2 = MatMul(EvecsA(:,1:NOccA),DiffM(1:NOccA,:))
      DiffM  = DenHF2-DenA
      DenA   = DenHF2
      RMS    = MaxVal(Abs(DiffM))

! Do the beta MOs
      If(NOccB == 0) GoTo 10
      DiffM  = Transpose(EvecsB)
      DenHF2 = MatMul(EvecsB(:,1:NOccB),DiffM(1:NOccB,:))
      DiffM  = DenHF2-DenB
      DenB   = DenHF2
      RMS    = Max(MaxVal(Abs(DiffM)),RMS)

! Check convergence and out!
10    NIter = NIter + 1
      If(NIter >= SCFMaxCyc) Stop 'SCF did not converge'

      Return
      End Subroutine NewDen






      Subroutine MKFock(OneH,ERI,Den1,Den2,Fock,NAO)
      Implicit None
      Integer,           Intent(In)  :: NAO
      Real (Kind=pr),    Intent(In)  :: OneH(NAO,NAO)
      Real (Kind=pr), Intent(In)  :: ERI(NAO,NAO,NAO,NAO)
      Real (Kind=pr),    Intent(In)  :: Den1(NAO,NAO), Den2(NAO,NAO)
      Real (Kind=pr),    Intent(Out) :: Fock(NAO,NAO)
      Integer :: Mu, Nu, Lam, Kap

!=================================================!
!  This forms the Fock matrix, obviously enough.  !
!=================================================!

      Fock = OneH
      Do Mu = 1,NAO
       Do Nu = Mu,NAO
        Do Lam = 1,NAO
         Do Kap  = Lam+1,NAO
          Fock(Mu,Nu) = Fock(Mu,Nu)                                        &
                      + ERI(Mu,Nu,Lam,Kap)*Two*(Den1(Lam,Kap)              &
                                              + Den2(Lam,Kap))             &
                      - ERI(Mu,Lam,Nu,Kap)*Den1(Lam,Kap)                   &
                      - ERI(Mu,Kap,Nu,Lam)*Den1(Lam,Kap)
         End Do
         Fock(Mu,Nu) = Fock(Mu,Nu)                                         &
                     + ERI(Mu,Nu,Lam,Lam)*(Den1(Lam,Lam) + Den2(Lam,Lam))  &
                     - ERI(Mu,Lam,Nu,Lam)*Den1(Lam,Lam)
        End Do   
        Fock(Nu,Mu) = Fock(Mu,Nu)
       End Do
      End Do

      Return
      End Subroutine MKFock






      Subroutine SCFEnergy(OneH,FockA,FockB,DenA,DenB,ENuc,ESCF,NAO)
      Implicit None
      Integer,        Intent(In)  :: NAO
      Real (Kind=pr), Intent(In)  :: ENuc, OneH(NAO*NAO)
      Real (Kind=pr), Intent(In)  :: FockA(NAO*NAO), DenA(NAO*NAO)
      Real (Kind=pr), Intent(In)  :: FockB(NAO*NAO), DenB(NAO*NAO)
      Real (Kind=pr), Intent(Out) :: ESCF

!======================================!
!  This gets the total energy from     !
!  E = 1/2*P(i,j) * [F(i,j) + H(i,j)]  ! 
!======================================!

      ESCF = ENuc + Dot_Product(DenA,FockA+OneH)/2   &
                  + Dot_Product(DenB,FockB+OneH)/2

      Return
      End Subroutine SCFEnergy






      Subroutine FockDiag(Fock,Evec,Evals,Trans,NAO,FockMO,Scr)
      Implicit None
      Integer,        Intent(In)  :: NAO
      Real (Kind=pr), Intent(In)  :: Trans(NAO,NAO), Fock(NAO,NAO)
      Real (Kind=pr), Intent(Out) :: Evec(NAO,NAO), Evals(NAO)
      Real (Kind=pr) :: FockMO(NAO,NAO), Scr(NAO,NAO)

!==================================================================!
!  This generates the MO basis Fock matrix, FockMO                 !
!  We diagonalize it in EIG (leaving eigenvalues on the diagonal)  !
!  In this EIG call, the eigenvectors are put in Scr               !
!  The real eigenvectors are then made by another transformation   ! 
!==================================================================!

      FockMO = MatMul(Fock,Trans)
      Scr    = Transpose(Trans)
      FockMO = MatMul(Scr,FockMO)
      Call DiagR(FockMO,Evals,Scr,NAO)
      Evec = MatMul(Trans,Scr)

      Return
      End Subroutine FockDiag






      Subroutine DrvDIIS(Olap,Trans,FockA,FockB,FockADIIS,FockBDIIS,FockAVec,FockBVec,  &
                         DenMatA,DenMatB,Scr1,Scr2,RMS,NAO,NIter)
      Implicit None
      Integer,        Intent(In)  :: NAO, NIter
      Real (Kind=pr), Intent(In)  :: RMS
      Real (Kind=pr), Intent(In)  :: Olap(NAO,NAO),  Trans(NAO,NAO)
      Real (Kind=pr), Intent(In)  :: FockA(NAO,NAO), DenMatA(NAO,NAO)
      Real (Kind=pr), Intent(In)  :: FockB(NAO,NAO), DenMatB(NAO,NAO)
      Real (Kind=pr), Intent(Out) :: FockADIIS(NAO,NAO), FockBDIIS(NAO,NAO)
      Real (Kind=pr) :: FockAVec(NAO,NAO,NDIIS), FockBVec(NAO,NAO,NDIIS)
      Real (Kind=pr) :: Scr1(NAO,NAO), Scr2(NAO,NAO)
! Local variables
      Integer :: I
      Real (Kind=pr) :: Coeffs(NDIIS)

!====================================!
!  This subroutine drives the DIIS.  !
!====================================!

! Step 1: Cycle the list of error vectors and Fock matrices
      Do I = 1,NDIIS-1  
        ErrVecs(:,I)    = ErrVecs(:,I+1)
        FockAVec(:,:,I) = FockAVec(:,:,I+1)
        FockBVec(:,:,I) = FockBVec(:,:,I+1)
      End Do

! Step 2: Append to the list of Fock matrices and error vectors
      FockAVec(:,:,NDIIS) = FockA
      FockBVec(:,:,NDIIS) = FockB
      Call GetErrVec(ErrVecs(:,NDIIS),Olap,Trans,FockA,FockB,DenMatA,DenMatB,   &
                     Scr1,Scr2,NAO)

! Step 3: Construct FockDIIS, the Fock matrix to be diagonalized next
      If(NIter >= StartDIIS) Then
        Call DoDIIS(Coeffs)
        FockADIIS = Zero
        FockBDIIS = Zero
        Do I = 1,NDIIS
          FockADIIS = FockADIIS + FockAVec(:,:,I)*Coeffs(I)
          FockBDIIS = FockBDIIS + FockBVec(:,:,I)*Coeffs(I)
        End Do
      Else
        FockADIIS = FockA           ! By default, we diagonalize the new Fock matrix
        FockBDIIS = FockB           ! By default, we diagonalize the new Fock matrix
      End If

      Return
      End Subroutine DrvDIIS






      Subroutine GetErrVec(Vec,S,X,FA,FB,PA,PB,Scr1,Scr2,NAO)
      Implicit None
      Integer,        Intent(In)  :: NAO
      Real (Kind=pr), Intent(In)  :: S(NAO,NAO), X(NAO,NAO)
      Real (Kind=pr), Intent(In)  :: FA(NAO,NAO), FB(NAO,NAO)
      Real (Kind=pr), Intent(In)  :: PA(NAO,NAO), PB(NAO,NAO)
      Real (Kind=pr), Intent(Out) :: Vec(2*NAO*NAO)
      Real (Kind=pr) :: Scr1(NAO,NAO), Scr2(NAO,NAO)

!=======================================================!
!  The error vector is the real part of the commutator  !
!  done in an orthogonalized basis                      !
!=======================================================!

! The commutators are FPS - SPF = C.  Do alpha first.
      Scr1 = MatMul(FA,PA)
      Scr2 = MatMul(Scr1,S)             ! C = F.P.S
      Scr1 = MatMul(PA,FA)
      Scr2 = Scr2 - MatMul(S,Scr1)      ! C = F.P.S - S.P.F
      Scr1 = MatMul(Scr2,X)             ! C.X
      Scr2 = MatMul(Transpose(X),Scr1)  ! X+.C.X
      Vec(1:NAO*NAO) = ReShape(Scr2, (/ NAO*NAO /))

! Do beta.
      Scr1 = MatMul(FB,PB)
      Scr2 = MatMul(Scr1,S)             ! C = F.P.S
      Scr1 = MatMul(PB,FB)
      Scr2 = Scr2 - MatMul(S,Scr1)      ! C = F.P.S - S.P.F
      Scr1 = MatMul(Scr2,X)             ! C.X
      Scr2 = MatMul(Transpose(X),Scr1)  ! X+.C.X
      Vec(NAO*NAO+1:2*NAO*NAO) = ReShape(Scr2, (/NAO*NAO /))

      Return
      End Subroutine GetErrVec






      Subroutine UHFMix(EvecsA,EvecsB,NAO,NOccA,NOccB)
      Implicit None
      Integer, Intent(In) :: NoccA, NOccB, NAO
      Real (Kind=pr)      :: EvecsA(NAO,NAO), EvecsB(NAO,NAO)
      Real (Kind=pr), Allocatable :: V1(:), V2(:)

!==================================================!
!  Mixes alpha HOMO and LUMO, beta HOMO and LUMO.  !
!==================================================!

      Allocate(V1(NAO), V2(NAO))
! Do alpha
      V1 = EvecsA(:,NOccA)
      V2 = EvecsA(:,NOccA+1)
      EvecsA(:,NOccA)   = (V1 + V2)/Sqrt(Two)
      EvecsA(:,NOccA+1) = (V1 - V2)/Sqrt(Two)
! Do beta
      If(NOccB == 0) GoTo 10
      V1 = EvecsB(:,NOccB)
      V2 = EvecsB(:,NOccB+1)
      EvecsB(:,NOccB)   = (V1 - V2)/Sqrt(Two)
      EvecsB(:,NOccB+1) = (V1 + V2)/Sqrt(Two)
10    DeAllocate(V1, V2)

      Return
      End Subroutine UHFMix



      Subroutine MKFock_MO(OneH,ERI,Fock,ESCF,NOcc,NAO)
      Implicit None
      Integer,           Intent(In)  :: NAO, NOcc
      Complex (Kind=pr), Intent(In)  :: OneH(NAO,NAO)
      Complex (Kind=pr), Intent(In)  :: ERI(NAO,NAO,NAO,NAO)
      Complex (Kind=pr), Intent(Out) :: Fock(NAO,NAO)
      Complex (Kind=pr), Intent(Out) :: ESCF 
      Integer :: p, q, i, j 

!=================================================!
!  This forms the Fock matrix, obviously enough.  !
!  ERI should be put in Dirac notation            !
!=================================================!

      Fock = OneH
      !$omp parallel do reduction(+:Fock)
        Do i = 1, NOcc
        Fock = Fock + ERI(:,i,:,i) 
        EndDo
      !$omp end parallel do
      ESCF = Zero
      !$omp parallel do reduction(+:ESCF)
      Do i = 1, NOcc
        ESCF = ESCF + OneH(i,i)
      EndDo
      !$omp end parallel do

      !$omp parallel do reduction(+:ESCF)
      Do i = 1, NOcc
      Do j = 1, NOcc
        ESCF = ESCF + ERI(i,j,i,j)/Two 
      EndDo 
      EndDo 
      !$omp end parallel do

      Return
      End Subroutine MKFock_MO




      Subroutine MKFock_AO(OneH,ERI,Fock,Ref,NOcc,NSO)
      Implicit None
      Integer,           Intent(In)  :: NSO, NOcc
      Complex (Kind=pr), Intent(In)  :: OneH(NSO,NSO)
      Complex (Kind=pr), Intent(In)  :: Ref(NSO,NSO)
      Complex (Kind=pr), Intent(In)  :: ERI(NSO,NSO,NSO,NSO)
      Complex (Kind=pr), Intent(Out) :: Fock(NSO,NSO)
      Complex (Kind=pr), Allocatable :: rho(:,:)
      Integer :: p, q, r, s, i, j 

!=================================================!
!  This forms the Fock matrix, obviously enough.  !
!  ERI should be put in Dirac notation            !
!=================================================!

      Allocate(rho(NSO,NSO))

      rho = Zero
      Do p = 1, NSO
      Do q = 1, NSO
      Do i = 1, NOcc
        rho(p,q) = rho(p,q) + ref(p,i)*ref(q,i)
      EndDo
      EndDo
      EndDo

      Fock = OneH
      Do p = 1, NSO
      Do q = 1, NSO
      Do r = 1, NSO 
      Do s = 1, NSO 
        Fock(p,q) = Fock(p,q) + ERI(p,r,q,s) * rho(r,s) 
      EndDo
      EndDo
      EndDo
      EndDo

      Deallocate(rho)

      Return
      End Subroutine MKFock_AO

	  Subroutine SemiCanon(OneH, ERI, Ref, NOcc, NSO, cMO, Roo, Rvv)
	  Implicit None
	  Integer,           Intent(In)  :: NSO, NOcc
	  Complex (Kind=pr), Intent(In)  :: OneH(NSO,NSO)
	  Complex (Kind=pr), Intent(In)  :: Ref(NSO,NSO)
	  Complex (Kind=pr), Intent(In)  :: ERI(NSO,NSO,NSO,NSO)
	  Complex (Kind=pr), Intent(Out) :: cMO(NSO,NSO)
	  Complex (Kind=pr), Intent(Out) :: Roo(NOcc,NOcc), Rvv(NSO-NOcc,NSO-NOcc)
	  Complex (Kind=pr) :: F(NSO,NSO)
	  Complex (Kind=pr) :: Foo(NOcc,NOcc), Fvv(NSO-NOcc,NSO-NOcc)
	  Complex (Kind=pr) :: Eo(NOcc), Ev(NSO-NOcc), Fock(NSO,NSO)
	  Integer :: I
	  Call MKFock_AO(OneH, ERI, Fock, Ref, NOcc, NSO)
	  Foo = Fock(1:NOcc,1:NOcc)
	  Fvv = Fock(NOcc+1:NSO,NOcc+1:NSO)
	  Call DiagC2(Foo, Eo, Roo, NOcc)
	  cMO(:,1:NOcc) = MatMul(Ref(:,1:NOcc), Roo)
	  Call DiagC2(Fvv, Ev, Rvv, NSO-NOcc)
	  cMO(:,NOcc+1:NSO) = MatMul(Ref(:,NOcc+1:NSO), Rvv)
	  Roo = Zero
	  Do I = 1, NOcc
	  	Roo(I,I) = Eo(I)
	  End Do
	  Rvv = Zero
	  Do I = 1, NSO-NOcc
	  	Rvv(I,I) = Ev(I)
	  End Do
	  Print *, "Test Foo", maxval(abs(Foo-Roo))
	  Print *, "Test Fvv", maxval(abs(Fvv-Rvv))
	  stop
	  End Subroutine SemiCanon

   End Module UHF

