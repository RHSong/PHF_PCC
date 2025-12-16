
      Subroutine DiagR(Matrix,Evals,Evecs,N)
      Use Precision
      Implicit None
      Integer,        Intent(In)  :: N
      Real (Kind=pr), Intent(In)  :: Matrix(N,N)
      Real (Kind=pr), Intent(Out) :: Evecs(N,N), Evals(N)
      Real (Kind=pr), Allocatable :: Work(:)
      Integer :: NBlocks, LWork, Info, ILAENV

!==========================================!
!  Non-destructive wrapper to DSyEv.       !
!  Eigenvalues sorted in ascending order.  !
!==========================================!

! Find LWork and Allocate
      NBlocks = ILAENV(1,"DSYTRD","U",N,-1,-1,-1) + 2
      LWork   = N*NBlocks
      Allocate(Work(LWork),  Stat=Info)
      If(Info /= 0) Stop "Could not allocate in DiagR"

! Diagonalize!
      Evecs = Matrix
      Call DSyEv("V","U",N,Evecs,N,Evals,Work,LWork,Info)
      If(Info /= 0) Stop "Error in LAPack!"

! Deallocate and exit safely
      DeAllocate(Work,  Stat=Info)
      If(Info /= 0) Stop "Could not deallocate in DiagR"

      Return
      End Subroutine DiagR






      Subroutine DiagC(Matrix,Evals,Evecs,N)
      Use Precision
      Implicit None
      Integer,           Intent(In)  :: N
      Complex (Kind=pr), Intent(In)  :: Matrix(N,N)
      Complex (Kind=pr), Intent(Out) :: Evecs(N,N)
      Real (Kind=pr),    Intent(Out) :: Evals(N)
      Complex (Kind=pr), Allocatable :: CWork(:)
      Real (Kind=pr),    Allocatable :: RWork(:)
      Integer :: NBlocks, LWork, Info, ILAENV

!==========================================!
!  Non-destructive wrapper to ZHeEv.       !
!  Eigenvalues sorted in ascending order.  !
!==========================================!

! Find LWork and Allocate
      NBlocks = ILAENV(1,"ZHETRD","U",N,-1,-1,-1) + 1
      LWork   = N*NBlocks
      Allocate(CWork(LWork), RWork(3*N-2), Stat=Info)
      If(Info /= 0) Stop "Could not allocate in DiagC"

! Diagonalize!
      Evecs = Matrix
      Call ZHeEv("V","U",N,Evecs,N,Evals,CWork,LWork,RWork,Info)
      If(Info /= 0) Stop "Error in LAPack!"

! Deallocate and exit safely
      DeAllocate(CWork, RWork,  Stat=Info)
      If(Info /= 0) Stop "Could not deallocate in DiagC"

      Return
      End Subroutine DiagC



      Subroutine DiagR2(Matrix,Evals,NDim)
      Use Precision
      Implicit None
      Integer, Intent(In)            :: NDim
      Real (Kind=pr),    Intent(In)  :: Matrix(NDim,NDim)
      Complex (Kind=pr), Intent(Out) :: Evals(NDim)
      Real (Kind=pr),    Allocatable :: EvalsR(:), EvalsI(:)
      Real (Kind=pr),    Allocatable :: LEvecs(:,:), REvecs(:,:)
      Real (Kind=pr),    Allocatable :: Scr(:,:), Work(:)
      Integer :: LenScr, Info

!===============================================!
!  This is a non-destructive wrapper to DGeEv.  !
!===============================================!

! Allocate workspace
      LenScr = 3*NDim
      Allocate(EvalsR(NDim), EvalsI(NDim), Work(LenScr),      &
               Scr(NDim,NDim), LEvecs(1,1), REvecs(1,1),      &
               Stat=Info)
      If(Info /= 0) Stop "Could not allocate in DiagR2!"


! Do the actual diagonalization
      Scr = Matrix
      Call DGeEv("N","N",NDim,Scr,NDim,EvalsR,EvalsI,      &
                 LEvecs,NDim,REvecs,NDim,Work,LenScr,Info)
      If(Info /= 0) Stop "Error in DGeEv!"

! Complexify the eigenvalues
      Evals = Cmplx(EvalsR,EvalsI,Kind=pr)

! Deallocate and exit safely
      DeAllocate(EvalsR, EvalsI, Work, Scr, LEvecs, REvecs, Stat=Info)
      If(Info /= 0) Stop "Could not deallocate in DiagR2!"

      Return
      End Subroutine DiagR2


      Subroutine DiagC2(Matrix,Evals,NDim)
      Use Precision
      Implicit None
      Integer, Intent(In)            :: NDim
      Complex (Kind=pr), Intent(In)  :: Matrix(NDim,NDim)
      Complex (Kind=pr), Intent(Out) :: Evals(NDim)
      Complex (Kind=pr), Allocatable :: LEvecs(:,:), REvecs(:,:)
      Complex (Kind=pr), Allocatable :: Scr(:,:), Work(:)
      Real (Kind=pr),    Allocatable :: RWork(:)
      Integer :: LenScr, Info

!===============================================!
!  This is a non-destructive wrapper to DGeEv.  !
!===============================================!

! Allocate workspace
      LenScr = 3*NDim
      Allocate(Work(LenScr), RWork(2*NDim),                 &
               Scr(NDim,NDim), LEvecs(1,1), REvecs(1,1),    &
               Stat=Info)
      If(Info /= 0) Stop "Could not allocate in DiagC2!"


! Do the actual diagonalization
      Scr = Matrix
      Call ZGeEv("N","N",NDim,Scr,NDim,Evals,      &
                 LEvecs,NDim,REvecs,NDim,Work,LenScr,RWork,Info)
      If(Info /= 0) Stop "Error in ZGeEv!"

! Deallocate and exit safely
      DeAllocate(Work, RWork, Scr, LEvecs, REvecs, Stat=Info)
      If(Info /= 0) Stop "Could not deallocate in DiagC2!"

      Return
      End Subroutine DiagC2



      Subroutine DeterminantR(det,Mat,N)
      Use Precision
      Use Constants
      Implicit None
      Integer,        Intent(In)  :: N
      Real (Kind=pr), Intent(In)  :: Mat(N,N)
      Real (Kind=pr), Intent(Out) :: det
      Complex (Kind=pr) :: Evals(N), detC
      Integer :: I

      Call DiagR2(Mat,Evals,N)
      detC = One
      Do I = 1,N
        detC = detC*Evals(I)
      End Do
      det = Real(detC,Kind=pr)


      Return
      End Subroutine DeterminantR


!     Subroutine DeterminantC(det,Mat,N)
!     Use Precision
!     Use Constants
!     Implicit None
!     Integer,           Intent(In)  :: N
!     Complex (Kind=pr), Intent(In)  :: Mat(N,N)
!     Complex (Kind=pr), Intent(Out) :: det
!     Complex (Kind=pr) :: Evals(N)
!     Integer :: I
!
!     Call DiagC2(Mat,Evals,N)
!     det = One
!     Do I = 1,N
!       det = det*Evals(I)
!     End Do
!
!     Return
!     End Subroutine DeterminantC


      Subroutine determinantC(det,Mat,N)
! Use LU decomposition to calculate determinant
      Use Precision
      Use Constants
      Implicit None
      Integer,           Intent(In)  :: N
      Complex (Kind=pr), Intent(In)  :: Mat(N,N)
      Complex (Kind=pr), Intent(Out) :: det
      Complex (Kind=pr), Allocatable :: Scr(:,:) 
      Integer,           Allocatable :: Piv(:) 
      Integer :: I, Info
      Allocate(Scr(N,N), Piv(N), Stat=Info)
      If(Info /= 0) Stop "Could not allocate in determinant"
! Perform LU decomposition with pivot
      Scr = Mat
      Call ZGETRF(N,N,Scr,N,PIV,Info)
!     Write(*,*) "Pivot", PIV
      If (Info>0) then
! Some U(i,i) is exactly zero
        det = Zero
      Else If (Info==0) then
        det = One
        Do I = 1, N
! Consider the sign introduced by pivoting
          If (PIV(I)==I) then
            det = det*Scr(I,I)
          Else
            det = -det*Scr(I,I)
          EndIf
        EndDo
      Else
        Stop 'Error in LU decomposition of SR determinant'
      EndIf
      Deallocate(Scr, Piv, Stat=Info)
      If(Info /= 0) Stop "Could not deallocate in determinant"
      Return
      End Subroutine determinantC


      Subroutine traceC(tr,Mat,N)
! Calculate trace
      Use Precision
      Use Constants
      Implicit None
      Integer,           Intent(In)  :: N
      Complex (Kind=pr), Intent(In)  :: Mat(N,N)
      Complex (Kind=pr), Intent(Out) :: tr 
      Integer :: I
      tr = Zero
      Do I = 1, N
        tr = tr + Mat(I,I)
      EndDo
      Return
      End Subroutine traceC


      Subroutine InvertR(Matrix,IMatrix,N)
      Use Precision
      Implicit None
      Integer,        Intent(In)  :: N
      Real (Kind=pr), Intent(In)  :: Matrix(N,N)
      Real (Kind=pr), Intent(Out) :: IMatrix(N,N)
      Real (Kind=pr), Allocatable :: Work(:)
      Integer :: Info, IPiv(N), LWork

!========================================!
!  Inverts a matrix, non-destructively.  !
!========================================!

! Copy the matrix and factorize it in whatever dgetrf does
      IMatrix = Matrix
      Call DGetRF(N,N,IMatrix,N,IPiv,Info)
      If(Info /= 0) Stop "Error in DGetRF"

! Do a work-space query
      LWork = N
      Allocate(Work(N), Stat=Info)
      If(Info /= 0) Stop "Could not allocate in Invert"
      Call DGetRI(N,IMatrix,N,IPiv,Work,-1,Info)
      LWork = NInt(Work(1))

! Deallocate and reallocate
      DeAllocate(Work, Stat=Info)
      If(Info /= 0) Stop "Could not deallocate in Invert"
      Allocate(Work(LWork), Stat=Info)
      If(Info /= 0) Stop "Could not allocate in Invert"

! Invert!
      Call DGetRI(N,IMatrix,N,IPiv,Work,LWork,Info)
      If(Info /= 0) Stop "Error in DGetRI"

! Deallocate and exit safely
      DeAllocate(Work)
      If(Info /= 0) Stop "Could not deallocate in Invert"

      Return
      End Subroutine InvertR



      Subroutine InvertC(Matrix,IMatrix,N)
      Use Precision
      Implicit None
      Integer,           Intent(In)  :: N
      Complex (Kind=pr), Intent(In)  :: Matrix(N,N)
      Complex (Kind=pr), Intent(Out) :: IMatrix(N,N)
      Complex (Kind=pr), Allocatable :: Work(:)
      Integer :: Info, IPiv(N), LWork

!========================================!
!  Inverts a matrix, non-destructively.  !
!========================================!

! Copy the matrix and factorize it in whatever dgetrf does
      IMatrix = Matrix
      Call ZGetRF(N,N,IMatrix,N,IPiv,Info)
      If(Info /= 0) Stop "Error in ZGetRF"

! Do a work-space query
      LWork = N
      Allocate(Work(N), Stat=Info)
      If(Info /= 0) Stop "Could not allocate in InvertC"
      Call ZGetRI(N,IMatrix,N,IPiv,Work,-1,Info)
      LWork = NInt(Real(Work(1)))

! Deallocate and reallocate
      DeAllocate(Work, Stat=Info)
      If(Info /= 0) Stop "Could not deallocate in InvertC"
      Allocate(Work(LWork), Stat=Info)
      If(Info /= 0) Stop "Could not allocate in InvertC"

! Invert!
      Call ZGetRI(N,IMatrix,N,IPiv,Work,LWork,Info)
      If(Info /= 0) Stop "Error in ZGetRI"

! Deallocate and exit safely
      DeAllocate(Work)
      If(Info /= 0) Stop "Could not deallocate in InvertC"

      Return
      End Subroutine InvertC


      Subroutine QRdecompR(Matrix,QMat,RMat,N)
      Use Precision
      Use Constants 
      Implicit None
      Integer,        Intent(In)  :: N
      Real (Kind=pr), Intent(In)  :: Matrix(N,N)
      Real (Kind=pr), Intent(Out) :: QMat(N,N), RMat(N,N)
      Real (Kind=pr), Allocatable :: Work(:), Scr(:,:), tau(:), v(:)
      Integer :: LWork, Info
      Integer :: I, J 

! Find LWork and Allocate
      LWork   = N
      Allocate(Work(LWork), Scr(N,N), tau(N), v(N), Stat=Info)
      If(Info /= 0) Stop "Could not allocate in qr"

! Diagonalize!
      Scr = Matrix
      Call DGEQR2(N,N,Scr,N,TAU,WORK,INFO)
      If(Info /= 0) Then
        Print *, "Error in LAPack: ", Info
        Stop
      End If
! Construct Q and R
      QMat = Zero
      Do I = 1, N
        QMat(I,I) = One
      EndDo
      Do I = 1, N
        RMat = Zero
        Do J = 1, N
          RMat(J,J) = One
        EndDo
        v        = Zero
        v(I)     = One
        v(I+1:N) = Scr(I+1:N,I)
        RMat = RMat - tau(I) *   &
          spread(v(1:N),dim=2,ncopies=N)*spread(v(1:n),dim=1,ncopies=n)
!       Write(*,*) "I, tau(I)",I, tau(I)
!       Call Outmat(6,N,1,1,v,"v")
!       Call Outmat(6,N,N,1,RMat,"1-tau*v*v'")
        QMat = Matmul(QMat,RMat)
      EndDo
      RMat = Zero
      Do I = 1, N
      Do J = I, N
        RMat(I,J) = Scr(I,J)
      EndDo
      EndDo

! Deallocate and exit safely
      DeAllocate(Work, Scr, tau, v,  Stat=Info)
      If(Info /= 0) Stop "Could not deallocate in qr"

      Return
      End Subroutine QRdecompR



      Subroutine QRdecompC(Matrix,QMat,RMat,N)
      Use Precision
      Use Constants 
      Implicit None
      Integer,           Intent(In)  :: N
      Complex (Kind=pr), Intent(In)  :: Matrix(N,N)
      Complex (Kind=pr), Intent(Out) :: QMat(N,N), RMat(N,N)
      Complex (Kind=pr), Allocatable :: Work(:), Scr(:,:), tau(:), v(:)
      Integer :: LWork, Info
      Integer :: I, J 

! Find LWork and Allocate
      LWork   = N
      Allocate(Work(LWork), Scr(N,N), tau(N), v(N), Stat=Info)
      If(Info /= 0) Stop "Could not allocate in qr"

! Diagonalize!
      Scr = Matrix
      Call ZGEQR2(N,N,Scr,N,TAU,WORK,INFO)
      If(Info /= 0) Then
        Print *, "Error in LAPack: ", Info
        Stop
      End If
! Construct Q and R
      QMat = Zero
      Do I = 1, N
        QMat(I,I) = One
      EndDo
      Do I = 1, N
        RMat = Zero
        Do J = 1, N
          RMat(J,J) = One
        EndDo
        v        = Zero
        v(I)     = One
        v(I+1:N) = Scr(I+1:N,I)
        RMat = RMat - tau(I) *   &
          spread(v(1:N),dim=2,ncopies=N)    &
        * spread(Conjg(v(1:N)),dim=1,ncopies=N)
!       Write(*,*) "I, tau(I)",I, tau(I)
!       Call Outmat(6,N,1,1,v,"v")
!       Call Outmat(6,N,N,1,RMat,"1-tau*v*v'")
        QMat = Matmul(QMat,RMat)
      EndDo
      RMat = Zero
      Do I = 1, N
      Do J = I, N
        RMat(I,J) = Scr(I,J)
      EndDo
      EndDo

! Deallocate and exit safely
      DeAllocate(Work, Scr, tau, v,  Stat=Info)
      If(Info /= 0) Stop "Could not deallocate in qr"

      Return
      End Subroutine QRdecompC

      Subroutine SetUpRandom()
      Implicit None
      Integer :: I, N, Clock
      Integer, Allocatable :: Seed(:)
! Sets up for random number stuff
      Call System_Clock(Count=Clock)
      Call Random_Seed(Size=N)
      Allocate(Seed(N))
! Seed random numbers for UHFMix
      Call System_Clock(Count=Clock)
      Seed = Clock + 37*(/(i-1,i=1,n)/)
      Call Random_Seed(Put=Seed)
      Deallocate(Seed) 
      End Subroutine SetUpRandom

      Subroutine SvdR(Matrix,U,S,V,M,N)
      Use Constants
      Use Precision
      Implicit None
      Integer,        Intent(In)  :: M, N
      Real (Kind=pr), Intent(In)  :: Matrix(M,N)
      Real (Kind=pr), Intent(Out) :: U(M,M), V(N,N), S(M)
      Real (Kind=pr), Allocatable :: Work(:), SVal(:), A(:,:)
      Integer :: LWork, Info, dimS, I

!==========================================!
!  Non-destructive wrapper to DGeSVD.      !
!  S values sorted in descending order.    !
!==========================================!
      dimS = Min(M,N)
      Allocate(SVal(dimS),      &
               A(M,N),          &
               Work(M),         &
               Stat=Info)
      If(Info /= 0) Stop "Could not allocate in svdR"

! Find LWork and Allocate
      A = Matrix
      Call DGeSVD("A","A",M,N,A,M,SVal,U,M,V,N,Work,-1,Info)
      If(Info /= 0) Stop "Could not get LWork in svdR"
      LWork   = Work(1)
      Deallocate(Work,       Stat=Info)
      Allocate(Work(LWork),  Stat=Info)
      If(Info /= 0) Stop "Could not allocate in svdR"

! Perform SVD 
      A = Matrix
      Call DGeSVD("A","A",M,N,A,M,SVal,U,M,V,N,Work,LWork,Info)
      If(Info /= 0) Stop "Error in DGeSVD!"
! Process S and V
      V = Transpose(V)
      S = Zero
      Do I = 1, dimS
        S(I) = SVal(I)
      EndDo

! Deallocate and exit safely
      DeAllocate(A, SVal, Work,  Stat=Info)
      If(Info /= 0) Stop "Could not deallocate in svdR"

      Return
      End Subroutine SvdR

      Subroutine SvdC(Matrix,U,S,V,M,N)
      Use Constants
      Use Precision
      Implicit None
      Integer,           Intent(In)  :: M, N
      Complex (Kind=pr), Intent(In)  :: Matrix(M,N)
      Complex (Kind=pr), Intent(Out) :: U(M,M), V(N,N)
      Real (Kind=pr),    Intent(Out) :: S(M)
      Complex (Kind=pr), Allocatable :: Work(:), A(:,:)
      Real (Kind=pr),    Allocatable :: RWork(:), SVal(:)
      Integer :: LWork, Info, dimS, I

!==========================================!
!  Non-destructive wrapper to ZGeSVD.      !
!  S values sorted in descending order.    !
!==========================================!
      dimS = Min(M,N)
      Allocate(SVal(dimS),      &
               A(M,N),          &
               Work(M),         &
               RWork(5*dimS),   &
               Stat=Info)
      If(Info /= 0) Stop "Could not allocate in svdC"

! Find LWork and Allocate
      A = Matrix
      Call ZGeSVD("A","A",M,N,A,M,SVal,U,M,V,N,Work,-1,RWork,Info)
      If(Info /= 0) Stop "Could not get LWork in svdC"
      LWork   = Work(1)
      Deallocate(Work,       Stat=Info)
      Allocate(Work(LWork),  Stat=Info)
      If(Info /= 0) Stop "Could not allocate in svdC"

! Perform SVD 
      A = Matrix
      Call ZGeSVD("A","A",M,N,A,M,SVal,U,M,V,N,Work,LWork,RWork,Info)
      If(Info /= 0) Stop "Error in ZGeSVD!"
! Process S and V
      V = Conjg(Transpose(V))
      S = Zero
      Do I = 1, dimS
        S(I) = SVal(I)
      EndDo

! Deallocate and exit safely
      DeAllocate(A, SVal, Work, RWork, Stat=Info)
      If(Info /= 0) Stop "Could not deallocate in svdC"

      Return
      End Subroutine SvdC

      Subroutine LinearSolveC(A,x,b,N)
! Solve Ax = b
      Use Precision
      Use Constants
      Implicit None
      Integer,           Intent(In)  :: N
      Complex (Kind=pr), Intent(In)  :: A(N,N), b(N)
      Complex (Kind=pr), Intent(Out) :: x(N)
      Complex (Kind=pr) :: U(N,N), V(N,N), Sinv(N,N)
      Real (Kind=pr) :: tol = 1e-8, S(N)
      Integer :: I
      Call SvdC(A,U,S,V,N,N)
      Sinv = Zero
      Do I = 1, N
         If (abs(real(S(i),kind=pr)) > tol) then
            Sinv(i,i) = One / S(i)
         Else
            Sinv(i,i) = Zero
         End If
      End Do
      V = MatMul(V,Sinv)
      V = MatMul(V,transpose(U))
      x = MatMul(V,b)
      Return
      End Subroutine LinearSolveC

      Subroutine EigenSolve(Matrix,Metric,N,EVecs,Evals)
      Use Precision
      Use Constants
! Solve HC = SCE, H and S should be real
      Implicit None
      Integer,           Intent(In)  :: N
      Real    (Kind=pr), Intent(In)  :: Matrix(N,N), Metric(N,N)
      Real    (Kind=pr), Intent(Out) :: Evals(N), Evecs(N,N)
      Real (Kind=pr) :: tol = 1e-8
      Real (Kind=pr), allocatable, Dimension(:,:) :: HTmp, Trans, VTmp
      Integer :: I, M
      Call DiagR(-Metric,EVals,EVecs,N)
      Evals = -Evals
      Do I = 1,N
        If(Evals(I) < Tol) Exit
      End Do
      M = I - 1
      Allocate(Trans(N,M), HTmp(M,M), VTmp(M,M))
      Do I = 1,M
        Trans(:,I) = Evecs(:,I)/Sqrt(Evals(I))
      End Do

      Evals = 0
      HTmp = MatMul(Transpose(Trans),MatMul(Matrix,Trans))
      Call DiagR(HTmp,Evals(1:M),VTmp,M)
      EVecs(1:N,1:M) = MatMul(Trans,VTmp)
      DeAllocate(Trans, HTmp, VTmp)
      End Subroutine EigenSolve

      Subroutine OuterProd(A,B,Prod,M,N)
      Use Constants
      Use Precision
      Implicit None
      Integer,           Intent(In)  :: M, N
      Complex (Kind=pr), Intent(In)  :: A(M), B(N)
      Complex (Kind=pr), Intent(Out) :: Prod(M,N)
      Integer :: I
      Do I = 1, M
        Prod(I,:) = A(I) * B(:)
      End Do
      Return
      End Subroutine OuterProd

	  Subroutine A_outer_B(A,B,C,MA,MB)
	  Use Constants
	  Use Precision
	  Implicit None
	  Integer, Intent(In)           ::  MA, MB
	  Complex (Kind=pr), Intent(In) ::  A(MA), B(MB)
	  Complex (Kind=pr), Intent(InOut) ::  C(MA,MB)
	  Call ZGERU(MA,MB,One,A,1,B,1,C,MA)
	  End Subroutine A_outer_B

      Subroutine A_dot_B(A,B,C,MA,NA,NB)
      Use Constants
      Use Precision
      Implicit None
      Integer, Intent(In)           ::  MA,NA,NB
      Complex (Kind=pr), Intent(In)    ::  A(MA,NA), B(NA,NB)
      Complex (Kind=pr), Intent(InOut) ::  C(MA,NB)
      integer, allocatable :: workload(:), istart(:), iend(:)
      integer :: m, nblk, blksize = 8
      Call ZGEMM('N','N',MA,NB,NA,ZOne,A,MA,B,NA,ZOne,C,MA)
      End Subroutine A_dot_B

      Subroutine A_dot_Bh(A,B,C,MA,NA,MB)
      Use Constants
      Use Precision
      Implicit None
      Integer, Intent(In)           ::  MA,NA,MB
      Complex (Kind=pr), Intent(In)    :: A(MA,NA), B(MB,NA)
      Complex (Kind=pr), Intent(InOut) :: C(MA,MB)
      Call ZGEMM('N','C',MA,MB,NA,ZOne,A,MA,B,MB,ZOne,C,MA)
      End Subroutine A_dot_Bh

      Subroutine A_dot_Bt(A,B,C,MA,NA,MB)
      Use Constants
      Use Precision
      Implicit None
      Integer, Intent(In)           ::  MA,NA,MB
      Complex (Kind=pr), Intent(In)    :: A(MA,NA), B(MB,NA)
      Complex (Kind=pr), Intent(InOut) :: C(MA,MB)
      Call ZGEMM('N','T',MA,MB,NA,ZOne,A,MA,B,MB,ZOne,C,MA)
      End Subroutine A_dot_Bt

      Subroutine Ah_dot_B(A,B,C,NA,MA,NB)
      Use Constants
      Use Precision
      Integer, Intent(In)           ::  MA,NA,NB
      Complex (Kind=pr), Intent(In)    :: A(NA,MA), B(NA,NB)
      Complex (Kind=pr), Intent(InOut) ::  C(MA,NB)
      Call ZGEMM('C','N',MA,NB,NA,ZOne,A,NA,B,NA,ZOne,C,MA)
      End Subroutine Ah_dot_B

      Subroutine At_dot_B(A,B,C,NA,MA,NB)
      Use Constants
      Use Precision
      Integer, Intent(In)           ::  MA,NA,NB
      Complex (Kind=pr), Intent(In)    :: A(NA,MA), B(NA,NB)
      Complex (Kind=pr), Intent(InOut) ::  C(MA,NB)
      Call ZGEMM('T','N',MA,NB,NA,ZOne,A,NA,B,NA,ZOne,C,MA)
      End Subroutine At_dot_B

      Subroutine Ah_dot_Bh(A,B,C,NA,MA,MB)
      Use Constants
      Use Precision
      Integer, Intent(In)           ::  MA,NA,MB
      Complex (Kind=pr), Intent(In)    :: A(NA,MA), B(MB,NA)
      Complex (Kind=pr), Intent(InOut) :: C(MA,MB)
      Call ZGEMM('C','C',MA,MB,NA,ZOne,A,NA,B,MB,ZOne,C,MA)
      End Subroutine Ah_dot_Bh

      Subroutine At_dot_Bt(A,B,C,NA,MA,MB)
      Use Constants
      Use Precision
      Integer, Intent(In)           ::  MA,NA,MB
      Complex (Kind=pr), Intent(In)    :: A(NA,MA), B(MB,NA)
      Complex (Kind=pr), Intent(InOut) :: C(MA,MB)
      Call ZGEMM('T','T',MA,MB,NA,ZOne,A,NA,B,MB,ZOne,C,MA)
      End Subroutine At_dot_Bt

      Subroutine mA_dot_mB(A,B,C,MA,NA)
      Use Constants
      Use Precision
      ! calculate C = Aij Bij
      Complex (Kind=pr), Intent(In)    :: A(MA,NA), B(MA,NA)
      Complex (Kind=pr), Intent(Out)   :: C(1,1)
      Integer :: N
      N = MA * NA
      C = Zero
      Call A_dot_B(A,B,C,1,N,1)
      End Subroutine mA_dot_mB

      Subroutine mA_dot_mBt(A,B,C,MA,NA)
      Use Constants
      Use Precision
      ! calculate C = Aij Bji
      Complex (Kind=pr), Intent(In)    :: A(MA,NA), B(NA,MA)
      Complex (Kind=pr), Intent(Out)   :: C(1,1)
      Integer :: N
      N = MA * NA
      C = Zero
      Call A_dot_B(A,transpose(B),C,1,N,1)
      End Subroutine mA_dot_mBt

      Subroutine Transpose1(A,M)
      Use Constants
      Use Precision
      Integer, Intent(In)            :: M
      Complex (Kind=pr), Intent(InOut) :: A(M,M)
      A = Transpose(A)
      End Subroutine Transpose1

      Subroutine Transpose2(A,B,M,N)
      Use Constants
      Use Precision
      Integer, Intent(In)            :: M,N
      Complex (Kind=pr), Intent(In)  :: A(M,N)
      Complex (Kind=pr), Intent(Out) :: B(N,M)
      B = Transpose(A)
      End Subroutine Transpose2

      Subroutine IDMat(M,N)
      Use Constants
      Use Precision
      Integer, Intent(In)            :: N
      Complex (Kind=pr), Intent(Out) :: M(N,N)
      Integer :: I
      M = Zero
      Do I = 1, N
        M(I,I) = ZOne
      End Do
      End Subroutine IDMat

      Subroutine Wigner(J,N,M,alpha,beta,gama,wig)
      Use Constants
      Use Precision
      Implicit None
      Integer         , Intent(In)  :: J,N,M
      Real(kind=pr)   , Intent(In)  :: alpha,beta,gama
      Complex(kind=pr), Intent(Out) :: wig
      Real(kind=pr) :: Djnm, numer, denom, fac
      Integer :: s, smin, smax, p1, p2, p3
      smax = min(j+m,j-n)
      smin = max(0,m-n)
      fac = One
      fac = fac * gamma(j+n+One) * gamma(j-n+One)
      fac = fac * gamma(j+m+One) * gamma(j-m+One)
      fac = sqrt(fac)
      Djnm = 0
      Do s = smin, smax
        p1 = n - m + s
        p2 = 2*j + m - n -2*s
        p3 = n - m + 2*s
        numer = (-1)**p1 * (cos(beta/2))**p2 * (sin(beta/2))**p3
        denom = One
        denom = denom * gamma(j+m-s+One) * gamma(s+One)
        denom = denom * gamma(n-m+s+One) * gamma(j-n-s+One)
        Djnm = Djnm + numer / denom
      End Do
      Djnm = Djnm * fac
      wig = exp(n*Cmplx(Zero,-alpha,kind=pr)) * Djnm * exp(m*Cmplx(Zero,-gama,kind=pr))
      End Subroutine Wigner

      Function Delta(i,j)
      Implicit None
      Integer, intent(in) :: i,j
      Integer :: Delta
      If (i==j) Then
        Delta = 1
      Else
        Delta = 0
      End If
      End Function

      Function Factorial(n)
      Implicit None
      Integer, intent(in) :: n
      Integer :: Factorial
      Integer :: I
      Factorial = 1
      Do I = 2, n
        Factorial = Factorial * I
      End Do
      End Function

      Subroutine systMemUsage()
      Implicit None
      character(len=80) :: line
      integer ::  ios, fu, valueRSS
      valueRSS = -1   
      open(newunit=fu, file='/proc/self/status', action='read')
      do
        read(fu, '(a)',iostat=ios ) line
        if(ios /=0) exit
        if(line(1:6) == 'VmRSS:') then
            read(line(7:), *) valueRSS
            exit
        endif
      enddo
      close(fu)
      Print *, "mem usage(GB)=", valueRSS / (1024 * 1024)
      end Subroutine


