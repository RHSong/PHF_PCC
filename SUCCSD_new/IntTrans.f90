
   Module IntTrans
   Use Precision
   Use Constants
   Use omp_lib
   Implicit None
   Integer                        :: nblk
   Integer                        :: blksize = 2

   Contains


      Subroutine IntTran2(Matrix,Evec,NSO)
      Implicit None
      Integer,           Intent(In)    :: NSO
      Complex (Kind=pr), Intent(In)    :: Evec(NSO,NSO)
      Complex (Kind=pr), Intent(InOut) :: Matrix(NSO,NSO)
      Complex (Kind=pr), Allocatable   :: Scr(:,:)
      Integer :: IAlloc

!=============================================!
!  Does the 2-index integral transformation.  !
!=============================================!

!     Write(6,*) "Doing 2-index transformation..."
      Allocate(Scr(NSO,NSO), Stat=IAlloc)
      If(IAlloc /= 0) Stop "Could not allocate in IntTran2"

      Scr = MatMul(Matrix,Evec)
      Matrix = MatMul(Transpose(Conjg(Evec)),Scr)

      DeAllocate(Scr, Stat=IAlloc)
      If(IAlloc /= 0) Stop "Could not deallocate in IntTran2"

      Return
      End Subroutine IntTran2


      Subroutine IntTran4(ERI,Evec,NSO)
      Integer,           Intent(In)    :: NSO
      Complex (Kind=pr), Intent(In)    :: Evec(NSO,NSO)
      Complex (Kind=pr), Intent(InOut) :: ERI(NSO,NSO,NSO,NSO)
      Complex (Kind=pr), Allocatable   :: ERI2(:,:,:,:)
      Integer :: IAlloc, P, Q, NSO2

!=============================================!
!  Does the 4-index integral transformation.  !
!  <pq|rs> = sum C*(mu,p) C*(nu,q) <mu nu|lam sig> C(lam,r) C(sig,s)
!  All four passes are ZGEMMs (the old code did passes 2 and 3 with scalar  !
!  Dot_Product loops over strided NSO^4 sections).                           !
!=============================================!

      Allocate(ERI2(NSO,NSO,NSO,NSO), Stat=IAlloc)
      If(IAlloc /= 0) Stop "Could not allocate in IntTran4"
      NSO2 = NSO*NSO
! index 4: ERI2(mu,nu,lam,s) = sum_sig ERI(mu,nu,lam,sig) C(sig,s)
      ERI2 = Zero
      Call PZGEMM('N','N',NSO*NSO2,NSO,NSO,ZOne,ERI,NSO*NSO2,Evec,NSO,ZOne,ERI2,NSO*NSO2)
! index 3: ERI(mu,nu,r,s) = sum_lam ERI2(mu,nu,lam,s) C(lam,r)
      !$omp parallel do schedule(static) private(P)
      Do P = 1, NSO
        Call ZGEMM('N','N',NSO2,NSO,NSO,ZOne,ERI2(1,1,1,P),NSO2,Evec,NSO,ZZero,ERI(1,1,1,P),NSO2)
      End Do
      !$omp end parallel do
! index 2: ERI2(mu,q,r,s) = sum_nu C*(nu,q) ERI(mu,nu,r,s)   (B = ERI(:,:,r,s) as (mu x nu), times conj(C))
      !$omp parallel do schedule(static) private(P,Q)
      Do P = 1, NSO
      Do Q = 1, NSO
        Call ZGEMM('N','N',NSO,NSO,NSO,ZOne,ERI(1,1,Q,P),NSO,Conjg(Evec),NSO,ZZero,ERI2(1,1,Q,P),NSO)
      End Do
      End Do
      !$omp end parallel do
! index 1: ERI(p,q,r,s) = sum_mu C*(mu,p) ERI2(mu,q,r,s)
      ERI = Zero
      Call PZGEMM('C','N',NSO,NSO*NSO2,NSO,ZOne,Evec,NSO,ERI2,NSO,ZOne,ERI,NSO)

      DeAllocate(ERI2, Stat=IAlloc)
      If(IAlloc /= 0) Stop "Could not deallocate in IntTran4"
      Return
      End Subroutine IntTran4

      Subroutine SaveT1(T1,T1Chk,Evec,Xinv,NOcc,NSO)
      Implicit None
      Integer,           Intent(In)    :: NOcc,NSO
      Complex (Kind=pr), Intent(In)    :: Evec(NSO,NSO)
      Complex (Kind=pr), Intent(In)    :: T1(NOcc+1:NSO,NOcc)
      Complex (Kind=pr), Intent(In)    :: Xinv(NSO,NSO)
      Complex (Kind=pr), Intent(Out)   :: T1Chk(NSO,NSO)
      Complex (Kind=pr) :: tmpVec(NSO,NSO)
      tmpVec = MatMul(Transpose(Conjg(Evec)),Xinv)
      T1Chk = 0
      T1Chk(NOcc+1:NSO,1:NOcc) = T1(NOcc+1:NSO,1:NOcc)
      Call IntTran2(T1Chk,tmpVec,NSO)
      End Subroutine SaveT1

      Subroutine SaveT2(T2,T2Chk,Evec,Xinv,NOcc,NSO)
      Implicit None
      Integer,           Intent(In)    :: NOcc,NSO
      Complex (Kind=pr), Intent(In)    :: Evec(NSO,NSO)
      Complex (Kind=pr), Intent(In)    :: T2(NOcc+1:NSO,NOcc+1:NSO,NOcc,NOcc)
      Complex (Kind=pr), Intent(In)    :: Xinv(NSO,NSO)
      Complex (Kind=pr), Intent(Out)   :: T2Chk(NSO,NSO,NSO,NSO)
      Complex (Kind=pr) :: tmpVec(NSO,NSO)
      tmpVec = MatMul(Transpose(Conjg(Evec)),Xinv)
      T2Chk = 0
      T2Chk(NOcc+1:NSO,NOcc+1:NSO,1:NOcc,1:NOcc) = T2(NOcc+1:NSO,NOcc+1:NSO,1:NOcc,1:NOcc)
      Call IntTran4(T2Chk,tmpVec,NSO)
      End Subroutine SaveT2

      Subroutine ReadT1(T1,T1Chk,Evec,X,NOcc,NSO)
      Implicit None
      Integer,           Intent(In)    :: NOcc,NSO
      Complex (Kind=pr), Intent(In)    :: Evec(NSO,NSO)
      Complex (Kind=pr), Intent(Out)   :: T1(NOcc+1:NSO,NOcc)
      Complex (Kind=pr), Intent(In)    :: T1Chk(NSO,NSO)
      Complex (Kind=pr), Intent(In)    :: X(NSO,NSO)
      Complex (Kind=pr) :: tmp(NSO,NSO), tmpVec(NSO,NSO)
      tmp = T1Chk
      tmpVec = MatMul(X,Evec)
      Call IntTran2(tmp,tmpVec,NSO)
      T1(NOcc+1:NSO,1:NOcc) = tmp(NOcc+1:NSO,1:NOcc)
      End Subroutine ReadT1

      Subroutine ReadT2(T2,T2Chk,Evec,X,NOcc,NSO)
      Implicit None
      Integer,           Intent(In)    :: NOcc,NSO
      Complex (Kind=pr), Intent(In)    :: Evec(NSO,NSO)
      Complex (Kind=pr), Intent(Out)   :: T2(NOcc+1:NSO,NOcc+1:NSO,NOcc,NOcc)
      Complex (Kind=pr), Intent(In)    :: T2Chk(NSO,NSO,NSO,NSO)
      Complex (Kind=pr), Intent(In)    :: X(NSO,NSO)
      Complex (Kind=pr) :: tmp(NSO,NSO,NSO,NSO), tmpVec(NSO,NSO)
      tmp = T2Chk
      tmpVec = MatMul(X,Evec)
      Call IntTran4(tmp,tmpVec,NSO)
      T2(NOcc+1:NSO,NOcc+1:NSO,1:NOcc,1:NOcc) = tmp(NOcc+1:NSO,NOcc+1:NSO,1:NOcc,1:NOcc)
      End Subroutine ReadT2


      Subroutine Mulliken2Dirac(ERI,NSO)
      Integer,           Intent(In)    :: NSO
      Complex (Kind=pr), Intent(InOut) :: ERI(NSO,NSO,NSO,NSO)
      Complex (Kind=pr), Allocatable   :: ERI2(:,:,:,:)
      Integer :: P, Q, R, S
      Integer :: IAlloc
! ERI at input is in Mulliken notation (PQ;RS) = v^PR_QS
! At output it is put in Dirac notation and antisymmetrized
      Allocate(ERI2(NSO,NSO,NSO,NSO), Stat=IAlloc)
      If(IAlloc /= 0) Stop "Could not allocate in Mulliken2Dirac"
      Do S = 1,NSO
      Do R = 1,NSO
      Do Q = 1,NSO
      Do P = 1,NSO
        ERI2(P,R,Q,S) = ERI(P,Q,R,S)
      End Do
      End Do
      End Do
      End Do
      Do S = 1,NSO
      Do R = 1,NSO
      Do Q = 1,NSO
      Do P = 1,NSO
        ERI(P,Q,R,S) = ERI2(P,Q,R,S) - ERI2(P,Q,S,R)
      End Do
      End Do
      End Do
      End Do
! Deallocate and exit safely
      DeAllocate(ERI2, Stat=IAlloc)
      If(IAlloc /= 0) Stop "Could not deallocate in Mulliken2Dirac"
      Return
      End Subroutine Mulliken2Dirac



      Subroutine IntTranEx4(ERI,T1,NOcc,NSO)
! Transform 2-el integral with exp(-T1) ERI exp(T1)
      Integer,           Intent(In)    :: NOcc, NSO
      Complex (Kind=pr), Intent(In)    :: T1(NOcc+1:NSO,NOcc)
      Complex (Kind=pr), Intent(InOut) :: ERI(NSO,NSO,NSO,NSO)
      Complex (Kind=pr), Allocatable   :: ERI1(:,:,:,:)
      Complex (Kind=pr), Allocatable   :: R1(:,:), R2(:,:), R1T(:,:)
      Integer :: IAlloc
      Integer :: Mu, Nu, Lam, Sig, P, Q, R, S, m, m1
	  integer, allocatable :: workload(:), istart(:), iend(:)
	  Integer         :: cnt1, cnt2, clock_rate, clock_max
      nblk = NSO / blksize
      Allocate(workload(nblk), istart(nblk), iend(nblk))
      workload = nso / nblk
      workload(1:mod(nso,nblk)) = workload(1:mod(nso,nblk)) + 1
      do m = 1, nblk
        istart(m) = sum(workload(1:m-1))
        iend(m) = istart(m) + workload(m)
      end do

!=============================================!
!  Does the 4-index integral transformation.  !
!  Important Note: Dot_Product(A,B) = A*.B    !
!=============================================!

      Allocate(ERI1(NSO,NSO,NSO,NSO), R1T(NSO,NSO),  &
               R1(NSO,NSO), R2(NSO,NSO), &
               Stat=IAlloc)
      If(IAlloc /= 0) Stop "Could not allocate in IntTran4"

      Call IDMat(R1,NSO)
      Call IDMat(R2,NSO)
      R1(NOcc+1:NSO,1:NOcc) = -T1
      R2(NOcc+1:NSO,1:NOcc) = T1
      R1T = transpose(R1)
      ERI1 = Zero
      !$omp parallel do schedule(static)
      do m = 1, nblk
        Do P  = istart(m)+ 1, iend(m)
            Call A_dot_B(ERI(:,:,P,:),R2,ERI1(:,:,P,:),NSO**2,NSO,NSO)
        End Do
      end do
      !$omp end parallel do
      ERI = Zero

! ERI are stored in Dirac notation at input
! Transform index 4 - in Dirac notation, we're builidng <Mu Nu | Lam S>
! S is transforms as exp(-U1^T) = [ 1 U_ia ]
!                                [ 0   1  ]


! Transform index 3 - in Dirac notation, we're building <Mu Nu | R S>
! R is transforms as exp(-Z^T) = [ 1 U_ia ]
!                                [ 0   1  ]
    !$omp parallel do schedule(static)
    do m = 1, nblk
        Do P = istart(m)+ 1, iend(m)
            Call A_dot_B(ERI1(:,:,:,P),R2,ERI(:,:,:,P),NSO**2,NSO,NSO)
        End Do
    end do
    !$omp end parallel do


! Transform index 2 - in Dirac notation, we're building <Mu Q | R S>
! Q is transforms as exp(Z) = [   1     0  ]
!                             [ -U_ai   1  ]
    ERI1 = zero
    !$omp parallel do schedule(static)
    do m = 1, nblk
    do m1 = 1, nblk
        Do P = istart(m)+ 1, iend(m)
        Do Q = istart(m1)+ 1, iend(m1)
            Call A_dot_B(ERI(:,:,Q,P),R1T,ERI1(:,:,Q,P),NSO,NSO,NSO)
        End Do
        End Do
    end do
    end do
    !$omp end parallel do

    ERI = Zero
    !$omp parallel do
    do m = 1, nblk
    do m1 = 1, nblk
        Do P = istart(m)+ 1, iend(m)
        Do Q = istart(m1)+ 1, iend(m1)
            Call A_dot_B(R1,ERI1(:,:,Q,P),ERI(:,:,Q,P),NSO,NSO,NSO)
        End Do
        End Do
    end do
    end do
    !$omp end parallel do


! Deallocate and exit safely
      DeAllocate(ERI1, R1,R2, Stat=IAlloc)
      If(IAlloc /= 0) Stop "Could not deallocate in IntTran4"

      Return
      End Subroutine IntTranEx4


      Subroutine IntTranDeex4(ERI,Z1,NOcc,NSO)
! Transform 2-el integral with exp(Z1) ERI exp(-Z1)
      Integer,           Intent(In)    :: NOcc, NSO
      Complex (Kind=pr), Intent(In)    :: Z1(NOcc,NOcc+1:NSO)
      Complex (Kind=pr), Intent(InOut) :: ERI(NSO,NSO,NSO,NSO)
      Complex (Kind=pr), Allocatable   :: ERI1(:,:,:,:)
      Integer :: IAlloc
      Integer :: P, Q, R, S
	  integer, allocatable :: workload(:), istart(:), iend(:)
	  integer :: m, m1
	  nblk = NSO / blksize
	  Allocate(workload(nblk), istart(nblk), iend(nblk))
	  workload = nso / nblk
	  workload(1:mod(nso,nblk)) = workload(1:mod(nso,nblk)) + 1
	  do m = 1, nblk
	      istart(m) = sum(workload(1:m-1))
		  iend(m) = istart(m) + workload(m)
	  end do

!=============================================!
!  Does the 4-index integral transformation.  !
!  Important Note: Dot_Product(A,B) = A*.B    !
!=============================================!

      Allocate(ERI1(NSO,NSO,NSO,NSO),   Stat=IAlloc)
      If(IAlloc /= 0) Stop "Could not allocate in IntTran4"

! ( exp(Z1) H exp(-Z1) )* = exp(-Z1*) H* exp(Z1*)
! IntTranEx4 can be used to transform the dagger
      !$omp parallel do
	  do m = 1, nblk
	  do m1 = 1, nblk
      Do Q = istart(m)+ 1, iend(m)
      Do P = istart(m1)+ 1, iend(m1)
        ERI1(P,Q,:,:) = Conjg(ERI(:,:,P,Q))
      EndDo
      EndDo
	  end do
	  end do
      !$omp end parallel do
      Call IntTranEx4(ERI1,Transpose(Conjg(Z1)),NOcc,NSO)
      !$omp parallel do
	  do m = 1, nblk
	  do m1 = 1, nblk
      Do Q = istart(m)+ 1, iend(m)
      Do P = istart(m1)+ 1, iend(m1)
        ERI(P,Q,:,:) = Conjg(ERI1(:,:,P,Q))
      EndDo
      EndDo
	  end do
	  end do
      !$omp end parallel do

! Deallocate and exit safely
      DeAllocate(ERI1, Stat=IAlloc)
      If(IAlloc /= 0) Stop "Could not deallocate in IntTran4"

      Return
      End Subroutine IntTranDeex4



      Subroutine IntTranEx2(OneH,T1,NOcc,NSO)
! Transform 2-el integral with exp(-T1) OneH exp(T1)
      Integer,           Intent(In)    :: NOcc, NSO
      Complex (Kind=pr), Intent(In)    :: T1(NOcc+1:NSO,NOcc)
      Complex (Kind=pr), Intent(InOut) :: OneH(NSO,NSO)
      Complex (Kind=pr), Allocatable   :: R1(:,:)
      Complex (Kind=pr), Allocatable   :: R2(:,:)
      Integer :: IAlloc
      Integer :: Mu, Nu, P, Q

!=============================================!
!  Does the 2-index integral transformation.  !
!  Important Note: Dot_Product(A,B) = A*.B    !
!=============================================!

      Allocate(R1(NSO,NSO),   &
               R2(NSO,NSO),   Stat=IAlloc)
      If(IAlloc /= 0) Stop "Could not allocate in IntTran2"

      Call IDMat(R1,NSO)
      Call IDMat(R2,NSO)
      R1(NOcc+1:NSO,1:NOcc) = -T1
      R2(NOcc+1:NSO,1:NOcc) = T1
      OneH = MatMul(OneH,R2)
      OneH = MatMul(R1,OneH)

! Deallocate and exit safely
      DeAllocate(R1, R2, Stat=IAlloc)
      If(IAlloc /= 0) Stop "Could not deallocate in IntTran2"

      Return
      End Subroutine IntTranEx2


      Subroutine IntTranDeex2(OneH,Z1,NOcc,NSO)
! Transform 2-el integral with exp(Z1) OneH exp(-Z1)
      Integer,           Intent(In)    :: NOcc, NSO
      Complex (Kind=pr), Intent(In)    :: Z1(NOcc,NOcc+1:NSO)
      Complex (Kind=pr), Intent(InOut) :: OneH(NSO,NSO)
      Complex (Kind=pr), Allocatable   :: OneH1(:,:)
      Integer :: IAlloc
      Integer :: P, Q

!=============================================!
!  Does the 2-index integral transformation.  !
!  Important Note: Dot_Product(A,B) = A*.B    !
!=============================================!

      Allocate(OneH1(NSO,NSO),   Stat=IAlloc)
      If(IAlloc /= 0) Stop "Could not allocate in IntTran2"

! ( exp(Z1) H exp(-Z1) )* = exp(-Z1*) H* exp(Z1*)
! IntTranEx2 can be used to transform the dagger
      OneH1 = Transpose(Conjg(OneH))
      Call IntTranEx2(OneH1,Transpose(Conjg(Z1)),NOcc,NSO)
      OneH1 = Transpose(Conjg(OneH1))
      OneH = OneH1

! Deallocate and exit safely
      DeAllocate(OneH1, Stat=IAlloc)
      If(IAlloc /= 0) Stop "Could not deallocate in IntTran2"

      Return
      End Subroutine IntTranDeex2

      Subroutine TransposeHam(OneH,ERI,NSO)
      Implicit None
      Integer,           Intent(In)    :: NSO
      Complex (Kind=pr), Intent(InOut) :: OneH(NSO,NSO)
      Complex (Kind=pr), Intent(InOut) :: ERI(NSO,NSO,NSO,NSO)
      Complex (Kind=pr), Allocatable   :: ERI2(:,:,:,:)
      Integer :: P, Q, R, S, IAlloc
      Allocate(ERI2(NSO,NSO,NSO,NSO), Stat=IAlloc)
      If(IAlloc /= 0) Stop "Could not allocate in TransposeHam"
! This SR transpose one and two electron integrals, putting
! them in the up-left, down right order
      OneH = Transpose(OneH)
      Do P = 1,NSO
      Do Q = 1,NSO
      Do R = 1,NSO
      Do S = 1,NSO
        ERI2(P,Q,R,S) = ERI(R,S,P,Q) 
      End Do
      End Do
      End Do
      End Do
      ERI = ERI2
! Deallocate and exit safely
      DeAllocate(ERI2, Stat=IAlloc)
      Return
      End Subroutine TransposeHam

      Subroutine IntTran2Sim(Matrix,RMat,NSO)
      Implicit None
      Integer,           Intent(In)    :: NSO
      Complex (Kind=pr), Intent(In)    :: RMat(NSO,NSO)
      Complex (Kind=pr), Intent(InOut) :: Matrix(NSO,NSO)
      Complex (Kind=pr), Allocatable   :: Scr(:,:), InvR(:,:)
      Integer :: IAlloc

!=============================================!
!  Does the 2-index similarity transformation.!
!              O' = R O R^-1                  !   
!  i.e., the ket indices transforms as R^-1   !
!  while the bra indices transforms as R      !
!=============================================!

!     Write(6,*) "Doing 2-index transformation..."
      Allocate(Scr(NSO,NSO), InvR(NSO,NSO), Stat=IAlloc)
      If(IAlloc /= 0) Stop "Could not allocate in IntTran2"

      Call InvertC(RMat,InvR,NSO)
      Scr = MatMul(Matrix,InvR)
      Matrix = MatMul(RMat,Scr)

      DeAllocate(Scr, InvR, Stat=IAlloc)
      If(IAlloc /= 0) Stop "Could not deallocate in IntTran2"

      Return
      End Subroutine IntTran2Sim

      Subroutine IntTran4Sim(ERI,RMat,NSO)
      Integer,           Intent(In)    :: NSO
      Complex (Kind=pr), Intent(In)    :: RMat(NSO,NSO)
      Complex (Kind=pr), Intent(InOut) :: ERI(NSO,NSO,NSO,NSO)
      Complex (Kind=pr), Allocatable   :: ERI2(:,:,:,:), InvR(:,:)
      Integer :: IAlloc
      Integer :: Mu, Nu, Lam, Sig, P, Q, R, S
	  integer :: m, m1
	  integer, allocatable :: workload(:), istart(:), iend(:)
	  Integer         :: cnt1, cnt2, clock_rate, clock_max
	  nblk = NSO / blksize
	  Allocate(workload(nblk), istart(nblk), iend(nblk))
	  workload = nso / nblk
	  workload(1:mod(nso,nblk)) = workload(1:mod(nso,nblk)) + 1
	  do m = 1, nblk
	  	istart(m) = sum(workload(1:m-1))
		iend(m) = istart(m) + workload(m)
	  end do


!=============================================!
!  Does the 4-index integral transformation.  !
!  Important Note: Dot_Product(A,B) = A*.B    !
!              O' = R O R^-1                  !   
!  i.e., the ket indices transforms as R^-1   !
!  while the bra indices transforms as R      !
!=============================================!

!     Write(6,*) "Doing 4-index transformation..."
      Allocate(ERI2(NSO,NSO,NSO,NSO), &
               InvR(NSO,NSO), Stat=IAlloc)
      If(IAlloc /= 0) Stop "Could not allocate in IntTran4"
      Call InvertC(RMat,InvR,NSO)
    ERI2 = Zero
!    Call A_dot_B(ERI,InvR,ERI2,NSO**3,NSO,NSO)
    !$omp parallel do
    do m = 1, nblk
        Do P  = istart(m)+ 1, iend(m)
            Call A_dot_B(ERI(:,:,P,:),InvR,ERI2(:,:,P,:),NSO**2,NSO,NSO)
        End Do
    end do
    !$omp end parallel do
    ERI = Zero

! Transform index 3 - in Dirac notation, we're building <Mu Nu | Q P>
    !$omp parallel do schedule(static)
    do m = 1, nblk
        Do P  = istart(m)+ 1, iend(m)
            Call A_dot_B(ERI2(:,:,:,P),InvR,ERI(:,:,:,P),NSO**2,NSO,NSO)
        End Do
    end do
    !$omp end parallel do

! Transform index 2 - in Dirac notation, we're building <Mu R | Q P>
    ERI2 = Zero

    !$omp parallel do schedule(static)
    do m = 1, nblk
    do m1 = 1, nblk
        Do P = istart(m)+ 1, iend(m)
        Do Q = istart(m1)+ 1, iend(m1)
            Call A_dot_B(ERI(:,:,Q,P),transpose(RMat),ERI2(:,:,Q,P),NSO,NSO,NSO)
        End Do
        End Do
    end do
    end do
    !$omp end parallel do

    ERI = Zero
    !$omp parallel do
    do m = 1, nblk
    do m1 = 1, nblk
        Do P = istart(m)+ 1, iend(m)
        Do Q = istart(m1)+ 1, iend(m1)
            Call A_dot_B(RMat,ERI2(:,:,Q,P),ERI(:,:,Q,P),NSO,NSO,NSO)
        End Do
        End Do
    end do
    end do
    !$omp end parallel do

! Deallocate and exit safely
      DeAllocate(ERI2, InvR, Stat=IAlloc)
      If(IAlloc /= 0) Stop "Could not deallocate in IntTran4"

      Return
      End Subroutine IntTran4Sim

      Subroutine ASymmH2(H2,NSO)
      Implicit None
      Integer,           Intent(In)    :: NSO
      Complex (Kind=pr), Intent(InOut) :: H2(NSO,NSO,NSO,NSO)
      Complex (Kind=pr) :: tmp(NSO,NSO,NSO,NSO)
      Integer :: I,J
      !$omp parallel do
      Do I = 1, NSO
      Do J = 1, NSO
        tmp(:,:,I,J) = H2(:,:,I,J) - H2(:,:,J,I)
      End Do
      End Do
      !$omp end parallel do
      !$omp parallel do
      Do I = 1, NSO
      Do J = 1, NSO
        H2(I,J,:,:) = tmp(I,J,:,:) - tmp(J,I,:,:)
      End Do
      End Do
      !$omp end parallel do
      H2 = H2 / 4
      End Subroutine ASymmH2


      Subroutine TransBlockT2(T2,o,v,Lp,nP,Lq,nQ,Mr,nR,Ms,nS,Out)
!=====================================================================!
!  Out(P,Q,R,S) = sum_{a,b,i,j} Lp(P,a) Lq(Q,b) T2(a,b,i,j) Mr(i,R) Ms(j,S)
!  i.e. the (P,Q,R,S) block of the 4-index similarity transform of a
!  tensor that is nonzero only in its vvoo block, built with four ZGEMMs
!  and never touching an NSO^4 array.  The bra indices are contracted
!  first so that no intermediate exceeds max(v^2 o^2, size(Out)).
!=====================================================================!
      Implicit None
      Integer,           Intent(In)  :: o, v, nP, nQ, nR, nS
      Complex (Kind=pr), Intent(In)  :: T2(v,v,o,o)
      Complex (Kind=pr), Intent(In)  :: Lp(nP,v), Lq(nQ,v), Mr(o,nR), Ms(o,nS)
      Complex (Kind=pr), Intent(Out) :: Out(nP,nQ,nR,nS)
      Complex (Kind=pr), Allocatable :: X1(:,:,:,:), X2(:,:,:,:), X3(:,:,:,:)
      Integer :: i, j, S, IAlloc
      Allocate(X1(nP,v,o,o), X2(nP,nQ,o,o), X3(nP,nQ,nR,o), Stat=IAlloc)
      If (IAlloc /= 0) Stop "Could not allocate in TransBlockT2"
! P:  X1(P,b,i,j) = sum_a Lp(P,a) T2(a,b,i,j)
      X1 = Zero
      Call PZGEMM('N','N',nP,v*o*o,v,ZOne,Lp,nP,T2,v,ZOne,X1,nP)
! Q:  X2(P,Q,i,j) = sum_b X1(P,b,i,j) Lq(Q,b)
      X2 = Zero
      !$omp parallel do schedule(static) private(i,j)
      Do j = 1, o
      Do i = 1, o
        Call ZGEMM('N','T',nP,nQ,v,ZOne,X1(1,1,i,j),nP,Lq,nQ,ZZero,X2(1,1,i,j),nP)
      End Do
      End Do
      !$omp end parallel do
! R:  X3(P,Q,R,j) = sum_i X2(P,Q,i,j) Mr(i,R)
      X3 = Zero
      !$omp parallel do schedule(static) private(j)
      Do j = 1, o
        Call ZGEMM('N','N',nP*nQ,nR,o,ZOne,X2(1,1,1,j),nP*nQ,Mr,o,ZZero,X3(1,1,1,j),nP*nQ)
      End Do
      !$omp end parallel do
! S:  Out(P,Q,R,S) = sum_j X3(P,Q,R,j) Ms(j,S)
      Out = Zero
      Call PZGEMM('N','N',nP*nQ*nR,nS,o,ZOne,X3,nP*nQ*nR,Ms,o,ZOne,Out,nP*nQ*nR)
      Deallocate(X1, X2, X3)
      End Subroutine TransBlockT2

   End Module IntTrans

