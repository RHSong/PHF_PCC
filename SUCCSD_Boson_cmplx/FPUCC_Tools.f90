
   Module FPUCC_Tools 
   Use Precision
   Use Constants
   Use IntTrans
   Use UHF
   Use CCRes 
   Use PCCCI
   Implicit None
   Logical :: SetPUCCTool = .False.

   Contains


!===============================================================!
!               Build V1 and derivatives of V1                  !
!===============================================================!

      Subroutine BuildV1(V1,Olap,R,NOcc,NSO)
! Build dV1 / dOmega, where Omega is one of the three Euler angles
! Euler is a triple containing the three Euler angles
! dir determines which direction to be differentiated
      Implicit None
      Integer,           Intent(In)  :: NOcc, NSO
      Complex (Kind=pr), Intent(In)  :: R(NSO,NSO)
      Complex (Kind=pr), Intent(Out) :: V1(NOcc,NOcc+1:NSO), Olap
      Complex (Kind=pr), Allocatable :: RotMat(:,:), tmp(:,:), U1(:,:)
      Integer                        :: IAlloc, I
! Allocate
      Allocate(RotMat(NSO,NSO), tmp(NSO,NSO), U1(NOcc+1:NSO,NOcc),      &
               Stat=IAlloc)
      If(IAlloc /= 0) Stop "Could not allocate in BuildV1"
! Build V1
      tmp = Zero
      Do I = 1, NSO
        tmp(I,I) = One
      EndDo
! <phi|e^V1 = <phi|R <=> R^dag|phi> = e^V1dag|phi>
      RotMat = Transpose(Conjg(R))
      Call BuildThoulessZ(NOcc,NSO,tmp,RotMat,Olap,U1)
      V1 = Transpose(Conjg(U1))
      Olap = Conjg(Olap)
! Deallocate
      Deallocate(RotMat, tmp, U1,   &
                 Stat=IAlloc)
      Return
      End Subroutine BuildV1


      Subroutine BuildThoulessZ(NOcc,NSO,MO1,MO2,Olap,Zph)
      Implicit None
      Integer,           Intent(In)    :: NOcc, NSO
      Complex (Kind=pr), Intent(In)    :: MO1(NSO,NSO)
      Complex (Kind=pr), Intent(In)    :: MO2(NSO,NSO)
      Complex (Kind=pr), Intent(Out)   :: Olap
      Complex (Kind=pr), Intent(Out)   :: Zph(NOcc+1:NSO,1:NOcc)
! Work
      Complex (Kind=pr), Allocatable   :: L(:,:), Y(:,:)
      Complex (Kind=pr), Allocatable   :: InvL(:,:)
      Complex (Kind=pr), Allocatable   :: Tran(:,:)
      Integer                          :: IAlloc

! |2> = e^Z |1> * Olap
! Obtain Thouless transformation coefficients
      Allocate(L(NOcc,NOcc),            &
               Y(NOcc+1:NSO,1:NOcc),    &
               Tran(NSO,NSO),           &
               InvL(NOcc,NOcc),         &
               Stat=IAlloc)
      If(IAlloc /= 0) Stop 'Could not allocate in BuildThoulessZ'
                
      Tran = MatMul(Transpose((Conjg(MO1))),MO2)
!     Call OutMat(6,NSO,NSO,1,Tran,"MO of det2 in the basis of det1")
      L    = Tran(1:NOcc,1:NOcc)
      Y    = Tran(NOcc+1:NSO,1:NOcc)

      Call InvertC(L,InvL,NOcc)
      Zph  = MatMul(Y,InvL)
! Overlap
      Call determinantC(Olap,L,NOcc)

! Deallocate
      Deallocate(L, Y, Tran, InvL, Stat=IAlloc)
      If(IAlloc /= 0) Stop 'Could not allocate in BuildThoulessZ'

      Return
      End Subroutine BuildThoulessZ


!===========================================!
!               Utilities                   !
!===========================================!


      Subroutine OrthAO(Ref,Olap,NAO,NSO)
      Implicit None
      Integer,             Intent(In)    :: NAO, NSO
      Complex (Kind=pr),   Intent(In)    :: Olap(NSO,NSO) 
      Complex (Kind=pr),   Intent(InOut) :: Ref(NSO,NSO)
      Complex (Kind=pr)    :: OlapAO(NAO,NAO), Evec(NAO,NAO)
      Real (Kind=pr)       :: Eval(NAO) 
      Complex (Kind=pr)    :: Trans(NSO,NSO)
      Integer              :: I
!     Call Outmat(6,NSO,NSO,1,Olap,"Overlap matrix")
      OlapAO = Olap(1:NAO,1:NAO)
      Call DiagC(OlapAO,Eval,Evec,NAO) 
!     Call Outmat(6,NAO,1,1,Eval,"eval of overlap matrix")
      If (Any(Eval.lt.Zero)) Stop 'Negative e-value in overlap matrix'
      Evec = Conjg(Transpose(Evec))
      Do I = 1, NAO
        OlapAO(I,:) = sqrt(Eval(I)) * Evec(I,:)
      EndDo 
      OlapAO = Matmul(Conjg(Transpose(Evec)),OlapAO)
!     Call Outmat(6,NAO,NAO,1,OlapAO,"S^1/2")
! Build the AO to Orth AO transformation
      Trans = Zero
      Trans(1:NAO,1:NAO) = OlapAO
      Trans(NAO+1:NSO,NAO+1:NSO) = OlapAO
      Ref = MatMul(Trans,Ref)
      Return
      End Subroutine OrthAO

      !Subroutine ProjSz(R1,R2,NOcc,NAO,NSO)
      !Implicit None
      !Integer,           Intent(In)    :: NOcc, NAO, NSO
      !Complex (Kind=pr), Intent(InOut) :: R1(NOcc+1:NSO,NOcc)
      !Complex (Kind=pr), Intent(InOut) :: R2(NOcc+1:NSO,NOcc+1:NSO,NOcc,NOcc)
	  !Integer :: I, J, A, B
	  !Integer :: Sz(NSO)
	  !Do I = 1, NOcc/2
	  !	Sz(I) = 1
	  !  Sz(I+NOcc/2) = -1
	  !End Do
	  !Do I = NOcc+1, (NSO+NOcc)/2
	  !	Sz(I) = 1
	  !  Sz(I+(NSO-NOcc)/2) = -1
	  !End Do
	  !Do I = 1, NOcc
	  !Do A = NOcc+1, NSO
	  !	If(Sz(I).ne.Sz(A)) R1(A,I) = Zero
	  !Do J = 1, NOcc
	  !Do B = NOcc+1, NSO
	  !	If( (Sz(I)+Sz(J)) .ne. (Sz(A)+Sz(B)) ) R2(A,B,I,J) = Zero
	  !End Do
	  !ENd Do
	  !End Do
	  !End Do
      !End Subroutine ProjSz

	  Subroutine ProjSz(R1,R2,NOcc,NAO,NSO)
      Implicit None
      Integer,           Intent(In)    :: NOcc, NAO, NSO
      Complex (Kind=pr), Intent(InOut) :: R1(NOcc+1:NSO,NOcc)
      Complex (Kind=pr), Intent(InOut) :: R2(NOcc+1:NSO,NOcc+1:NSO,NOcc,NOcc)
      Integer :: NOA, NVA
      Complex (Kind=pr) :: R2tmp(NOcc+1:NSO,NOcc+1:NSO,NOcc,NOcc)
      NOA = NOcc / 2
      NVA = NAO - NOA
      R1(NOcc+1:NOcc+NVA,NOA+1:NOcc) = Zero
      R1(NOcc+1+NVA:NSO,1:NOA) = Zero
      R2tmp = Zero
      R2tmp(NOcc+1:NOcc+NVA,NOcc+1:NOcc+NVA,1:NOA,1:NOA) = R2(NOcc+1:NOcc+NVA,NOcc+1:NOcc+NVA,1:NOA,1:NOA)!aaaa
      R2tmp(NOcc+1+NVA:NSO,NOcc+1+NVA:NSO,NOA+1:NOcc,NOA+1:NOcc) = R2(NOcc+1+NVA:NSO,NOcc+1+NVA:NSO,NOA+1:NOcc,NOA+1:NOcc) !bbbb
      R2tmp(NOcc+1:NOcc+NVA,NOcc+1+NVA:NSO,1:NOA,NOA+1:NOcc) = R2(NOcc+1:NOcc+NVA,NOcc+1+NVA:NSO,1:NOA,NOA+1:NOcc)!abab
      R2tmp(NOcc+1:NOcc+NVA,NOcc+1+NVA:NSO,NOA+1:NOcc,1:NOA) = R2(NOcc+1:NOcc+NVA,NOcc+1+NVA:NSO,NOA+1:NOcc,1:NOA) !abba
      R2tmp(NOcc+1+NVA:NSO,NOcc+1:NOcc+NVA,NOA+1:NOcc,1:NOA) = R2(NOcc+1+NVA:NSO,NOcc+1:NOcc+NVA,NOA+1:NOcc,1:NOA) !baba
      R2tmp(NOcc+1+NVA:NSO,NOcc+1:NOcc+NVA,1:NOA,NOA+1:NOcc) = R2(NOcc+1+NVA:NSO,NOcc+1:NOcc+NVA,1:NOA,NOA+1:NOcc) !baab
      R2 = R2tmp
      End Subroutine ProjSz

      Subroutine ProjSz1(R1,NOcc,NAO,NSO)
      Implicit None
      Integer,           Intent(In)    :: NOcc, NAO, NSO
      Complex (Kind=pr), Intent(InOut) :: R1(NSO,NSO)
      Integer :: NOA, NVA
      NOA = NOcc / 2
      NVA = NAO - NOA
      R1(NOcc+1:NOcc+NVA,NOA+1:NOcc) = Zero ! ab
      R1(NOcc+1+NVA:NSO,1:NOA) = Zero  ! ba
      End Subroutine ProjSz1

      Subroutine ProjSz2(R2,NOcc,NAO,NSO)
      Implicit None
      Integer,           Intent(In)    :: NOcc, NAO, NSO
      Complex (Kind=pr), Intent(InOut) :: R2(NSO,NSO,NSO,NSO)
      Complex (Kind=pr) :: R2tmp(NOcc+1:NSO,NOcc+1:NSO,NOcc,NOcc)
      R2tmp = R2(NOcc+1:NSO,NOcc+1:NSO,:NOcc,:NOcc)
      Call ProjSz2VO(R2tmp,NOcc,NAO,NSO)
      R2(NOcc+1:NSO,NOcc+1:NSO,:NOcc,:NOcc) = R2tmp
      End Subroutine ProjSz2

      !Subroutine ProjSz2VO(R2,NOcc,NAO,NSO)
      !Implicit None
      !Integer,           Intent(In)    :: NOcc, NAO, NSO
      !Complex (Kind=pr), Intent(InOut) :: R2(NOcc+1:NSO,NOcc+1:NSO,NOcc,NOcc)
      !Integer :: I, J, A, B
	  !Integer :: Sz(NSO)
	  !Do I = 1, NOcc/2
	  !	Sz(I) = 1
	  !  Sz(I+NOcc/2) = -1
	  !End Do
	  !Do I = NOcc+1, (NSO+NOcc)/2
	  !	Sz(I) = 1
	  !  Sz(I+(NSO-NOcc)/2) = -1
	  !End Do
	  !Do I = 1, NOcc
	  !Do J = 1, NOcc
	  !Do A = NOcc+1, NSO
	  !Do B = NOcc+1, NSO
	  !	If( (Sz(I)+Sz(J)) .ne. (Sz(A)+Sz(B)) ) R2(A,B,I,J) = Zero
	  !End Do
	  !End Do
	  !End Do
	  !End Do
      !End Subroutine ProjSz2VO

	  Subroutine ProjSz2VO(R2,NOcc,NAO,NSO)
      Implicit None
      Integer,           Intent(In)    :: NOcc, NAO, NSO
      Complex (Kind=pr), Intent(InOut) :: R2(NOcc+1:NSO,NOcc+1:NSO,NOcc,NOcc)
      Integer :: NOA, NVA
      Complex (Kind=pr) :: R2tmp(NOcc+1:NSO,NOcc+1:NSO,NOcc,NOcc)
      NOA = NOcc / 2
      NVA = NAO - NOA
      R2tmp = Zero
      R2tmp(NOcc+1:NOcc+NVA,NOcc+1:NOcc+NVA,1:NOA,1:NOA) = R2(NOcc+1:NOcc+NVA,NOcc+1:NOcc+NVA,1:NOA,1:NOA)!aaaa
      R2tmp(NOcc+1+NVA:NSO,NOcc+1+NVA:NSO,NOA+1:NOcc,NOA+1:NOcc) = R2(NOcc+1+NVA:NSO,NOcc+1+NVA:NSO,NOA+1:NOcc,NOA+1:NOcc) !bbbb
      R2tmp(NOcc+1:NOcc+NVA,NOcc+1+NVA:NSO,1:NOA,NOA+1:NOcc) = R2(NOcc+1:NOcc+NVA,NOcc+1+NVA:NSO,1:NOA,NOA+1:NOcc)!abab
      R2tmp(NOcc+1:NOcc+NVA,NOcc+1+NVA:NSO,NOA+1:NOcc,1:NOA) = R2(NOcc+1:NOcc+NVA,NOcc+1+NVA:NSO,NOA+1:NOcc,1:NOA) !abba
      R2tmp(NOcc+1+NVA:NSO,NOcc+1:NOcc+NVA,NOA+1:NOcc,1:NOA) = R2(NOcc+1+NVA:NSO,NOcc+1:NOcc+NVA,NOA+1:NOcc,1:NOA) !baba
      R2tmp(NOcc+1+NVA:NSO,NOcc+1:NOcc+NVA,1:NOA,NOA+1:NOcc) = R2(NOcc+1+NVA:NSO,NOcc+1:NOcc+NVA,1:NOA,NOA+1:NOcc) !baab
      R2 = R2tmp
      End Subroutine ProjSz2VO


      Subroutine ProjOrthZ(R0,R1,R2,Z0,Z1,Z2,NOcc,NSO)
! Make R0, R1, R2 orthogonal to Z0, Z1, Z2
      Implicit None
      Integer,           Intent(In)    :: NOcc, NSO
      Complex (Kind=pr), Intent(In)    :: Z0, Z1(NOcc,NOcc+1:NSO)
      Complex (Kind=pr), Intent(In)    :: Z2(NOcc,NOcc,NOcc+1:NSO,NOcc+1:NSO)
      Complex (Kind=pr), Intent(InOut) :: R0, R1(NOcc+1:NSO,NOcc)
      Complex (Kind=pr), Intent(InOut) :: R2(NOcc+1:NSO,NOcc+1:NSO,NOcc,NOcc)
      Complex (Kind=pr)                :: norm, ovl
      Integer   :: I, J, A, B
! Get overlap between R and Z, as well as the norm of Z
      norm = Z0 * Z0
      ovl  = R0 * Z0
      Do I = 1, NOcc
      Do A = NOcc+1, NSO
        norm = norm + Z1(I,A)*Z1(I,A) 
        ovl  = ovl  + Z1(I,A)*R1(A,I) 
      EndDo
      EndDo
      Do I = 1, NOcc
      Do A = NOcc+1, NSO
      Do J = 1, NOcc
      Do B = NOcc+1, NSO
        norm = norm + Z2(I,J,A,B)*Z2(I,J,A,B)/Four 
        ovl  = ovl  + Z2(I,J,A,B)*R2(A,B,I,J)/Four 
      EndDo
      EndDo
      EndDo
      EndDo
! Project R
      R0 = R0 - Z0 * ovl / norm
      Do I = 1, NOcc
      Do A = NOcc+1, NSO
        R1(A,I) = R1(A,I) - Z1(I,A) * ovl / norm
      EndDo
      EndDo
      Do I = 1, NOcc
      Do A = NOcc+1, NSO
      Do J = 1, NOcc
      Do B = NOcc+1, NSO
        R2(A,B,I,J) = R2(A,B,I,J) - Z2(I,J,A,B) * ovl / norm 
      EndDo
      EndDo
      EndDo
      EndDo
      Return
      End Subroutine ProjOrthZ


!===============================================================!
!               Build Kernels                                   !
!===============================================================!



      Subroutine BuildKernSD(Ovl1,Ovl2,EOvl0,EOvl1,EOvl2,   &
                 Fock,ERI,ESCF,U1t,U2t,T1,T2,NOcc,NSO,SP,TrunODE)
      Implicit None
      Integer,           Intent(In)    :: NOcc, NSO, TrunODE,SP
      Complex (Kind=pr), Intent(In)    :: Fock(NSO,NSO)
      Complex (Kind=pr), Intent(In)    :: ERI(NSO,NSO,NSO,NSO)
      Complex (Kind=pr), Intent(In)    :: ESCF 
      Complex (Kind=pr), Intent(In)    :: U1t(NOcc+1:NSO,NOcc)
      Complex (Kind=pr), Intent(In)    :: U2t(NOcc+1:NSO,NOcc+1:NSO,NOcc,NOcc)
      Complex (Kind=pr), Intent(In)    :: T1(NOcc+1:NSO,NOcc) 
      Complex (Kind=pr), Intent(In)    :: T2(NOcc+1:NSO,NOcc+1:NSO,NOcc,NOcc) 
      Complex (Kind=pr), Intent(Out)   :: EOvl0 
      Complex (Kind=pr), Intent(Out)   :: Ovl1(NSO,NSO), EOvl1(NSO,NSO) 
      Complex (Kind=pr), Intent(Out)   :: Ovl2(NSO,NSO,NSO,NSO) 
      Complex (Kind=pr), Intent(Out)   :: EOvl2(NSO,NSO,NSO,NSO) 
      Complex (Kind=pr), Allocatable   :: tmpEx1(:,:)
      Complex (Kind=pr), Allocatable   :: tmpEx2(:,:,:,:)
      Complex (Kind=pr), Allocatable   :: Ex1(:,:)
      Complex (Kind=pr), Allocatable   :: Ex2(:,:,:,:)
      Integer        :: i, j, a, b, p, q, r, s, NAO
      Integer        :: IAlloc
      NAO = NSO / 2

      If((TrunODE.le.0).or.(TrunODE.ge.4)) Stop 'ODE truncation type not recognized'
! Allocate space for everything.
      Allocate(tmpEx2(NOcc+1:NSO,NOcc+1:NSO,NOcc,NOcc),            &
               tmpEx1(NOcc+1:NSO,NOcc), &
               Ex2(NOcc+1:NSO,NOcc+1:NSO,NOcc,NOcc),            &
               Ex1(NOcc+1:NSO,NOcc), &
               Stat=IAlloc)
      If(IAlloc /= 0) Stop 'Could not allocate in BuildKernSD'
! Fill Ovl1
      Ovl1 = Zero
      !$omp parallel do
      Do I = 1, NOcc
        Ovl1(I,I) = One
      End Do
      !$omp end parallel do
      Ovl1(NOcc+1:NSO,:NOcc) = U1t
      if (SP == 2) then
        Call ProjSz1(Ovl1,NOcc,NAO,NSO)
      endif
! Fill Ovl2
      Ovl2 = Zero
      !$omp parallel do
      Do I = 1, NOcc
      Do J = 1, NOcc 
        If(J.ne.I) then
        Ovl2(I,J,I,J) =  One
        Ovl2(I,J,J,I) = -One
        Do A = NOcc+1, NSO
          Ovl2(A,J,I,J) =  Ovl1(A,I)
          Ovl2(J,A,J,I) =  Ovl1(A,I)
          Ovl2(J,A,I,J) = -Ovl1(A,I)
          Ovl2(A,J,J,I) = -Ovl1(A,I)
        EndDo
        EndIf
        Do A = NOcc+1, NSO
        Do B = NOcc+1, NSO
          Ovl2(A,B,I,J) = U2t(A,B,I,J) + U1t(A,I)*U1t(B,J) - U1t(B,I)*U1t(A,J)
        EndDo
        EndDo
      EndDo
      EndDo
      !$omp end parallel do
      if (SP == 2) then
        Call ProjSz2(Ovl2,NOcc,NAO,NSO)
      end if
! Fill EOvl0
!     Call BuildRes0(EOvl0,U1t,U2t,Fock,ERI,NOcc,NSO)
      Call CCEnergy(U1t,U2t,Fock,ERI,EOvl0,NOcc,NSO)
      EOvl0 = EOvl0 + ESCF
! Build residuals
!      Call BuildRes1SD(Ex1,U1t,U2t,Fock,ERI,NOcc,NSO)
!      Call BuildRes2SD(Ex2,U1t,U2t,Fock,ERI,NOcc,NSO)
      Call CCRes12(Fock,ERI,U1t,U2t,tmpEx1,tmpEx2,NOcc,NSO-NOcc,NSO)

!      Call SetUpCCScr(NOcc,NOcc+1,NSO)
!      Call GetG2Ints(U1t,U2t,ERI,NOcc,NSO)
!      Call GetG1(tmpEx1,U1t,U2t,Fock,ERI,NOcc,NSO)
!      Call GetG2(tmpEx2,U1t,U2t,Fock,ERI,NOcc,NSO)
!      Call GetRes2(U1t,U2t,tmpEx1,tmpEx2,Fock,NOcc,NSO)
!      Call ShutDownCCScr
!      Call ResCCSD(Fock,ERI,U1t,U2t,Ex1,Ex2,NOcc,NSO)
! Fill EOvl1
      EOvl1 = Zero
      EOvl1(NOcc+1:NSO,1:NOcc) = tmpEx1
! Fill EOvl2
      EOvl2 = Zero
      EOvl2(NOcc+1:NSO,NOcc+1:NSO,1:NOcc,1:NOcc) = tmpEx2
! Approximated higher contribution
!      If(TrunODE.eq.3) then
!        Call BuildRes_SDtq(tmpEx1,tmpEx2,U1t,U2t,V1,T1,T2,Fock,ERI,NOcc,NSO)
!        EOvl1(NOcc+1:NSO,1:NOcc) = EOvl1(NOcc+1:NSO,1:NOcc) + tmpEx1
!        EOvl2(NOcc+1:NSO,NOcc+1:NSO,1:NOcc,1:NOcc) =  &
!           EOvl2(NOcc+1:NSO,NOcc+1:NSO,1:NOcc,1:NOcc) + tmpEx2
!      EndIf
      tmpEx1 = EOvl1(NOcc+1:NSO,1:NOcc)
      tmpEx2 = EOvl2(NOcc+1:NSO,NOcc+1:NSO,1:NOcc,1:NOcc) 
! finished vvoo block
      EOvl1 = Zero
      EOvl1(NOcc+1:NSO,:NOcc) = tmpEx1 + U1t * EOvl0
      !$omp parallel do
      Do I = 1, NOcc
        EOvl1(I,I) = EOvl0 
      EndDo
      !$omp end parallel do
      if (SP == 2) then
        Call ProjSz1(EOvl1,NOcc,NAO,NSO)
      end if
      EOvl2 = Zero
      !$omp parallel do
      Do I = 1, NOcc
      Do J = 1, NOcc 
        If(J.ne.I) then
        EOvl2(I,J,I,J) =  EOvl0
        EOvl2(I,J,J,I) = -EOvl0
        Do A = NOcc+1, NSO
          EOvl2(A,J,I,J) =  EOvl1(A,I)
          EOvl2(J,A,J,I) =  EOvl1(A,I)
          EOvl2(J,A,I,J) = -EOvl1(A,I)
          EOvl2(A,J,J,I) = -EOvl1(A,I)
        EndDo
        EndIf
        Do A = NOcc+1, NSO
        Do B = NOcc+1, NSO
          EOvl2(A,B,I,J) =          &
              tmpEx2(A,B,I,J)       &
            + Ovl2(A,B,I,J)*EOvl0   &
            + U1t(A,I)*tmpEx1(B,J) - U1t(B,I)*tmpEx1(A,J) &
            + U1t(B,J)*tmpEx1(A,I) - U1t(A,J)*tmpEx1(B,I) 
        EndDo
        EndDo
      EndDo
      EndDo
      !$omp end parallel do
      if (SP == 2) then
        Call ProjSz2(EOvl2,NOcc,NAO,NSO)
      end if
!     Call OutMat(6,NSO**2,NSO**2,2,EOvl2,"EOvl 2")
! Deallocate 
      Deallocate(tmpEx2, tmpEx1, &
                 Stat=IAlloc)
      If(IAlloc /= 0) Stop 'Could not deallocate in BuildKernSD'
      Return 
      End Subroutine BuildKernSD

      Subroutine OrbTransTAmp(MOGues,T1Gues,T2Gues,MO,T1,T2,NOcc,NSO)
      Implicit None
      Integer,           Intent(In)    :: NOcc, NSO
      Complex (Kind=pr), Intent(In)    :: MOGues(NSO,NSO), MO(NSO,NSO) 
      Complex (Kind=pr), Intent(In)    :: T1Gues(NOcc+1:NSO,NOcc)
      Complex (Kind=pr), Intent(In)    :: T2Gues(NOcc+1:NSO,NOcc+1:NSO,NOcc,NOcc)
      Complex (Kind=pr), Intent(Out)   :: T1(NOcc+1:NSO,NOcc)
      Complex (Kind=pr), Intent(Out)   :: T2(NOcc+1:NSO,NOcc+1:NSO,NOcc,NOcc)
      Complex (Kind=pr), Allocatable   :: Tran(:,:), Op1(:,:), Op2(:,:,:,:)
      Integer   :: I, J, A, B, IAlloc
      Allocate(Tran(NSO,NSO), Op1(NSO,NSO), OP2(NSO,NSO,NSO,NSO),  &
               Stat=IAlloc)
      Tran = MatMul(Transpose(Conjg(MOGues)),MO)
      Op1 = Zero
      Op2 = Zero
      Op1(NOcc+1:NSO,1:NOcc) = T1Gues
      Do J = 1, NOcc
      Do I = 1, NOcc
      Do B = NOcc+1, NSO
      Do A = NOcc+1, NSO
        Op2(A,B,I,J) = T2Gues(A,B,I,J)
      EndDo
      EndDo
      EndDo
      EndDo
      Call IntTran2(Op1,Tran,NSO)
      Call IntTran4(Op2,Tran,NSO)
      T1 = Op1(NOcc+1:NSO,1:NOcc) 
      Do J = 1, NOcc
      Do I = 1, NOcc
      Do B = NOcc+1, NSO
      Do A = NOcc+1, NSO
        T2(A,B,I,J) = Op2(A,B,I,J) 
      EndDo
      EndDo
      EndDo
      EndDo
      Deallocate(Tran, Op1, Op2, Stat=IAlloc)
      Return
      End Subroutine OrbTransTAmp


      Subroutine ReOrderMO(Ref,OneH,ERI,NOcc,NAO,NSO)
      Implicit None
      Integer,           Intent(In)    :: NOcc, NAO, NSO
      Complex (Kind=pr), Intent(In)    :: OneH(NSO,NSO)
      Complex (Kind=pr), Intent(In)    :: ERI(NSO,NSO,NSO,NSO)
      Complex (Kind=pr), Intent(InOut) :: Ref(NSO,NSO)
      Real (Kind=pr),    Parameter     :: Tol = 1.0E-6
      Real (Kind=pr)                   :: Occ 
      Integer                          :: I, J, A, B, cnt
      Integer                          :: Sz(NSO), sorted(NSO)
      Complex (Kind=pr)                :: Fock(NSO,NSO)
! Determine wether MO I is alpha or beta
      Do I = 1, NSO
        Occ = Dot_Product(Ref(1:NAO,I),Ref(1:NAO,I))
        If(Occ.gt.Tol) Then
          Sz(I) = 1
        Else
          Sz(I) = -1
        EndIf
      EndDo
      Call MKFock_AO(OneH,ERI,Fock,Ref,NOcc,NSO)
      Fock = MatMul(Transpose(Ref),MatMul(Fock,Ref))
!     Call outmat(6,NSO,NSO,2,Fock,"Fock from mkfock_ao") 
      Do I = 1, NSO
        Write(*,*) "I, e_MO", I, Fock(I,I)
      EndDo
! Sort orbitals
      sorted = 0
      Do I = 1, NOcc
        cnt = 1
        Do J = 1, I-1
          If((Sz(J).eq.Sz(I)).and.(Real(Fock(I,I)).gt.(Real(Fock(J,J))-Tol))) cnt=cnt+1
        EndDo
        If (Sz(I).eq.1) then
          Do J = NOcc/2, cnt+1, -1
            sorted(J) = sorted(J-1)
          EndDo
          sorted(cnt) = I
        Else
          Do J = NOcc, NOcc/2+cnt+1, -1
            sorted(J) = sorted(J-1)
          EndDo
          sorted(cnt+NOcc/2) = I
        EndIf
      EndDo
      Do I = NOcc+1, NSO
        cnt = 1
        Do J = NOcc+1, I-1
          If((Sz(J).eq.Sz(I)).and.(Real(Fock(I,I)).gt.(Real(Fock(J,J))-Tol))) cnt=cnt+1
        EndDo
        If (Sz(I).eq.1) then
          Do J = (NOcc+NSO)/2, NOcc+cnt+1, -1
            sorted(J) = sorted(J-1)
          EndDo
          sorted(cnt+NOcc) = I
        Else
          Do J = NSO, (NOcc+NSO)/2+cnt+1, -1
            sorted(J) = sorted(J-1)
          EndDo
          sorted(cnt+(NOcc+NSO)/2) = I
        EndIf
      EndDo
!      Call OutMat(6,1,NSO,1,Real(sorted,Kind=pr),"sorted MO")
! Reorder orbitals      
      Fock = Zero
      Do I = 1, NSO
        Fock(:,I) = Ref(:,sorted(I))
      EndDo
      Ref = Fock
      Return
      End Subroutine ReOrderMO

      Subroutine BuildU(Z,T1,T2,U0,U1,U2,NOcc,NSO)
      Implicit None
      Integer,             Intent(In)  :: NOcc,NSO
      Complex(kind=pr),    Intent(In)  :: Z(NOcc,NOcc+1:NSO)
      Complex(kind=pr),    Intent(In)  :: T1(NOcc+1:NSO,NOcc)
      Complex(kind=pr),    Intent(In)  :: T2(NOcc+1:NSO,NOcc+1:NSO,NOcc,NOcc)
      Complex(kind=pr),    Intent(Out) :: U0, U1(NOcc+1:NSO,NOcc)
      Complex(kind=pr),    Intent(Out) :: U2(NOcc+1:NSO,NOcc+1:NSO,NOcc,NOcc)
      Integer :: a,b,c,d,i,j,k,l
      U2 = T2
      U1 = T1
      U0 = 0
      Do i = 1,NOcc
      Do a = NOcc+1,NSO
! 1st order T1
        U0 = U0 + T1(a,i)*Z(i,a)
      End Do
      End Do

      Do i = 1,NOcc
      Do j = 1,NOcc
      Do a = NOcc+1,NSO
      Do b = NOcc+1,NSO
        if (a==b .or. i==j) cycle
! 1st order T2
        U1(a,i) = U1(a,i) + 2*T2(a,b,i,j)/4*Z(j,b)
        U1(a,j) = U1(a,j) - 2*T2(a,b,i,j)/4*Z(i,b)
        U0 = U0 + T2(a,b,i,j)/4 * (Z(i,a)*Z(j,b) - Z(j,a)*Z(i,b))
! 2nd order T1^2/2
        U1(a,j) = U1(a,j) - T1(a,i) * T1(b,j) * Z(i,b)
        U0 = U0 - T1(a,i) * T1(b,j) * Z(j,a)*Z(i,b)/2
      End Do
      End Do
      End Do
      End Do

      Do i = 1,NOcc
      Do j = 1,NOcc
      Do k = 1,NOcc
      Do a = NOcc+1,NSO
      Do b = NOcc+1,NSO
      Do c = NOcc+1,NSO
! 2nd order T1*T2
        if (a==b .or. a==c .or. b==c) cycle
        if (i==j .or. i==k .or. j==k) cycle
        U0 = U0 - 2*T1(c,k)*T2(a,b,i,j)/4*Z(k,a)*Z(i,c)*Z(j,b)
        U0 = U0 + 2*T1(c,k)*T2(a,b,i,j)/4*Z(k,a)*Z(i,b)*Z(j,c)
        U1(c,i) = U1(c,i) - 2*T1(c,k)*T2(a,b,i,j)/4*Z(k,a)*Z(j,b)
        U1(a,k) = U1(a,k) - 2*T1(c,k)*T2(a,b,i,j)/4*Z(i,c)*Z(j,b)
        U1(c,j) = U1(c,j) + 2*T1(c,k)*T2(a,b,i,j)/4*Z(k,a)*Z(i,b)
        U1(b,k) = U1(b,k) + 2*T1(c,k)*T2(a,b,i,j)/4*Z(i,c)*Z(j,a)
        U2(c,b,i,j) = U2(c,b,i,j) - 2*T1(c,k)*T2(a,b,i,j)/4*Z(k,a)
        U2(a,b,k,j) = U2(a,b,k,j) - 2*T1(c,k)*T2(a,b,i,j)/4*Z(i,c)
      End Do
      End Do
      End Do
      End Do
      End Do
      End Do

      Do i = 1,NOcc
      Do j = 1,NOcc
      Do k = 1,NOcc
      Do l = 1,NOcc
      Do a = NOcc+1,NSO
      Do b = NOcc+1,NSO
      Do c = NOcc+1,NSO
      Do d = NOcc+1,NSO
! 2nd order T2^2/2
        if (a==b .or. a==c .or. a==d .or. b==c .or. b==d .or. c==d) cycle
        if (i==j .or. i==k .or. i==l .or. j==k .or. j==l .or. k==l) cycle
        U0 = U0 - T2(a,b,i,j)/4 * T2(c,d,k,l)/4 * Z(i,a)*Z(l,d)*Z(k,b)*Z(j,c)
        U0 = U0 - T2(a,b,i,j)/4 * T2(c,d,k,l)/4 * Z(i,b)*Z(j,c)*Z(k,d)*Z(l,a)
        U0 = U0 + T2(a,b,i,j)/4 * T2(c,d,k,l)/4 * Z(k,b)*Z(j,c)*Z(i,d)*Z(l,a)
        U1(a,i) = U1(a,i) - T2(a,b,i,j)/4 * T2(c,d,k,l)/4 * Z(j,c)*Z(k,b)*Z(l,d)
        U1(b,k) = U1(b,k) - T2(a,b,i,j)/4 * T2(c,d,k,l)/4 * Z(i,a)*Z(l,d)*Z(j,c)
        U1(a,l) = U1(a,l) + T2(a,b,i,j)/4 * T2(c,d,k,l)/4 * Z(i,d)*Z(j,c)*Z(k,b)
        U1(a,l) = U1(a,l) - T2(a,b,i,j)/4 * T2(c,d,k,l)/4 * Z(i,b)*Z(j,c)*Z(k,d)
        U2(a,d,i,l) = U2(a,d,i,l) - T2(a,b,i,j)/4 * T2(c,d,k,l)/4 * Z(j,c)*Z(k,b)
        U2(a,b,i,k) = U2(a,b,i,k) - T2(a,b,i,j)/4 * T2(c,d,k,l)/4 * Z(l,d)*Z(j,c)
        U2(a,c,i,j) = U2(a,c,i,j) - T2(a,b,i,j)/4 * T2(c,d,k,l)/4 * Z(l,d)*Z(k,b)
        U2(a,b,k,l) = U2(a,b,k,l) + T2(a,b,i,j)/4 * T2(c,d,k,l)/4 * Z(i,d)*Z(j,c)
        U2(a,c,j,l) = U2(a,c,j,l) + T2(a,b,i,j)/4 * T2(c,d,k,l)/4 * Z(i,d)*Z(k,b)
        U2(a,d,k,l) = U2(a,d,k,l) + T2(a,b,i,j)/4 * T2(c,d,k,l)/4 * Z(i,b)*Z(j,c)
        U2(d,c,l,i) = U2(d,c,l,i) + T2(a,b,i,j)/4 * T2(c,d,k,l)/4 * Z(j,a)*Z(k,b)
      End Do
      End Do
      End Do
      End Do
      End Do
      End Do
      End Do
      End Do
      End Subroutine

!      Subroutine ProjT1(Z,T1,W0,W1,NOcc,NSO)
!      Implicit None
!      Integer         ,    Intent(In)  :: NOcc,NSO
!      Complex(kind=pr),    Intent(In)  :: Z(NOcc,NOcc+1:NSO)
!      Complex(kind=pr),    Intent(In)  :: T1(NOcc+1:NSO,NOcc)
!      Complex(kind=pr),    Intent(Out) :: W0, W1(NOcc+1:NSO,NOcc)
!      Complex(kind=pr) :: Prod
!      Complex(kind=pr) :: MO1(NSO,NSO), MO2(NSO,NSO)
!      Integer :: a,i
!!      Prod = 1
!!      Do i = 1, NOcc
!!      Do a = NOcc+1, NSO
!!        Prod = Prod * (1 + Z(i,a)*T1(a,i))
!!        W1(a,i) = T1(a,i) / (1 + Z(i,a)*T1(a,i))
!!      End Do
!!      End Do
!      MO1 = Zero
!      MO2 = Zero
!      Do i = 1, NSO
!        MO1(i,i) = One
!        MO2(i,i) = One
!      End Do
!      MO2(1:NOcc,1:NOcc) = MO2(1:NOcc,1:NOcc) + MatMul(Z,T1)
!      MO2(NOcc+1:NSO,1:NOcc) = T1
!      MO2(1:NOcc,NOcc+1:NSO) = Z
!      Call BuildThoulessZ(NOcc,NSO,MO1,MO2,Prod,W1)
!      W0 = Prod
!      End Subroutine

	  Subroutine ProjT1(R,T1,W0,W1,NOcc,NSO)
	  Implicit None
	  Integer         ,    Intent(In)  :: NOcc,NSO
	  Complex(kind=pr),    Intent(In)  :: R(NSO,NSO)
	  Complex(kind=pr),    Intent(In)  :: T1(NOcc+1:NSO,NOcc)
	  Complex(kind=pr),    Intent(Out) :: W0, W1(NOcc+1:NSO,NOcc)
	  Complex(kind=pr) :: Prod
	  Complex(kind=pr) :: MO1(NSO,NSO), MO2(NSO,NSO)
	  Integer :: a,i
	  MO1 = Zero
	  MO2 = Zero
	  Do i = 1, NSO
	  	MO1(i,i) = One
		MO2(i,i) = One
	  End Do
	  MO2(NOcc+1:NSO,1:NOcc) = T1
	  MO2 = matmul(R, MO2)
	  Call BuildThoulessZ(NOcc,NSO,MO1,MO2,Prod,W1)
	  W0 = Prod
	  End Subroutine

      Subroutine DSTT2(R,T1,T2,W0,W1,W2,NOcc,NSO)
      Implicit None
      Integer         ,    Intent(In)  :: NOcc,NSO
      Complex(kind=pr),    Intent(In)  :: R(NSO,NSO)
      Complex(kind=pr),    Intent(In)  :: T1(NOcc+1:NSO,NOcc)
      Complex(kind=pr),    Intent(In)  :: T2(NOcc+1:NSO,NOcc+1:NSO,NOcc,NOcc)
      Complex(kind=pr),    Intent(Out) :: W0, W1(NOcc+1:NSO,NOcc)
      Complex(kind=pr),    Intent(Out) :: W2(NOcc+1:NSO,NOcc+1:NSO,NOcc,NOcc)
      Complex(kind=pr) :: U0, U1(NOcc+1:NSO,NOcc)
      Complex(kind=pr) :: c0, c1(NOcc+1:NSO,NOcc)
      Complex(kind=pr) :: c2(NOcc+1:NSO,NOcc+1:NSO,NOcc,NOcc)
      Integer :: A, I
      Call ProjT1(R,T1,U0,U1,NOcc,NSO)
      Call BuildW(R,U1,T2,W0,W1,W2,NOcc,NSO)
!      Call Kernel4W(Z,U1,T2,W0,W1,W2,NOcc,NSO)
!      Call STC1(Z,U1,T2,c0,c1,c2,NOcc,NSO-NOcc)
!      Call CI2CC(c0,c1,c2,w0,w1,w2,NOcc,NSO)

!      Call STT2(Z,U1,T2,W0,W1,W2,NOcc,NSO-NOcc)
!      W2 = W2 * 4
!      Call ASymm2(W2,NOcc,NSO)

!      Call BuildW2(Z,U1,T2,W0,W1,W2,NOcc,NSO-NOcc)
      W0 = W0 * U0
      W1 = W1 + U1
      End Subroutine

	  Subroutine SemiCanonTrans(Fock,H1,H2,T1,T2,NOcc,NSO)
	  Implicit None
	  Integer          , Intent(In)  :: NOcc,NSO
	  Complex (Kind=pr), Intent(InOut)  :: T1(NOcc+1:NSO,NOcc)
	  Complex (Kind=pr), Intent(InOut)  :: T2(NOcc+1:NSO,NOcc+1:NSO,NOcc,NOcc)
	  Complex (Kind=pr), Intent(InOut)  :: Fock(NSO,NSO), H1(NSO,NSO)
	  Complex (Kind=pr), Intent(InOut)  :: H2(NSO,NSO,NSO,NSO)
	  Complex (Kind=pr) :: tmp(NSO,NSO,NSO,NSO)
	  Complex (Kind=pr) :: MO(NSO,NSO), cMO(NSO,NSO)
	  Complex (Kind=pr) :: Roo(NOcc,NOcc), Rvv(NSO-NOcc,NSO-NOcc)
	  Integer :: I
	  Call IDMat(MO,NSO)
	  Call SemiCanon(H1, H2, MO, NOcc, NSO, cMO, Roo, Rvv)
	  Call IntTran2(H1,cMO,NSO)
	  Call IntTran4(H2,cMO,NSO)
	  Call MKFock_AO(H1,H2,Fock,MO,NOcc,NSO)
	  T1 = MatMul(conjg(Transpose(Rvv)), T1)
	  T1 = MatMul(T1, Roo)
	  tmp = zero
	  tmp(NOcc+1:NSO,NOcc+1:NSO,1:NOcc,1:NOcc) = T2
	  Call IntTran4(tmp,cMO,NSO)
	  T2 = tmp(NOcc+1:NSO,NOcc+1:NSO,1:NOcc,1:NOcc)
	  End Subroutine

!       Subroutine BuildExactKern(Ovl0,Ovl1,Ovl2,EOvl0,EOvl1,EOvl2,V1,  &
!                                 Fock,ERI,T1,T2,NDet,NOcc,NOccA,NOccB,NAO,NSO)
!       Implicit None
!       Integer,        Intent(In)  :: NDet, NOcc, NOccA, NOccB, NAO, NSO
!       Complex (Kind=pr), Intent(In)  :: V1(NOcc,NOcc+1:NSO)
!       Complex (Kind=pr), Intent(In)  :: T1(NOcc+1:NSO,NOcc)
!       Complex (Kind=pr), Intent(In)  :: T2(NOcc+1:NSO,NOcc+1:NSO,NOcc,NOcc)
!       Complex (Kind=pr), Intent(In)  :: Fock(NSO,NSO)
!       Complex (Kind=pr), Intent(In)  :: ERI(NSO,NSO,NSO,NSO)
!       Complex (Kind=pr), Intent(Out) :: Ovl0,  Ovl1(NSO,NSO)
!       Complex (Kind=pr), Intent(Out) :: EOvl0, EOvl1(NSO,NSO)
!       Complex (Kind=pr), Intent(Out) :: Ovl2(NSO,NSO,NSO,NSO)
!       Complex (Kind=pr), Intent(Out) :: EOvl2(NSO,NSO,NSO,NSO)
!       Complex (Kind=pr), Allocatable :: T1full(:,:), V1Full(:,:)
!       Complex (Kind=pr), Allocatable :: Psi1(:), Psi2(:), Psi3(:), Psi4(:)
!       Complex (Kind=pr), Allocatable :: tmp1(:,:), tmp2(:,:,:,:)
! ! Spin integration stuff
!       Integer                     :: IAlloc
!       Integer                     :: I, J, A, B 
! ! Allocate
!       Allocate(Psi1(NDet), Psi2(NDet), Psi3(NDet),          &
!                tmp1(NOcc+1:NSO,NOcc), Psi4(NDet),           &
!                tmp2(NOcc+1:NSO,NOcc+1:NSO,NOcc,NOcc),       &
!                T1full(NSO,NSO), V1full(NSO,NSO), &
!                Stat=IAlloc)
!       If(IAlloc /= 0) Stop "Could not allocate in BuildExactKern"
! ! Build <ex| R e^(U2+U1) |phi> 
!       T1full = Zero
!       T1full(NOcc+1:NSO,1:NOcc) = T1
!       V1Full = Zero
!       V1Full(1:NOcc,NOcc+1:NSO) = V1
!       Psi1    = Zero
!       Psi1(1) = One
!       ! exp(U1) |phi>
!       Call GetExpT1FullPsi(Psi1,Psi2,T1full,NDet,NOccA,NOccB,NOcc,NAO,NSO)  
!       ! exp(U1+U2) |phi>
!       Call GetExpT2Psi(Psi2,Psi1,T2,NDet,NOccA,NOccB,NOcc,NAO,NSO)  
!       ! exp(V) exp(U1+U2) |phi>
!       Call GetExpT1FullPsi(Psi1,Psi2,V1Full,NDet,NOccA,NOccB,NOcc,NAO,NSO)
!       ! H exp(V) exp(U1+U2) |phi>
!       Call BuildHC(Psi2,Psi1,Fock,ERI,NOcc,NSO,NDet)
! ! Get overlaps
! ! Ovl0 = exp(W0)
!       Ovl0  = Psi2(1)
!       Psi2 = Psi2 / Ovl0
!       Psi1 = Psi1 / Ovl0
!       EOvl0 = Psi1(1)
!       Ovl1 = Zero
!       Ovl2 = Zero
!       EOvl1 = Zero
!       EOvl2 = Zero
!       Do I = 1, NOcc
!        Ovl1(I,I) = One
!        EOvl1(I,I) = EOvl0
!       Do A = NOcc+1, NSO
!         Psi3      = Zero
!         Psi3(1)   = One
!         tmp1      = Zero
!         tmp1(A,I) = One
!         T1full = Zero
!         T1full(NOcc+1:NSO,1:NOcc) = tmp1
!         Call GetT1FullPsi(Psi3,Psi4,T1full,NDet,NOccA,NOccB,NOcc,NAO,NSO)
!          Ovl1(A,I) = dot_product(Psi4,Psi2)
!         EOvl1(A,I) = dot_product(Psi4,Psi1)
!       End Do
!       End Do
!
!
!       Do I = 1, NOcc
!       Do J = 1, NOcc
!        If(J.ne.I) then
!            Ovl2(I,J,I,J) =  One
!            Ovl2(I,J,J,I) = -One
!            EOvl2(I,J,I,J) =  EOvl0
!            EOvl2(I,J,J,I) = -EOvl0
!        Do A = NOcc+1, NSO
!            Ovl2(A,J,I,J) =  Ovl1(A,I)
!            Ovl2(J,A,J,I) =  Ovl1(A,I)
!            Ovl2(J,A,I,J) = -Ovl1(A,I)
!            Ovl2(A,J,J,I) = -Ovl1(A,I)
!            EOvl2(A,J,I,J) =  EOvl1(A,I)
!            EOvl2(J,A,J,I) =  EOvl1(A,I)
!            EOvl2(J,A,I,J) = -EOvl1(A,I)
!            EOvl2(A,J,J,I) = -EOvl1(A,I)
!        End Do
!        End IF
!       Do A = NOcc+1, NSO
!       Do B = NOcc+1, NSO
!         Psi3      = Zero
!         Psi3(1)   = One
!         tmp2      = Zero
!         tmp2(A,B,I,J) = One
!         Call ASymm2(tmp2,NOcc,NSO)
!         Call GetT2Psi(Psi3,Psi4,tmp2,NDet,NOccA,NOccB,NOcc,NAO,NSO)
!          Ovl2(A,B,I,J) = dot_product(Psi4,Psi2)
!         EOvl2(A,B,I,J) = dot_product(Psi4,Psi1)
!       EndDo
!       EndDo
!       EndDo
!       EndDo
!!       EOvl0 = EOvl0 / Ovl0
!!       EOvl1 = EOvl1 / Ovl0
!!       EOvl2 = EOvl2 / Ovl0
!!       Ovl1 = Ovl1 / Ovl0
!!       Ovl2 = Ovl2 / Ovl0
! ! Deallocate
!       Deallocate(Psi1, Psi2, Psi3, Psi4, tmp1, tmp2, T1full, &
!                  Stat=IAlloc)
!       If(IAlloc /= 0) Stop "Could not deallocate in BuildExactKern"
!       Return
!       End Subroutine BuildExactKern

    End Module FPUCC_Tools
