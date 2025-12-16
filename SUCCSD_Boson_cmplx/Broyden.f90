    Module Broyden
! Implementation of modified Broyden's methods
! Ref: PRC 78, 014318; PRB, 38, 12807 
    Use Precision
    Use Constants
    Use UHF
    Use FPUCC
    Use FPUCC_Tools
!   Use InvSy
    Implicit None
    Private
    Public  :: BroydenIter, GuessVec
! nBrd: # of unknowns, pBrd: # of Broyden step stored
    Integer                        :: nBrd, pBrd
! Local copies of PUCC variables   
    Integer                        :: nsBrd, ndBrd
    Integer                        :: NOccBrd
    Integer                        :: NAOBrd, NSOBrd
    Integer                        :: QJBrd, SPBrd
	Integer                        :: ncispBrd, ncipgBrd, ncikBrd
    Integer                        :: TrunODEBrd, NPointsBrd(2), NpgBrd
    Real   (Kind=pr), Allocatable  :: RootaBrd(:), RootbBrd(:), RootyBrd(:)
    Real   (Kind=pr), Allocatable  :: WeightspBrd(:,:), WeightpgBrd(:,:,:)
    Complex(Kind=pr), Allocatable  :: fspBrd(:),fpgBrd(:),fkBrd(:)
    Complex(Kind=pr), Allocatable  :: HOneBrd(:,:)
    Complex(Kind=pr), Allocatable  :: RefBrd(:,:)
	Complex(Kind=pr), Allocatable  :: R1Brd(:,:,:), R2Brd(:,:,:)
	Complex(Kind=pr), Allocatable  :: RpgBrd(:,:,:), RkBrd(:,:,:)
    Complex(Kind=pr), Allocatable  :: MixVec(:)
    Real (Kind=pr)                 :: ENucBrd
    Complex(Kind=pr)               :: EneBrd
	Real (Kind=pr)                 :: damp = 1.0_pr
    Logical                        :: DoCCD = .False.
    Integer                        :: commBrd

    Contains

    Subroutine BroydenIter(HOne,HTwo,Ref,Xinv,Ene,T1,T2,NOcc,NAO,NSO,ENuc,QJ, &
                           CC,SP,NPoints,Npg, ncisp, ncipg, ncik, R1, R2, Rpg, Rk, &
						   Roota,Rootb,Rooty,Weightsp,Weightpg, fsp,fpg,fk, &
						   TrunODE,pBrdIn,DoCCDIn,comm)
    Implicit None
    Integer,           Intent(In)    :: NOcc, NAO, NSO
    Integer,           Intent(In)    :: NPoints(2), Npg
	Integer,           Intent(In)    :: ncisp, ncipg, ncik
	Integer,           Intent(In)    :: TrunODE, pBrdIn, QJ,CC,SP
    Complex (Kind=pr), Intent(In)    :: HOne(NSO,NSO)
    Complex (Kind=pr), Intent(In)    :: HTwo(NSO,NSO,NSO,NSO)
    Complex (Kind=pr), Intent(In)    :: Ref(NSO,NSO), Xinv(NSO,NSO)
    Real (Kind=pr),    Intent(In)    :: ENuc
    Real (Kind=pr),    Intent(In)    :: Roota(:), Rootb(:), Rooty(:)
    Real (Kind=pr),    Intent(In)    :: Weightsp(:,:), Weightpg(:,:,:)
	Complex (Kind=pr), Intent(In)    :: fsp(ncisp),fpg(ncipg),fk(ncik)
	Complex (Kind=pr), Intent(In)    :: R1(:,:,:), R2(:,:,:)
	Complex (Kind=pr), Intent(In)    :: Rpg(npg,NSO,NSO), Rk(ncik,NSO,NSO)
    Complex (Kind=pr), Intent(InOut) :: T1(NOcc+1:NSO,NOcc)
    Complex (Kind=pr), Intent(InOut) :: T2(NOcc+1:NSO,NOcc+1:NSO,NOcc,NOcc)
    Complex (Kind=pr), Intent(InOut) :: Ene 
    Logical,           Intent(In)    :: DoCCDIn
    Integer,           Intent(In)    :: comm
! CCSD Amp    
    Complex (Kind=pr), Allocatable   :: T1Chk(:,:)
    Complex (Kind=pr), Allocatable   :: T2Chk(:,:,:,:)
! Mixing vector stuff
    Complex (Kind=pr), Allocatable   :: Fock(:,:)
    Complex (Kind=pr)                :: ESCF, Denom 
    Integer                          :: I, J, A, B
! Broyden iteration stuff
    Complex (Kind=pr), Allocatable   :: xold(:), xnew(:), x0(:)
    Complex (Kind=pr), Allocatable   :: fold(:), fnew(:), F0(:)
    Complex (Kind=pr), Allocatable   :: dF(:,:), dx(:,:)
    Complex (Kind=pr), Allocatable   :: dF2(:,:), dx2(:,:)
    Real (Kind=pr),    Allocatable   :: w(:)
    Real (Kind=pr)                   :: w0, mix
    Real (Kind=pr),    Parameter     :: TolMax = 5.0E-5_pr
    Integer,           Parameter     :: CycMax = 300
    Integer                          :: NIter, IAlloc
    Real (Kind=pr)                   :: dT, dTold 
    Real (Kind=pr)                   :: ResNew, ResOld
! DIIS
    Integer                          :: NDIIS = 10
    Complex (Kind=pr), Allocatable   :: ErrVecs(:,:), Veclist(:,:)
	Real (Kind=pr), Allocatable      :: tmp(:)
! timer
    Integer         :: cnt1, cnt2, clock_rate, clock_max
    Real (Kind=pr)  :: cput1, cput2
! Set things up 
    Write(*,*) "Enter Broyden"
    DoCCD = DoCCDIn
    Call SetUpBroyden(HOne,HTwo,Ref,NOcc,NAO,NSO,ENuc,QJ, SP, &
					  NPoints,Npg, ncisp, ncipg, ncik, R1, R2, Rpg, Rk, &
					  Roota,Rootb,Rooty,Weightsp,Weightpg,fsp,fpg,fk,TrunODE,comm)
!   Call SetUpInvSy(Ref,NOcc,NAO,NSO)
    pBrd = pBrdIn
! Allocate
    Allocate(xold(nBrd), xnew(nBrd), fold(nBrd), fnew(nBrd),    &
             Fock(NSO,NSO), tmp(nBrd), x0(nBrd),F0(nBrd), &
             Stat=IAlloc)
    If(IAlloc/=0) Stop "Could not allocate in BroydenIter"
    If(pBrd>0) then
      Allocate(dF(nBrd,pBrd), dx(nBrd,pBrd), w(pBrd),             &
               dF2(nBrd,pBrd), dx2(nBrd,pBrd), &
               Stat=IAlloc)
      If(IAlloc/=0) Stop "Could not allocate in BroydenIter"
    EndIf
    Write(*,*) "Entered PbarHbar Feedback"
    Open(8,File='Output',Position='Append')
    Write(8,1020)
    Write(8,1010)
    Write(8,1020)
    Write(8,1030)
! Initialize 
    If(DoCCD) T1 = Zero ! CCD!
    Call Mat2Vec(xold,Zero*Ene,T1,T2)
    dT = maxval(abs(xold))
!	if (dT < TolMax) then
!		continue
!	else
!		Call BuildSRes(xold)
!		dT = maxval(abs(xold)) / dT
!		xold = xold / dT
!	end if
    Call EvalF(xold,fold,HTwo,NSO)
    dT = sqrt(Dot_Product(fold,fold))
    Call BuildSRes(fold) 
    damp = 2.0_pr * sqrt(Dot_Product(fold,fold)) / dT
    fold = fold / damp
    Print *, "damping=", damp
    NIter = 0
    dT = 1.0_pr
    dF = Zero
    dx = Zero 
    dF2 = Zero
    dx2 = Zero
    w  = One
    w0 = 0.01_pr
    mix= 0.50_pr
    If (CC == 0) then
        fold = Cmplx(Real(fold,kind=pr),Zero,kind=pr)
    EndIF
    Write(8,1040) Real(EneBrd),NIter,dT
    Print *, "E(start)=", EneBrd
! Build the mixing vector
    Call MKFock_MO(HOne,HTwo,Fock,ESCF,NOcc,NSO)
    Do I = 1, NOcc
    Do A = NOcc+1, NSO
      Denom = Fock(I,I) - Fock(A,A)
      T1(A,I) = One/Denom 
    Do J = 1, NOcc
    Do B = NOcc+1, NSO
      Denom = Fock(I,I) + Fock(J,J) - Fock(A,A) - Fock(B,B)
      T2(A,B,I,J) = One/Denom 
    EndDo
    EndDo
    EndDo
    EndDo
    Ene = Zero
    Call Mat2Vec(MixVec,Ene,T1,T2)
!	MixVec = One
!	Call ProjectNull(NullVec,MixVec,NNull)
! Start iteration
    xnew = xold
    Do While(dT >= TolMax)
      NIter = NIter + 1
      Call BroydenStep(xnew,xold,fold,mix,w0,w,dF,dx,NIter)
      Call EvalF(xnew,fnew,HTwo,NSO)
	  Call BuildSRes(fnew)
!      if (NIter .ge. 10) then
!        Call BuildSRes(fnew)
!        if (NIter == 0) then
!            dF = Zero
!            dx = Zero
!        endif
!	  Call ProjectNull(NullVec,fnew,NNull)
!      endif
      If (CC == 0) then
        fnew = Cmplx(Real(fnew,kind=pr),Zero,kind=pr)
      EndIF
      Call UpdateBroyden(dF,dx,fnew,fold,xnew,xold)
!      ResNew = Sqrt(Dot_Product(fnew,fnew)/real(nBrd,kind=pr))
!      dT = ResNew
!      dTold = dT
      dT = MaxVal(Abs(fnew)) * damp
!     dT = MaxVal(Abs(xnew-xold))
      Write(8,1040) Real(EneBrd),NIter,dT
      If (NIter > CycMax) Exit
! Prepare for the next step
      xold = xnew
      fold = fnew
      flush(6)
      flush(8)
! Debug: save T1 T2 at each iteration
      Call Vec2Mat(xnew,Ene,T1,T2)
      Open(70,File="Chk_PUCC_Ene",status="Replace",form="unformatted")
      Write(70) EneBrd 
      Close(70)
      Open(70,File="Chk_PUCC_T1",status="Replace",form="unformatted")
      Write(70) T1
      Close(70)
      Open(70,File="Chk_PUCC_T2",status="Replace",form="unformatted")
      Write(70) T2
      Close(70)
! End of Debug
    EndDo
!	Do While (dT > 1e-5 .and. NIter < 140)
!		Call DoDIISPCC(ErrVecs,Veclist,nbrd,xnew)
!		Call EvalF(xnew,fnew,HTwo,NSO)
!		Call UpdateDIISPCC(xnew,xnew-xold,Veclist,ErrVecs,nbrd)
!		xold = xnew
!		fold = fnew
!		dT = Sqrt(Dot_Product(fnew,fnew)/real(nBrd,kind=pr))
!		Print *, "DIIS dT = ", dT, xnew(1)
!		NIter = NIter + 1
!	End Do
    Call Vec2Mat(xnew,Ene,T1,T2)
!    allocate(T1Chk(NSO,NSO), T2Chk(NSO,NSO,NSO,NSO))
!    Call SaveT1(T1,T1Chk,Ref,Xinv,NOcc,NSO)
!    Call SaveT2(T2,T2Chk,Ref,Xinv,NOcc,NSO)
    Write(*,*) "max|T1|", maxval(abs(T1))
    Write(*,*) "max|T2|", maxval(abs(T2))
    Write(*,*) "max|imag(T1)|", maxval(abs(aimag(T1)))
    Write(*,*) "max|imag(T2)|", maxval(abs(aimag(T2)))
    Open(70,File="Chk_PUCC_T1",status="Replace",form="unformatted")
    Write(70) T1
    Close(70)
    Open(70,File="Chk_PUCC_T2",status="Replace",form="unformatted")
    Write(70) T2
    Close(70)
    Write(8,1020)
    Write(8,1050) NIter
    Write(8,1070) Real(EneBrd)
    Write(8,1080) Aimag(EneBrd)
    Write(8,1000)
    Close(8)
    Write(*,*) "VAP PUCC energy from Broyden:", EneBrd
!    deallocate(T1Chk,T2Chk)
    Ene = EneBrd
! Outputs
1000  Format(14x,'**************************************************')
1010  Format(14X,'*       PbarHbar FeedBack summary follows        *')
1020  Format(14x,'*------------------------------------------------*')
1030  Format(14X,'* PUCCSD Energy    Iteration    Biggest Res      *')
1040  Format(14X,'* ',F15.8,2x,I5,7X,F14.10,4x,'*')
1050  Format(14x,'* PUCC has converged in ',I3,' iterations',11x,'*')
1070  Format(14x,'* Final PUCC Energy is (real) ',F15.9,'a.u.*')
1080  Format(14x,'* Final PUCC Energy is (imag) ',F15.9,'a.u.*')
! deallocate
    Deallocate(xold, xnew, fold, fnew, dF, dx, Fock, &
               Stat=IAlloc)
    If(IAlloc/=0) Stop "Could not deallocate in BroydenIter"
    Call ShutDownBroyden
!    Call ShutDownDIISCC
!   Call ShutDownInvSy
    Return
    End Subroutine BroydenIter

    Subroutine SetUpBroyden(HOne,HTwo,Ref,NOcc,NAO,NSO,ENuc,QJ, SP, &
							NPoints,Npg, ncisp, ncipg, ncik, R1, R2, Rpg, Rk, &
							Roota,Rootb,Rooty,Weightsp,Weightpg,fsp,fpg,fk,TrunODE,comm)
    Implicit None
    Integer,           Intent(In)    :: NOcc, NAO, NSO
    Integer,           Intent(In)    :: NPoints(2), Npg
	Integer,           Intent(In)    :: TrunODE, QJ, SP, comm
	Integer,           Intent(In)    :: ncisp, ncipg, ncik
    Real (Kind=pr),    Intent(In)    :: Roota(:), Rootb(:), Rooty(:)
    Real (Kind=pr),    Intent(In)    :: Weightsp(:,:), Weightpg(:,:,:)
	Complex (Kind=pr), Intent(In)    :: fsp(ncisp),fpg(ncipg),fk(ncik)
	Complex (Kind=pr), Intent(In)    :: R1(:,:,:), R2(:,:,:)
	Complex (Kind=pr), Intent(In)    :: Rpg(npg,NSO,NSO), Rk(ncik,NSO,NSO)
    Complex (Kind=pr), Intent(In)    :: HOne(NSO,NSO)
    Complex (Kind=pr), Intent(In)    :: HTwo(NSO,NSO,NSO,NSO)
    Complex (Kind=pr), Intent(In)    :: Ref(NSO,NSO)
    Real (Kind=pr),    Intent(In)    :: ENuc 
    Integer     :: IAlloc, I, J
! Set up local copies of necessary variables used in derivative evaluation
    NOccBrd  = NOcc
    NAOBrd   = NAO
    NSOBrd   = NSO
    ENucBrd  = ENuc
    QJBrd = QJ
    SPBrd = SP
    TrunODEBrd = TrunODE
    NPointsBrd = NPoints
	NpgBrd = Npg
	ncispBrd = ncisp
	ncipgBrd = ncipg
	ncikBrd = ncik
    I = NPoints(1)
    J = NPoints(2)
! y contains T1, T2
    nsBrd = NOcc*(NSO-NOcc)
    ndBrd = NOcc*(NOcc-1)*(NSO-NOcc)*(NSO-NOcc-1)/4
    nBrd  = 1 + nsBrd + ndBrd 
! Allocate and continue copying variables
    Allocate(HOneBrd(NSO,NSO), &
             RefBrd(NSO,NSO), MixVec(nBrd), R1Brd(I,NSO,NSO),        &
			 R2Brd(J,NSO,NSO),RpgBrd(Npg,NSO,NSO),RkBrd(ncik,NSO,NSO), &
			 fspBrd(ncisp),fpgBrd(ncipg),fkBrd(ncik), WeightpgBrd(Npg,ncipg,ncipg), &
             RootaBrd(I), RootbBrd(I), RootyBrd(J), WeightspBrd(I,J), &
             Stat=IAlloc)
    If(IAlloc /= 0) Stop "Could not allocate in SetUpBroyden"
    HOneBrd = HOne 
    RefBrd  = Ref 
    RootaBrd = Roota
    RootbBrd = Rootb
    RootyBrd = Rooty
    WeightspBrd = Weightsp
	WeightpgBrd = Weightpg
    fspBrd = fsp
	fpgBrd = fpg
	fkBrd = fk
	R1Brd = R1
	R2Brd = R2
	RpgBrd = Rpg
	RkBrd = Rk
    commBrd = comm
    Return
    End Subroutine SetUpBroyden 

    Subroutine ShutDownBroyden
    Implicit None
    Integer       :: IAlloc
    Deallocate(HOneBrd, RefBrd, RootaBrd, RootbBrd, RootyBrd, WeightspBrd, MixVec,&
				R1Brd, R2Brd, RpgBrd, RkBrd, Stat=IAlloc)
    If(IAlloc /= 0) Stop "Could not deallocate in ShutDownBroyden"
    Return
    End Subroutine ShutDownBroyden

    Subroutine Mat2Vec(y,Ene,T1,T2)
    Implicit None
    Complex (Kind=pr),   Intent(In)  :: Ene, T1(NOccBrd+1:NSOBrd,NOccBrd)
    Complex (Kind=pr),   Intent(In)  :: T2(NOccBrd+1:NSOBrd,NOccBrd+1:NSOBrd,NOccBrd,NOccBrd)
    Complex (Kind=pr),   Intent(Out) :: y(nBrd)
    Integer       :: IamNum
    Integer       :: I, J, A, B
    IamNum = 0
! Put Ene into y
    IamNum = IamNum + 1
    y(IamNum) = Ene 
! Put T1 into y
    Do I = 1, NOccBrd
    Do A = NOccBrd+1, NSOBrd
      IamNum    = IamNum + 1
      y(IamNum) = T1(A,I)
    EndDo
    EndDo
! Put T2tilde into y
    Do J = 2, NOccBrd
    Do I = 1, J-1 
    Do B = NOccBrd+2, NSOBrd
    Do A = NOccBrd+1, B-1 
      IamNum    = IamNum + 1
      y(IamNum) = T2(A,B,I,J)
    EndDo
    EndDo
    EndDo
    EndDo
! Check compatibility between dimension of y and Utildes
    If(nBrd.ne.IamNum) Stop "Dimension of y is not compatible within Broyden"
    Return
    End Subroutine Mat2Vec

    Subroutine Vec2Mat(y,Ene,T1,T2)
    Implicit None
    Complex (Kind=pr),   Intent(In)  :: y(nBrd)
    Complex (Kind=pr),   Intent(Out) :: Ene, T1(NOccBrd+1:NSOBrd,NOccBrd)
    Complex (Kind=pr),   Intent(Out) :: T2(NOccBrd+1:NSOBrd,NOccBrd+1:NSOBrd,NOccBrd,NOccBrd)
    Complex (Kind=pr)    :: tmp 
    Integer              :: IamNum
    Integer              :: I, J, A, B
    IamNum = 0
! Get Ene from y
    IamNum = IamNum + 1
    Ene = y(IamNum) 
! Get T1 from y
    Do I = 1, NOccBrd
    Do A = NOccBrd+1, NSOBrd
      IamNum  = IamNum + 1
      T1(A,I) = y(IamNum) 
    EndDo
    EndDo
! Get T2 from y
    T2 = Zero
    Do J = 2, NOccBrd
    Do I = 1, J-1 
    Do B = NOccBrd+2, NSOBrd
    Do A = NOccBrd+1, B-1 
      IamNum = IamNum + 1
      tmp    = y(IamNum)
      T2(A,B,I,J) =  tmp
      T2(B,A,J,I) =  tmp
      T2(A,B,J,I) = -tmp
      T2(B,A,I,J) = -tmp
    EndDo
    EndDo
    EndDo
    EndDo
! Check compatibility between dimension of y and Utildes
    If(nBrd.ne.IamNum) Stop "Dimension of y is not compatible within Broyden "
    Return
    End Subroutine Vec2Mat


    Subroutine EvalF(x,fx,HTwo,NSO)
    Implicit None
    Complex (Kind=pr),   Intent(In)  :: x(nBrd)
    Integer,             Intent(In)  :: NSO
    Complex (Kind=pr),   Intent(In)  :: HTwo(NSO,NSO,NSO,NSO)
    Complex (Kind=pr),   Intent(Out) :: fx(nBrd)
    Complex (Kind=pr),   Allocatable :: T1(:,:), T2(:,:,:,:)
    Complex (Kind=pr),   Allocatable :: Res1(:,:), Res2(:,:,:,:)
    Complex (Kind=pr),   Allocatable :: OlapEx1(:,:), OlapEx2(:,:,:,:)
    Complex (Kind=pr),   Allocatable :: EOlapEx1(:,:), EOlapEx2(:,:,:,:)
    Complex (Kind=pr)                :: Olap0, EOlap0
    Complex (Kind=pr)                :: Ene, Res0
    Complex (Kind=pr)                :: Norm1, Norm2, scal
    Integer                          :: IAlloc
! Allocate
    Allocate(T1(NOccBrd+1:NSOBrd,NOccBrd), Res1(NOccBrd+1:NSOBrd,NOccBrd), &
             T2(NOccBrd+1:NSOBrd,NOccBrd+1:NSOBrd,NOccBrd,NOccBrd),        &
             Res2(NOccBrd+1:NSOBrd,NOccBrd+1:NSOBrd,NOccBrd,NOccBrd),      &
             OlapEx1(NOccBrd+1:NSOBrd,NOccBrd),     & 
             EOlapEx1(NOccBrd+1:NSOBrd,NOccBrd),    &
             OlapEx2(NOccBrd+1:NSOBrd,NOccBrd+1:NSOBrd,NOccBrd,NOccBrd),   &
             EOlapEx2(NOccBrd+1:NSOBrd,NOccBrd+1:NSOBrd,NOccBrd,NOccBrd),  &
             Stat=IAlloc)
    If(IAlloc /= 0) Stop "Could not allocate in EvalF"
! Determine what T1 and T2 are
    Call Vec2Mat(x,Ene,T1,T2)
! Build Residuals 
    If(DoCCD) T1 = Zero ! CCD!
    Call PbarHbarOlap(Olap0,OlapEx1,OlapEx2,EOlap0,EOlapEx1,        &
                      EOlapEx2,HOneBrd,HTwo,RefBrd,T1,T2,        &
                      NOccBrd,NAOBrd,NSOBrd,ENucBrd,QJBrd,SPBrd,  &
					  NPointsBrd, NpgBrd, ncispBrd, ncipgBrd, ncikBrd, &
					  R1Brd, R2Brd, RpgBrd,RkBrd,RootaBrd,RootbBrd,RootyBrd, &
					  WeightspBrd,WeightpgBrd,fspBrd,fpgBrd,fkBrd,TrunODEBrd,commBrd)
!    Call PbarHbarResFromOlap2(Res0,Res1,Res2,RefBrd,Ene,Olap0,      &
!                      OlapEx1,OlapEx2,EOlap0,EOlapEx1,EOlapEx2,     &
!                      NOccBrd,NAOBrd,NSOBrd)
    Call PbarHbarResFromOlap(Ene,Res0,Res1,Res2,  &
	     RefBrd,Olap0,OlapEx1,OlapEx2,EOlap0,&
		 EOlapEx1,EOlapEx2,NOccBrd,NAOBrd,NSOBrd,SPBrd)
    EneBrd = Ene
    If(.not.DoCCD) then
!     Call ProjZeroMode(Res0,Res1,Res2,NOccBrd,NSOBrd) 
    EndIf
!   Call CheckInvSyRes(Res0,Res1,Res2,RefBrd,NOccBrd,NAOBrd,NSOBrd)
    If(DoCCD) Res1 = Zero ! CCD!
! Put dUtilde into dy
    Call Mat2Vec(fx,Res0,Res1,Res2)
!    fx = Cmplx(Real(fx,kind=pr),kind=pr)
! Deallocate
    Deallocate(T1, Res1, T2, Res2, OlapEx1, EOlapEx1, OlapEx2, EOlapEx2, &
               Stat=IAlloc)
    If(IAlloc /= 0) Stop "Could not deallocate in EvalF"
    Return
    End Subroutine EvalF

    Subroutine BroydenStep(xnew,xold,fold,mix,w0,w,dF,dx,NIter)
    Implicit None
    Integer,           Intent(In)  :: NIter
    Complex (Kind=pr), Intent(In)  :: xold(nBrd), fold(nBrd) 
    Complex (Kind=pr), Intent(In)  :: dF(nBrd,pBrd), dx(nBrd,pBrd) 
    Real (Kind=pr),    Intent(In)  :: w(pBrd) 
    Real (Kind=pr),    Intent(In)  :: mix, w0
    Complex (Kind=pr), Intent(Out) :: xnew(nBrd) 
    Complex (Kind=pr), Allocatable :: AMat(:,:)
    Complex (Kind=pr), Allocatable :: BetaMat(:,:), GammaVec(:)
    Complex (Kind=pr), Allocatable :: UMat(:,:), CVec(:)
    Integer         :: I, J 
    Integer         :: IAlloc
    xnew = xold 
    If (pBrd>0) then
! Allocate
    Allocate(AMat(pBrd,pBrd), UMat(nBrd,pBrd), CVec(pBrd),  &
             BetaMat(pBrd,pBrd), GammaVec(pBrd),            &
             Stat=IAlloc)
    If(IAlloc/=0) Stop "Could not allocate in BroydenStep"
! Build a matrix
! a_kn = wk wn <dF_n|dF_k>
    Do I = 1, pBrd
    Do J = 1, pBrd
      AMat(I,J) = w(I)*w(J)*Dot_Product(dF(:,J),dF(:,I)) 
    EndDo
    EndDo
! beta = (w0^2*I + a)^-1
    Do I = 1, pBrd
      AMat(I,I) = AMat(I,I) + w0*w0
    EndDo
    Call InvertC(AMat,BetaMat,pBrd)
! c_k = wk <dF_k|Fold>
    Do I = 1, pBrd
      CVec(I) = w(I)*Dot_Product(dF(:,I),fold)
    EndDo
! gamma_n = c_k beta_kn
    GammaVec = Zero
    Do J = 1, pBrd
    Do I = 1, pBrd
      GammaVec(J) = GammaVec(J) + CVec(I)*BetaMat(I,J) 
    EndDo
    EndDo
! u_n = mix dF_n + dx_n
    Do I = 1, pBrd
    Do J = 1, nBrd
      UMat(J,I) = mix*MixVec(J)*dF(J,I) + dx(J,I)
    EndDo
    EndDo
! x_new = x_old + mix F_old - w_n u_n gamma_n
!   Call outmat(6,pBrd,1,2,GammaVec,"Gamma Vector in Broyden")
! For the first few steps, only do mixing
    Do I = 1, pBrd
      xnew = xnew - w(I)*GammaVec(I)*UMat(:,I)
    EndDo
! Deallocate
    Deallocate(AMat, UMat, CVec, BetaMat, GammaVec, &
               Stat=IAlloc)
    If(IAlloc/=0) Stop "Could not deallocate in BroydenStep"
    EndIf
! Simple mixing if pBrd==0
    Do J = 1, nBrd
      xnew(J) = xnew(J) + mix * MixVec(J) * fold(J)
    EndDo
    Return
    End Subroutine BroydenStep

    Subroutine UpdateBroyden(dF,dx,fnew,fold,xnew,xold)
    Implicit None
    Complex(Kind=pr), Intent(In)    :: fnew(nBrd), fold(nBrd)
    Complex(Kind=pr), Intent(In)    :: xnew(nBrd), xold(nBrd)
    Complex(Kind=pr), Intent(InOut) :: dF(nBrd,pBrd), dx(nBrd,pBrd) 
    Real(Kind=pr)                   :: norm
    Integer                         :: I
    If (pBrd>0) then
! Put all dF and dx vec one index to the left
    Do I = 1, pBrd-1
      dF(:,I) = dF(:,I+1)
      dx(:,I) = dx(:,I+1)
    EndDo
! Build norm = |fnew-fold| 
    norm = sqrt(Dot_Product(fnew-fold,fnew-fold))
	if (norm < 1e-10) Stop "Norm too small"
! dx = (xnew - xold) / norm 
    dx(:,pBrd) = (xnew-xold)/norm
! df = (fnew - fold) / norm 
    dF(:,pBrd) = (fnew-fold)/norm
    EndIf
    Return
    End Subroutine UpdateBroyden

	Subroutine ProjectNull(NullVec,Vec,M)
	Implicit None
	Integer,          Intent(In)    :: M
	Complex(Kind=pr), Intent(In)    :: NullVec(nBrd,M)
	Complex(Kind=pr), Intent(InOut) :: Vec(nBrd)
	Complex(Kind=pr) :: tmp(M), tmp1
	tmp = matmul(transpose(conjg(NullVec)), Vec)
    tmp1 = Vec(1)
	Vec = Vec - matmul(NullVec, tmp)
    Vec(1) = tmp1
	End Subroutine ProjectNull

	Subroutine GuessVec(T1,T2,NullVec,NOcc,NSO)
	Implicit None
	Integer,          Intent(In)    :: NSO,NOcc
	Complex(Kind=pr), Intent(In)    :: T1(NOcc+1:NSO,NOcc)
	Complex(Kind=pr), Intent(In)    :: T2(NOcc+1:NSO,NOcc+1:NSO,NOcc,NOcc)
	Complex(Kind=pr), Intent(Out)   :: NullVec(nBrd,nBrd)
	Real(Kind=pr) :: MaxT1, MaxT2, MaxReal, MaxImag
	Integer :: A,B,I,J,M,N
	MaxT1 = 1e-3*maxval(abs(T1))
	MaxT2 = 1e-3*maxval(abs(T2))
	NullVec = Zero
	M = 1
	N = 0
	Do I = 1, NOccBrd
	Do A = NOccBrd+1, NSOBrd
		M = M + 1
		If (abs(T1(A,I)) > MaxT1) then
			N = N + 1
			NullVec(M,N) = One
		End If
	End Do
	End Do
	Do J = 2, NOccBrd
	Do I = 1, J-1
	Do B = NOccBrd+2, NSOBrd
	Do A = NOccBrd+1, B-1
		M = M + 1
		If (abs(T2(A,B,I,J)) > MaxT2) then
			N = N + 1
			NullVec(M,N) = One
		End If
	End Do
	End Do
	End Do
	End Do
	End Subroutine GuessVec

	Subroutine BuildSRes(Res)
	Implicit None
	Complex (Kind=pr), Intent(InOut)  :: Res(nBrd)
	Complex (Kind=pr)              :: r0tmp, r1tmp(NOccBrd+1:NSOBrd,NOccBrd)
	Complex (Kind=pr)              :: r2tmp(NOccBrd+1:NSOBrd,NOccBrd+1:NSOBrd,NOccBrd,NOccBrd)
    Complex (Kind=pr)              :: t1(NOccBrd+1:NSOBrd,NOccBrd)
    Complex (Kind=pr)              :: t2(NOccBrd+1:NSOBrd,NOccBrd+1:NSOBrd,NOccBrd,NOccBrd)
	Complex (Kind=pr)              :: tmp
    Integer :: i,j,a,b
	tmp = Res(1)
	Call Vec2Mat(Res,r0tmp,r1tmp,r2tmp)
    Call BuildSC_test(t1,t2,r1tmp,r2tmp,NOccBrd,NSOBrd,NPointsBrd,NpgBrd,QJBrd,SPBrd, &
                      ncispBrd,ncipgBrd,ncikBrd,R1Brd, R2Brd, RpgBrd, RkBrd, &
                      RootaBrd,RootbBrd,RootyBrd,WeightspBrd,WeightpgBrd,fspBrd,fpgBrd,fkBrd,commBrd)
	Call Mat2Vec(Res,tmp,t1,t2)
	Res = Res / damp
	End Subroutine BuildSRes



    End Module Broyden
