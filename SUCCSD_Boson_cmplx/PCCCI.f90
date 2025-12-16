
    Module PCCCI
    Use Precision
    Use Constants
    Use CCRes
    Use IntTrans
    Use UHF

! Trans CISD to CCSD
! C0 + C1 + C2 = exp(T0+T1+T2) = exp(T0)(1+T1+T2+T1^2/2)

    Contains

    Subroutine BuildW(R,T1,T2,W0,W1,W2,NOcc,NSO)
    Implicit None
    Integer         , Intent(In)  :: NOcc, NSO
    Complex(kind=pr), Intent(In)  :: R(NSO,NSO)
    Complex(kind=pr), Intent(In)  :: T1(NOcc+1:NSO,NOcc)
    Complex(kind=pr), Intent(In)  :: T2(NOcc+1:NSO,NOcc+1:NSO,NOcc,NOcc)
    Complex(kind=pr), Intent(Out) :: W0, W1(NOcc+1:NSO,NOcc)
    Complex(kind=pr), Intent(Out) :: W2(NOcc+1:NSO,NOcc+1:NSO,NOcc,NOcc)
    Complex(kind=pr) :: C0, C1(NOcc+1:NSO,NOcc)
    Complex(kind=pr) :: C2(NOcc+1:NSO,NOcc+1:NSO,NOcc,NOcc)
    Complex(kind=pr) :: a0, a1(NSO-NOcc,NOcc)
    Complex(kind=pr) :: a2(NSO-NOcc,NSO-NOcc,NOcc,NOcc)
    Integer :: I
! Drudge ST T2
!    Call TransQuadT2(Z,T1,T2,C0,C1,C2,NOcc,NSO-NOcc)
! new idea
    Call TransT2(R,T1,T2,a0,a1,a2,NOcc,NSO)
!    Call ASymm2(a2,NOcc,NSO)
!    a2 = a2 / 4
!    if (abs(c0-a0) > 1e-6) Print *, "c0 wrong"
!    if (maxval(abs(c1(NOcc+1:NSO,:NOcc)-a1)) > 1e-6) Print *, "c1 wrong"
!    if (maxval(abs(c2(NOcc+1:NSO,NOcc+1:NSO,:NOcc,:NOcc)-a2)) > 1e-6) Print *, "c2 wrong"
    C0 = a0
    C1(NOcc+1:NSO,:NOcc) = a1
    C2(NOcc+1:NSO,NOcc+1:NSO,:NOcc,:NOcc) = a2
    Call CI2CC(W0,W1,W2,C0,C1,C2,NOcc,NSO)
    End Subroutine BuildW

    Subroutine CC2CI(T0,T1,T2,C0,C1,C2,NOcc,NSO)
    Implicit None
    Integer         , Intent(In)  :: NOcc, NSO
    Complex(kind=pr), Intent(In)  :: T0, T1(NOcc+1:NSO,NOcc)
    Complex(kind=pr), Intent(In)  :: T2(NOcc+1:NSO,NOcc+1:NSO,NOcc,NOcc)
    Complex(kind=pr), Intent(Out) :: C0, C1(NOcc+1:NSO,NOcc)
    Complex(kind=pr), Intent(Out) :: C2(NOcc+1:NSO,NOcc+1:NSO,NOcc,NOcc)
    Integer :: a,b,c,d,i,j,k,l
    C1 = Zero
    C2 = Zero
    C0 = Exp(T0)
    C1 = C0 * T1
    C2 = T2 / 4
    !$omp parallel do
    Do a = NOcc+1, NSO
    Do b = NOcc+1, NSO
    Do i = 1, NOcc
    Do j = 1, NOcc
        if (a==b .or. i==j) cycle
        C2(a,b,i,j) = C2(a,b,i,j) + T1(a,i) * T1(b,j) / 2
    End Do
    End Do
    End Do
    End Do
    !$omp end parallel do
    C2 = C0 * C2
    Call ASymm2(C2,NOcc,NSO)
    C2 = C2 / 4
    End Subroutine CC2CI

    Subroutine CI2CC(T0,T1,T2,C0,C1,C2,NOcc,NSO)
    Implicit None
    Integer         , Intent(In)  :: NOcc, NSO
    Complex(kind=pr), Intent(In)  :: C0, C1(NOcc+1:NSO,NOcc)
    Complex(kind=pr), Intent(In)  :: C2(NOcc+1:NSO,NOcc+1:NSO,NOcc,NOcc)
    Complex(kind=pr), Intent(Out) :: T0, T1(NOcc+1:NSO,NOcc)
    Complex(kind=pr), Intent(Out) :: T2(NOcc+1:NSO,NOcc+1:NSO,NOcc,NOcc)
    Integer :: a,b,c,d,i,j,k,l
    T1 = Zero
    T2 = Zero
    T0 = C0
    T1 = C1 / C0
    T2 = C2 / C0
    !$omp parallel do
    Do a = NOcc+1, NSO
    Do b = NOcc+1, NSO
    Do i = 1, NOcc
    Do j = 1, NOcc
        if (a==b .or. i==j) cycle
        T2(a,b,i,j) = T2(a,b,i,j) - T1(a,i) * T1(b,j) / 2
    End Do
    End Do
    End Do
    End Do
    !$omp end parallel do
    Call ASymm2(T2,NOcc,NSO)
    End Subroutine CI2CC
    
    Subroutine TransT2(R,T1,T2,C0,C1,C2,NOcc,NSO)
! double ST of 1 + T2 + 1/2 T2^2
! U T U^-1, U = e^-T1 R
    Implicit None
    Integer         , Intent(In)  :: NOcc, NSO
    Complex(kind=pr), Intent(In)  :: R(NSO,NSO), T1(NOcc+1:NSO,NOcc)
    Complex(kind=pr), Intent(In)  :: T2(NOcc+1:NSO,NOcc+1:NSO,NOcc,NOcc)
    Complex(kind=pr), Intent(Out) :: C0, C1(NSO-NOcc,NOcc)
    Complex(kind=pr), Intent(Out) :: C2(NSO-NOcc,NSO-NOcc,NOcc,NOcc)
    Complex(kind=pr) :: E0, F(NSO,NSO)
    Complex(kind=pr) :: U(NSO,NSO,NSO,NSO)
    Complex(kind=pr) :: V(NSO,NSO)
    Complex(kind=pr), allocatable :: fov(:,:), fvo(:,:)
    Complex(kind=pr), allocatable :: foo(:,:), fvv(:,:)
    Complex(kind=pr), allocatable :: vvoo(:,:,:,:), oovv(:,:,:,:)
    Complex(kind=pr), allocatable :: ooov(:,:,:,:), ovvv(:,:,:,:)
    Complex(kind=pr), allocatable :: oooo(:,:,:,:), vvvv(:,:,:,:)
    Complex(kind=pr), allocatable :: ovoo(:,:,:,:), ovov(:,:,:,:)
    Complex(kind=pr), allocatable :: vvov(:,:,:,:), C2tmp(:,:,:,:)
    Integer :: NVrt, nv2,no2,ov, a,b,c,d,i,j,k,l
    Integer :: nvsize = 2, nosize = 4, m, m1, nblkv, nblko
    integer, allocatable :: loadv(:), isv(:), iev(:)
    Integer         :: cnt1, cnt2, clock_rate, clock_max
    NVrt = NSO - NOcc
    nv2 = NVrt**2
    no2 = NOcc**2
    ov = NVrt*NOcc
    nblkv = NVrt / nvsize
    Allocate(loadv(nblkv), isv(nblkv), iev(nblkv))
    loadv = NVrt / nblkv
    loadv(1:mod(nvrt,nblkv)) = loadv(1:mod(nvrt,nblkv)) + 1
    do m = 1, nblkv
        isv(m) = sum(loadv(1:m-1))
        iev(m) = isv(m) + loadv(m)
    end do
    C0 = One
    C1 = Zero
    C2 = Zero
    F = Zero
    ! build Transformed T2 and normal order
    Call IDMat(V, NSO)
    V(NOcc+1:NSO,:NOcc) = -T1
	V = matmul(V,R)
    U = Zero
    U(NOcc+1:NSO,NOcc+1:NSO,:NOcc,:NOcc) = T2
    Call IntTran4Sim(U, V, NSO)
    Call MKFock_MO(Zero*V,U,F,E0,NOcc,NSO)
    ! 1st order
    C0 = C0 + E0
    C1 = C1 + F(NOcc+1:NSO,:NOcc)
    C2 = C2 + U(NOcc+1:NSO,NOcc+1:NSO,:NOcc,:NOcc) / 4
    ! 2nd order
    allocate(fov(NOcc,NVrt), fvo(NVrt,NOcc), fvv(NVrt,NVrt), foo(NOcc,NOcc))
    fov = F(:NOcc,NOcc+1:NSO)
    fvo = F(NOcc+1:NSO,:NOcc)
    foo = F(:NOcc,:NOcc)
    fvv = F(NOcc+1:NSO,NOcc+1:NSO)
    allocate(vvoo(NVrt,NVrt,NOcc,NOcc), oovv(NOcc,NOcc,NVrt,NVrt))
    vvoo = U(NOcc+1:NSO,NOcc+1:NSO,:NOcc,:NOcc)
    oovv = U(:NOcc,:NOcc,NOcc+1:NSO,NOcc+1:NSO)
    ! c0,c1
    C0 = C0 + E0**2 / 2 + sum(fov * transpose(fvo)) / 2
    Call tA_dot_tB(oovv/8,vvoo,C0,no2,nv2)
    C1 = C1 + E0 * fvo
    C1 = C1 - matmul(fvo,foo)/2 + matmul(fvv,fvo)/2
    allocate(ooov(NOcc,NOcc,NOcc,NVrt))
    ooov = U(:NOcc,:NOcc,:NOcc,NOcc+1:NSO)
    Call Transpose2(vvoo,oovv,nv2,no2)
    !$omp parallel do reduction(+:c1)
    Do m = 1, nblkv
    Do b = isv(m)+ 1, iev(m)
!        Call A_dot_B(-vvoo(:,b,:,:)/4,ooov(:,:,:,b),C1,NVrt,no2,NOcc)
        Call At_dot_B(-oovv(:,:,:,b)/4,ooov(:,:,:,b),C1,no2,NVrt,NOcc)
    End Do
    End Do
    !$omp end parallel do
    deallocate(ooov)
    allocate(vvov(NVrt,NVrt,NOcc,NVrt))
    Call Transpose2(U(:NOcc,NOcc+1:NSO,NOcc+1:NSO,NOcc+1:NSO),vvov,ov,nv2)
    !$omp parallel do reduction(+:c1)
    Do j = 1, NOcc
        Call At_dot_B(-vvov(:,:,j,:)/4,vvoo(:,:,:,j),C1,nv2,NVrt,NOcc)
    End Do
    !$omp end parallel do
    deallocate(vvov)
    allocate(ovov(NOcc,NVrt,NOcc,NVrt))
    ovov = U(:NOcc,NOcc+1:NSO,:NOcc,NOcc+1:NSO)
    !$omp parallel do
    Do a = 1, NVrt
    Do i = 1, NOcc
        C1(a,i) = C1(a,i) + sum(vvoo(a,:,i,:) * transpose(fov)) / 2
        C1(a,i) = C1(a,i) - sum(ovov(:,a,i,:) * transpose(fvo)) / 2
    End Do
    End Do
    !$omp end parallel do
    ! c2
    C2 = C2 + E0 / 4 * vvoo
    !$omp parallel do
    Do a = 1, NVrt
    Do b = 1, NVrt
    Do i = 1, NOcc
    Do j = 1, NOcc
        C2(a,b,i,j) = C2(a,b,i,j) - fvo(a,j) * fvo(b,i) / 2
        C2(a,b,i,j) = C2(a,b,i,j) + sum(oovv(j,:,a,:) * ovov(:,b,i,:)) / 2
    End Do
    End Do
    End Do
    End Do
    !$omp end parallel do
    deallocate(ovov)
    allocate(C2tmp(NOcc,NOcc,NVrt,NVrt))
    C2tmp = Zero
    allocate(vvvv(NVrt,NVrt,NVrt,NVrt))
    vvvv = U(NOcc+1:NSO,NOcc+1:NSO,NOcc+1:NSO,NOcc+1:NSO)
    Call Transpose1(vvvv,nv2,nv2)
    !$omp parallel do
    Do m = 1, nblkv
    Do b = isv(m)+ 1, iev(m)
        Call At_dot_B(vvoo/16,vvvv(:,:,:,b),C2tmp(:,:,:,b),nv2,no2,nvrt)
    End Do
    End Do
    !$omp end parallel do
    deallocate(vvvv)
    allocate(oooo(NOcc,NOcc,NOcc,NOcc))
    oooo = U(:NOcc,:NOcc,:NOcc,:NOcc)
    !$omp parallel do
    Do m = 1, nblkv
    Do b = isv(m)+ 1, iev(m)
        Call At_dot_B(oooo/16,oovv(:,:,:,b),C2tmp(:,:,:,b),no2,no2,nvrt)
    End Do
    End Do
    !$omp end parallel do
    deallocate(oooo)
    !$omp parallel do
    Do m = 1, nblkv
    Do b = isv(m)+ 1, iev(m)
    Do a = 1, NVrt
        Call A_dot_B(oovv(:,:,a,b),-foo/4,C2tmp(:,:,a,b),nocc,nocc,NOcc)
        !C2tmp(:,:,a,b) = C2tmp(:,:,a,b) - matmul(oovv(:,:,a,b),foo) / 2
    End Do
    End Do
    End Do
    !$omp end parallel do
    !$omp parallel do
    Do m = 1, nblkv
    Do b = isv(m)+ 1, iev(m)
        Call A_dot_Bt(oovv(:,:,:,b),fvv/4,C2tmp(:,:,:,b),no2,NVrt,NVrt)
    End Do
    End Do
    !$omp end parallel do
    deallocate(oovv)
    allocate(ooov(NOcc,NOcc,NOcc,NVrt))
    Call Transpose2(U(:NOcc,NOcc+1:NSO,:NOcc,:NOcc),ooov,ov,no2)
    !$omp parallel do
    Do m = 1, nblkv
    Do b = isv(m)+ 1, iev(m)
        Call A_dot_Bt(ooov(:,:,:,b),-fvo/4,C2tmp(:,:,:,b),no2,NOcc,NVrt)
    End Do
    End Do
    !$omp end parallel do
    deallocate(ooov)
    allocate(ovvv(NOcc,NVrt,NVrt,NVrt))
    Call Transpose2(U(NOcc+1:NSO,NOcc+1:NSO,:NOcc,NOcc+1:NSO),ovvv,nv2,ov)
    !$omp parallel do
    Do m = 1, nblkv
    Do b = isv(m)+ 1, iev(m)
    Do a = 1, NVrt
        !Call A_dot_B(ovvv(:,:,a,b),fvo/2,C2tmp(:,:,a,b),NOcc,NVrt,NOcc)
        C2tmp(:,:,a,b) = C2tmp(:,:,a,b) + matmul(ovvv(:,:,a,b),fvo) / 4
    End Do
    End Do
    End Do
    !$omp end parallel do
    deallocate(ovvv)
    Call Transpose2(C2tmp,vvoo,no2,nv2)
    C2 = C2 + vvoo
    End Subroutine TransT2

    Subroutine tA_dot_tB(A,B,C,MA,NA)
    Implicit None
    Integer         , Intent(In)     :: MA, NA
    Complex (Kind=pr), Intent(in)    :: A(MA,NA), B(NA,MA)
    Complex (Kind=pr), Intent(InOut) :: C
    C = C + sum(A * transpose(B))
    End Subroutine

    End Module

