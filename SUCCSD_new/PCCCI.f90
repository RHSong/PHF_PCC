
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
!=====================================================================!
! double ST of 1 + T2 + 1/2 T2^2 :  U T U^-1, U = e^-T1 R                !
! Same algebra as the original routine, but the transformed tensor      !
! U' = V U V^-1 (U nonzero only in the vvoo block) is never built as an !
! NSO^4 array: every block that the second-order terms need is formed   !
! directly with TransBlockT2 (largest object o v^3), and the vvvv block !
! is never formed at all - its contraction with vvoo' is done           !
! "contract first" (O(o^2 v^3) instead of O(o^2 v^4 + v^5)).            !
!=====================================================================!
    Implicit None
    Integer         , Intent(In)  :: NOcc, NSO
    Complex(kind=pr), Intent(In)  :: R(NSO,NSO), T1(NOcc+1:NSO,NOcc)
    Complex(kind=pr), Intent(In)  :: T2(NOcc+1:NSO,NOcc+1:NSO,NOcc,NOcc)
    Complex(kind=pr), Intent(Out) :: C0, C1(NSO-NOcc,NOcc)
    Complex(kind=pr), Intent(Out) :: C2(NSO-NOcc,NSO-NOcc,NOcc,NOcc)
    Complex(kind=pr), Allocatable :: VV(:,:), Vi(:,:), Lo(:,:), Lvv(:,:), Moo(:,:), Mov(:,:), Ido(:,:)
    Complex(kind=pr), Allocatable :: Nt(:,:), G(:,:), F(:,:), Lv(:,:), Mo(:,:)
    Complex(kind=pr), Allocatable :: fov(:,:), fvo(:,:), foo(:,:), fvv(:,:)
    Complex(kind=pr), Allocatable :: vvoo(:,:,:,:), oovv(:,:,:,:), ooov(:,:,:,:), ovov(:,:,:,:)
    Complex(kind=pr), Allocatable :: oooo(:,:,:,:), ovoo(:,:,:,:), vvov(:,:,:,:), ovvv(:,:,:,:)
    Complex(kind=pr), Allocatable :: X(:,:,:,:), Y(:,:,:,:), Z(:,:,:,:), W(:,:,:,:)
    Complex(kind=pr) :: E0, c0t
    Integer :: o, v, nv2, no2, ov, a, b, c, i, j, k, IAlloc
    o = NOcc
    v = NSO - NOcc
    nv2 = v*v; no2 = o*o; ov = o*v
!---- V = (1 - T1) R  and its inverse; the four sub-blocks that act on the vvoo block
    Allocate(VV(NSO,NSO), Vi(NSO,NSO), Lo(o,v), Lvv(v,v), Moo(o,o), Mov(o,v), Ido(o,o), &
             Lv(NSO,v), Mo(o,NSO), Nt(o,v), G(v,o), F(NSO,NSO), &
             fov(o,v), fvo(v,o), foo(o,o), fvv(v,v), Stat=IAlloc)
    If (IAlloc /= 0) Stop "Could not allocate in TransT2"
    Call IDMat(VV, NSO)
    VV(o+1:NSO,1:o) = -T1
    VV = matmul(VV,R)
    Call InvertC(VV,Vi,NSO)
    Lv  = VV(:,o+1:NSO)          ! bra indices:  U'(p,...) = sum_a V(p,a) ...
    Lo  = VV(1:o,o+1:NSO)
    Lvv = VV(o+1:NSO,o+1:NSO)
    Mo  = Vi(1:o,:)             ! ket indices:  ... sum_i U(..,i) Vi(i,r)
    Moo = Vi(1:o,1:o)
    Mov = Vi(1:o,o+1:NSO)
    Call IDMat(Ido, o)
!---- E0 and the one-body part F(p,q) = sum_i U'(p,i,q,i),  E0 = 1/2 sum_ij U'(i,j,i,j)
    Nt = matmul(Moo,Lo)         ! Nt(j,b) = sum_i Vi(j,i) V(i,b)
    E0 = Zero
    G = Zero
    !$omp parallel do private(a,b,i,j) reduction(+:E0)
    Do i = 1, o
    Do a = 1, v
      Do j = 1, o
      Do b = 1, v
        G(a,i) = G(a,i) + T2(o+a,o+b,i,j)*Nt(j,b)
      End Do
      End Do
      E0 = E0 + F12*G(a,i)*Nt(i,a)
    End Do
    End Do
    !$omp end parallel do
    F = matmul(Lv, matmul(G, Mo))
    fov = F(1:o,o+1:NSO)
    fvo = F(o+1:NSO,1:o)
    foo = F(1:o,1:o)
    fvv = F(o+1:NSO,o+1:NSO)
!---- blocks of U'
    Allocate(vvoo(v,v,o,o), oovv(o,o,v,v), ooov(o,o,o,v), ovov(o,v,o,v), oooo(o,o,o,o), &
             ovoo(o,v,o,o), Stat=IAlloc)
    If (IAlloc /= 0) Stop "Could not allocate blocks in TransT2"
    Call TransBlockT2(T2,o,v,Lvv,v,Lvv,v,Moo,o,Moo,o,vvoo)
    Call TransBlockT2(T2,o,v,Lo ,o,Lo ,o,Mov,v,Mov,v,oovv)
    Call TransBlockT2(T2,o,v,Lo ,o,Lo ,o,Moo,o,Mov,v,ooov)
    Call TransBlockT2(T2,o,v,Lo ,o,Lvv,v,Moo,o,Mov,v,ovov)
    Call TransBlockT2(T2,o,v,Lo ,o,Lo ,o,Moo,o,Moo,o,oooo)
    Call TransBlockT2(T2,o,v,Lo ,o,Lvv,v,Moo,o,Moo,o,ovoo)
!======================================================================!
!  first order
    C0 = One + E0
    C1 = fvo
    C2 = F14*vvoo
!  C0, second order
    C0 = C0 + F12*E0*E0 + F12*sum(fov*transpose(fvo))
    c0t = Zero
    !$omp parallel do private(a,b,i,j) reduction(+:c0t)
    Do j = 1, o
    Do i = 1, o
    Do b = 1, v
    Do a = 1, v
      c0t = c0t + oovv(i,j,a,b)*vvoo(a,b,i,j)
    End Do
    End Do
    End Do
    End Do
    !$omp end parallel do
    C0 = C0 + c0t/8
!  C1, second order
    C1 = C1 + E0*fvo - F12*matmul(fvo,foo) + F12*matmul(fvv,fvo)
    ! - 1/4 sum_{k,l,b} vvoo(a,b,k,l) ooov(k,l,i,b)
    Do b = 1, v
      Call ZGEMM('N','N',v,o,no2,-F14*ZOne,vvoo(1,b,1,1),nv2,ooov(1,1,1,b),no2,ZOne,C1,v)
    End Do
    ! - 1/4 sum_{c,d,j} ovvv(j,a,c,d) vvoo(c,d,i,j)  with X(c,d,a,j) = ovvv(j,a,c,d)
    Allocate(ovvv(o,v,v,v), X(v,v,v,o), Stat=IAlloc)
    If (IAlloc /= 0) Stop "Could not allocate ovvv in TransT2"
    Call TransBlockT2(T2,o,v,Lo ,o,Lvv,v,Mov,v,Mov,v,ovvv)
    !$omp parallel do private(a,c,i,j)
    Do j = 1, o
    Do a = 1, v
    Do c = 1, v
      X(c,:,a,j) = ovvv(j,a,c,:)
    End Do
    End Do
    End Do
    !$omp end parallel do
    Deallocate(ovvv)
    Do j = 1, o
      Call ZGEMM('T','N',v,o,nv2,-F14*ZOne,X(1,1,1,j),nv2,vvoo(1,1,1,j),nv2,ZOne,C1,v)
    End Do
    Deallocate(X)
    ! + 1/2 sum_{b,j} vvoo(a,b,i,j) fov(j,b) - 1/2 sum_{k,b} ovov(k,a,i,b) fvo(b,k)
    !$omp parallel do private(a,i,b,j,k)
    Do i = 1, o
    Do a = 1, v
      Do j = 1, o
      Do b = 1, v
        C1(a,i) = C1(a,i) + F12*vvoo(a,b,i,j)*fov(j,b) - F12*ovov(j,a,i,b)*fvo(b,j)
      End Do
      End Do
    End Do
    End Do
    !$omp end parallel do
!  C2, second order (all terms written directly in (a,b,i,j) order)
    !$omp parallel do private(a,b,i,j)
    Do j = 1, o
    Do i = 1, o
    Do b = 1, v
    Do a = 1, v
      C2(a,b,i,j) = C2(a,b,i,j) + F14*E0*vvoo(a,b,i,j) - F12*fvo(a,j)*fvo(b,i)
    End Do
    End Do
    End Do
    End Do
    !$omp end parallel do
    ! + 1/2 sum_{k,c} vvoo(a,c,j,k) ovov(k,b,i,c) :  X(a,j,c,k) = vvoo(a,c,j,k), Y(c,k,b,i) = ovov(k,b,i,c)
    Allocate(X(v,o,v,o), Y(v,o,v,o), Z(v,o,v,o))
    !$omp parallel do private(a,c,j,k)
    Do k = 1, o
    Do c = 1, v
    Do j = 1, o
    Do a = 1, v
      X(a,j,c,k) = vvoo(a,c,j,k)
    End Do
    End Do
    End Do
    End Do
    !$omp end parallel do
    !$omp parallel do private(b,c,i,k)
    Do i = 1, o
    Do b = 1, v
    Do k = 1, o
    Do c = 1, v
      Y(c,k,b,i) = ovov(k,b,i,c)
    End Do
    End Do
    End Do
    End Do
    !$omp end parallel do
    Z = Zero
    Call PZGEMM('N','N',ov,ov,ov,ZOne,X,ov,Y,ov,ZOne,Z,ov)      ! Z(a,j,b,i)
    !$omp parallel do private(a,b,i,j)
    Do j = 1, o
    Do i = 1, o
    Do b = 1, v
    Do a = 1, v
      C2(a,b,i,j) = C2(a,b,i,j) + F12*Z(a,j,b,i)
    End Do
    End Do
    End Do
    End Do
    !$omp end parallel do
    Deallocate(X, Y, Z)
    ! + 1/16 sum_{cd} vvoo(c,d,i,j) vvvv(a,b,c,d)   (contract first, no vvvv block)
    Allocate(Y(o,o,o,o), Z(v,v,o,o), W(v,v,o,o))
    Call TransBlockT2(vvoo,o,v,Mov,o,Mov,o,Ido,o,Ido,o,Y)        ! Y(k,l,i,j) = sum_cd Mov(k,c) Mov(l,d) vvoo(c,d,i,j)
    Z = Zero
    Call PZGEMM('N','N',nv2,no2,no2,ZOne,T2,nv2,Y,no2,ZOne,Z,nv2)  ! Z(e,f,i,j) = sum_kl T2(e,f,k,l) Y(k,l,i,j)
    Call TransBlockT2(Z,o,v,Lvv,v,Lvv,v,Ido,o,Ido,o,W)           ! W(a,b,i,j) = sum_ef Lvv(a,e) Lvv(b,f) Z(e,f,i,j)
    C2 = C2 + W/16
    Deallocate(Y, Z, W)
    ! + 1/16 sum_{kl} oooo(k,l,i,j) vvoo(a,b,k,l)
    Call PZGEMM('N','N',nv2,no2,no2,ZOne/16,vvoo,nv2,oooo,no2,ZOne,C2,nv2)
    ! - 1/4 sum_k vvoo(a,b,i,k) foo(k,j)
    Call PZGEMM('N','N',nv2*o,o,o,-F14*ZOne,vvoo,nv2*o,foo,o,ZOne,C2,nv2*o)
    ! + 1/4 sum_c fvv(a,c) vvoo(c,b,i,j)
    Call PZGEMM('N','N',v,v*no2,v,F14*ZOne,fvv,v,vvoo,v,ZOne,C2,v)
    ! - 1/4 sum_k fvo(a,k) ovoo(k,b,i,j)
    Call PZGEMM('N','N',v,v*no2,o,-F14*ZOne,fvo,v,ovoo,o,ZOne,C2,v)
    ! + 1/4 sum_c vvov(a,b,i,c) fvo(c,j)
    Allocate(vvov(v,v,o,v), Stat=IAlloc)
    If (IAlloc /= 0) Stop "Could not allocate vvov in TransT2"
    Call TransBlockT2(T2,o,v,Lvv,v,Lvv,v,Moo,o,Mov,v,vvov)
    Call PZGEMM('N','N',nv2*o,o,v,F14*ZOne,vvov,nv2*o,fvo,v,ZOne,C2,nv2*o)
    Deallocate(vvov)
    Deallocate(VV, Vi, Lo, Lvv, Moo, Mov, Ido, Lv, Mo, Nt, G, F, fov, fvo, foo, fvv, &
               vvoo, oovv, ooov, ovov, oooo, ovoo)
    End Subroutine TransT2

    Subroutine tA_dot_tB(A,B,C,MA,NA)
    Implicit None
    Integer         , Intent(In)     :: MA, NA
    Complex (Kind=pr), Intent(in)    :: A(MA,NA), B(NA,MA)
    Complex (Kind=pr), Intent(InOut) :: C
    C = C + sum(A * transpose(B))
    End Subroutine

    End Module

