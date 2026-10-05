    Module CCRes
    Use Precision
    Use Constants
    Use ERIBlocks
	Implicit None
	Integer :: blksize = 2

    Contains

      Subroutine CCEnergy(T1,T2,Fock,ERI,ECorr,NOcc,NSO)
      Implicit None
      Integer,           Intent(In)  :: NOcc,NSO
      Complex (Kind=pr), Intent(In)  :: Fock(NSO,NSO)
      Complex (Kind=pr), Intent(In)  :: ERI(NSO,NSO,NSO,NSO)
      Complex (Kind=pr), Intent(In)  :: T1(NOcc+1:NSO,NOcc)
      Complex (Kind=pr), Intent(In)  :: T2(NOcc+1:NSO,NOcc+1:NSO,NOcc,NOcc)
      Complex (Kind=pr), Intent(Out) :: ECorr
      Integer :: I, J, A, B
      Complex (Kind=pr) :: E2
      ECorr = sum(Fock(1:NOcc,NOcc+1:NSO) * transpose(T1))
      E2 = Zero
      !$omp parallel do reduction(+:E2) private(I,J,A,B)
      Do J = 1,NOcc
      Do I = 1,NOcc
        Do B = NOcc+1,NSO
        Do A = NOcc+1,NSO
            E2 = E2 + ERI(I,J,A,B)*F14*(T2(A,B,I,J)           &
                + T1(A,I)*T1(B,J) - T1(B,I)*T1(A,J))
        End Do
        End Do
      End Do
      End Do
      !$omp end parallel do
      ECorr = ECorr + E2
      End Subroutine CCEnergy

      Subroutine CCEnergyBlocks(T1,T2,Fock,ECorr,NOcc,NSO,ConjH)
!  E = f_ia t_ai + 1/4 <ij||ab> (t_ij^ab + t_a^i t_b^j - t_b^i t_a^j), <ij||ab> from ERIBlocks
      Implicit None
      Integer,           Intent(In)  :: NOcc,NSO
      Complex (Kind=pr), Intent(In)  :: Fock(NSO,NSO)
      Complex (Kind=pr), Intent(In)  :: T1(NOcc+1:NSO,NOcc)
      Complex (Kind=pr), Intent(In)  :: T2(NOcc+1:NSO,NOcc+1:NSO,NOcc,NOcc)
      Logical,           Intent(In)  :: ConjH
      Complex (Kind=pr), Intent(Out) :: ECorr
      Integer :: I, J, A, B
      Complex (Kind=pr) :: E2
      If (.not. ERIBlocksReady .or. (ERIBlocksConj .neqv. ConjH)) Stop "CCEnergyBlocks: ERI blocks not set up"
      ECorr = sum(Fock(1:NOcc,NOcc+1:NSO) * transpose(T1))
      E2 = Zero
      !$omp parallel do reduction(+:E2) private(I,J,A,B)
      Do J = 1,NOcc
      Do I = 1,NOcc
        Do B = NOcc+1,NSO
        Do A = NOcc+1,NSO
            E2 = E2 + Voovv(I,J,A-NOcc,B-NOcc)*F14*(T2(A,B,I,J)           &
                + T1(A,I)*T1(B,J) - T1(B,I)*T1(A,J))
        End Do
        End Do
      End Do
      End Do
      !$omp end parallel do
      ECorr = ECorr + E2
      End Subroutine CCEnergyBlocks

      Subroutine ASymm2(T2,NOcc,NSO)
      Implicit None
      Integer,           Intent(In)    :: NOcc, NSO
      Complex (Kind=pr), Intent(InOut) :: T2(NOcc+1:NSO,NOcc+1:NSO,NOcc,NOcc)
      Complex (Kind=pr), Allocatable   :: X2(:,:,:,:)
      Integer   :: I, J, A, B
      Integer   :: IAlloc
      Allocate(X2(NOcc+1:NSO,NOcc+1:NSO,NOcc,NOcc), Stat=IAlloc)
      If(IAlloc /= 0) Stop "Could not allocate in ASymm2"
      !$omp parallel do private(I,J,A,B)
      Do J = 1, NOcc
      Do I = 1, NOcc
      Do B = NOcc+1, NSO
      Do A = NOcc+1, NSO
        X2(A,B,I,J) = T2(A,B,I,J) - T2(A,B,J,I)
      EndDo
      EndDo
      EndDo
      EndDo
      !$omp end parallel do
      !$omp parallel do private(I,J,A,B)
      Do J = 1, NOcc
      Do I = 1, NOcc
      Do B = NOcc+1, NSO
      Do A = NOcc+1, NSO
        T2(A,B,I,J) = X2(A,B,I,J) - X2(B,A,I,J)
      EndDo
      EndDo
      EndDo
      EndDo
      !$omp end parallel do
      Deallocate(X2, Stat=IAlloc)
      If(IAlloc /= 0) Stop "Could not deallocate in ASymm2"
      Return
      End Subroutine ASymm2


      Subroutine ASymm3(T3,NOcc,NSO)
      Implicit None
      Integer,           Intent(In)    :: NOcc, NSO
      Complex (Kind=pr), Intent(InOut) :: T3(NOcc+1:NSO,NOcc+1:NSO,NOcc+1:NSO,NOcc,NOcc,NOcc)
      Complex (Kind=pr), Allocatable   :: X3(:,:,:,:,:,:)
      Integer   :: I, J, K, A, B, C
      Integer   :: IAlloc
      Allocate(X3(NOcc+1:NSO,NOcc+1:NSO,NOcc+1:NSO,NOcc,NOcc,NOcc), Stat=IAlloc)
      If(IAlloc /= 0) Stop "Could not allocate in ASymm3"
      Do I = 1, NOcc
      Do J = 1, NOcc
      Do K = 1, NOcc
      Do A = NOcc+1, NSO
      Do B = NOcc+1, NSO
      Do C = NOcc+1, NSO
        X3(A,B,C,I,J,K) =   &
            T3(A,B,C,I,J,K) + T3(A,B,C,J,K,I) + T3(A,B,C,K,I,J) &
          - T3(A,B,C,I,K,J) - T3(A,B,C,K,J,I) - T3(A,B,C,J,I,K)
      EndDo
      EndDo
      EndDo
      EndDo
      EndDo
      EndDo
      Do I = 1, NOcc
      Do J = 1, NOcc
      Do K = 1, NOcc
      Do A = NOcc+1, NSO
      Do B = NOcc+1, NSO
      Do C = NOcc+1, NSO
        T3(A,B,C,I,J,K) =   &
            X3(A,B,C,I,J,K) + X3(B,C,A,I,J,K) + X3(C,A,B,I,J,K) &
          - X3(A,C,B,I,J,K) - X3(C,B,A,I,J,K) - X3(B,A,C,I,J,K)
      EndDo
      EndDo
      EndDo
      EndDo
      EndDo
      EndDo
      Deallocate(X3, Stat=IAlloc)
      If(IAlloc /= 0) Stop "Could not deallocate in ASymm3"
      Return
      End Subroutine ASymm3

      Subroutine APerm(M,sort)
!     sort=1: antisymmetrise the first two indices, sort=2: the last two
      Implicit None
      Integer,           Intent(In)    :: sort
      Complex (Kind=pr), Intent(InOut) :: M(:,:,:,:)
      Integer   :: I, J, NA, NB
      Complex (kind=pr), allocatable :: tmp(:,:,:,:)
      NA = size(M, dim=1)
      NB = size(M, dim=4)
      Allocate(tmp(NA,NA,NB,NB))
      tmp = M
      If (sort == 1) then
        !$omp parallel do private(I,J)
        Do J = 1, NB
        Do I = 1, NB
            M(:,:,I,J) = tmp(:,:,I,J) - transpose(tmp(:,:,I,J))
        End Do
        End Do
        !$omp end parallel do
      Else if (sort == 2) then
        !$omp parallel do private(I,J)
        Do J = 1, NB
        Do I = 1, NB
            M(:,:,I,J) = tmp(:,:,I,J) - tmp(:,:,J,I)
        End Do
        End Do
        !$omp end parallel do
      End If
      Deallocate(tmp)
      End Subroutine APerm

!=======================================================================!
!  CCSD residual (Stanton-Gauss form), same equations as the original    !
!  CCRes12/CCItms but:                                                   !
!    * integral blocks come from ERIBlocks (extracted once per run),     !
!    * no Wvvvv: the ladder 1/2 tau_ij^ef <ab||ef> is contracted straight !
!      from H2 (strided ZGEMM, no copy), the -P(ab) t_b^m <am||ef> part  !
!      goes through Z_am^ij = <am||ef> tau_ij^ef (o^3 v^3), and the      !
!      1/8 tau<mn||ef>tau part is folded into Woooo (1/2 instead of 1/4),!
!    * the O(o^3 v^3) ring terms are ZGEMMs on permuted copies,          !
!    * all temporaries are o^2 v^2 or smaller.                           !
!=======================================================================!
      Subroutine CCRes12(Fock,H2,T1,T2,res1,res2,NOcc,NVrt,NSO,ConjH)
!  ConjH = .True. means the Hamiltonian to use is conjg(H2): the ERIBlocks
!  must then hold the conjugated blocks (SetUpERIBlocks(...,Conj=.True.)),
!  and the ladder reads conjg(<ab||ef>) = <ef||ab> from H2 (Hermiticity).
      Implicit None
      Logical, Optional, Intent(In)    :: ConjH
      Integer,           Intent(In)    :: NOcc,NVrt,NSO
      Complex (Kind=pr), Intent(In)    :: Fock(NSO,NSO), H2(NSO,NSO,NSO,NSO)
      Complex (Kind=pr), Intent(In)    :: T1(NVrt,NOcc)
      Complex (Kind=pr), Intent(In)    :: T2(NVrt,NVrt,NOcc,NOcc)
      Complex (Kind=pr), Intent(Out)   :: res1(NVrt,NOcc)
      Complex (Kind=pr), Intent(Out)   :: res2(NVrt,NVrt,NOcc,NOcc)
      Integer :: o, v, no2, nv2, ov, a, b, c, e, f, i, j, m, n, IAlloc
      Complex (Kind=pr), Allocatable :: tau(:,:,:,:), tb(:,:,:,:), itm4(:,:,:,:)
      Complex (Kind=pr), Allocatable :: Fvv(:,:), Foo(:,:), Fov(:,:), Xvv(:,:), Xoo(:,:)
      Complex (Kind=pr), Allocatable :: Woooo(:,:,:,:), Wovvo(:,:,:,:)
      Complex (Kind=pr), Allocatable :: P1(:,:,:,:), P2(:,:,:,:), Y(:,:,:,:), G(:,:,:,:)
      Complex (Kind=pr), Allocatable :: Z(:,:,:,:), Q(:,:,:,:)
      Complex (Kind=pr) :: x
      Logical :: CH
      o = NOcc; v = NVrt; no2 = o*o; nv2 = v*v; ov = o*v
      CH = .False.
      If (Present(ConjH)) CH = ConjH
      If (.not. ERIBlocksReady .or. nOB /= o .or. nSB /= NSO .or. (ERIBlocksConj .neqv. CH)) &
        Call SetUpERIBlocks(H2,NOcc,NSO,Conj=CH)

      Allocate(tau(v,v,o,o), tb(v,v,o,o), itm4(v,v,o,o), Fvv(v,v), Foo(o,o), Fov(o,v), &
               Xvv(v,v), Xoo(o,o), Woooo(o,o,o,o), Wovvo(o,v,v,o), Stat=IAlloc)
      If (IAlloc /= 0) Stop "Could not allocate in CCRes12"

!---- tau = T2 + t_a^i t_b^j - t_b^i t_a^j ; tb = T2 + (..)/2 ; itm4 = T2/2 + t_a^i t_b^j
      !$omp parallel do private(a,b,i,j,x)
      Do j = 1, o
      Do i = 1, o
      Do b = 1, v
      Do a = 1, v
        x = T1(a,i)*T1(b,j) - T1(b,i)*T1(a,j)
        tau(a,b,i,j)  = T2(a,b,i,j) + x
        tb(a,b,i,j)   = T2(a,b,i,j) + F12*x
        itm4(a,b,i,j) = F12*T2(a,b,i,j) + T1(a,i)*T1(b,j)
      End Do
      End Do
      End Do
      End Do
      !$omp end parallel do

!---- Fvv(a,e) = f_ae - 1/2 t_a^m f_me + t_f^m <am||ef> - 1/2 tb_mn^af <mn||ef>
      Fvv = Fock(o+1:NSO,o+1:NSO) - F12*matmul(T1,Fock(1:o,o+1:NSO))
      !$omp parallel do private(e,m,f)
      Do e = 1, v
        Do m = 1, o
        Do f = 1, v
          Fvv(:,e) = Fvv(:,e) + T1(f,m)*Vvovv(:,m,e,f)
        End Do
        End Do
      End Do
      !$omp end parallel do
      ! tb viewed as (a | f,m,n), Vfmne as (f,m,n | e)
      Call PZGEMM('N','N',v,v,v*no2,-F12*ZOne,tb,v,Vfmne,v*no2,ZOne,Fvv,v)

!---- Foo(m,i) = f_mi + 1/2 t_e^i f_me + t_e^n <mn||ie> + 1/2 tb_in^ef <mn||ef>
      Foo = Fock(1:o,1:o) + F12*matmul(Fock(1:o,o+1:NSO),T1)
      Do e = 1, v
      Do n = 1, o
        Do i = 1, o
          Foo(:,i) = Foo(:,i) + T1(e,n)*Vooov(:,n,i,e)
        End Do
      End Do
      End Do
      Allocate(P1(o,v,v,o))          ! P1(n,e,f,i) = tb(e,f,i,n)
      !$omp parallel do private(i,n,e,f)
      Do i = 1, o
      Do f = 1, v
      Do e = 1, v
      Do n = 1, o
        P1(n,e,f,i) = tb(e,f,i,n)
      End Do
      End Do
      End Do
      End Do
      !$omp end parallel do
      ! Voovv as (m | n,e,f), P1 as (n,e,f | i)
      Call PZGEMM('N','N',o,o,o*nv2,F12*ZOne,Voovv,o,P1,o*nv2,ZOne,Foo,o)
      Deallocate(P1)

!---- Fov(m,e) = f_me + t_f^n <mn||ef>
      Fov = Fock(1:o,o+1:NSO)
      Do f = 1, v
      Do n = 1, o
        Do e = 1, v
          Fov(:,e) = Fov(:,e) + T1(f,n)*Voovv(:,n,e,f)
        End Do
      End Do
      End Do

!---- Woooo(m,n,i,j) = <mn||ij> + P(ij) t_e^j <mn||ie> + 1/2 tau_ij^ef <mn||ef>
!     (1/2 instead of 1/4: absorbs the 1/8 tau<mn||ef>tau piece of the ladder)
      Allocate(Y(o,o,o,o))
      Y = Zero
      Call ZGEMM('N','N',o*no2,o,v,ZOne,Vooov,o*no2,T1,v,ZOne,Y,o*no2)   ! Y(m,n,i,j) = <mn||ie> t_e^j
      !$omp parallel do private(i,j)
      Do j = 1, o
      Do i = 1, o
        Woooo(:,:,i,j) = Voooo(:,:,i,j) + Y(:,:,i,j) - Y(:,:,j,i)
      End Do
      End Do
      !$omp end parallel do
      Deallocate(Y)
      Call PZGEMM('N','N',no2,no2,nv2,F12*ZOne,Voovv,no2,tau,nv2,ZOne,Woooo,no2)

!---- Wovvo(m,b,e,j) = <mb||ej> + t_f^j <mb||ef> - t_b^n <mn||ej>
!                      - (1/2 t_jn^fb + t_j^f t_n^b) <mn||ef>
      Wovvo = Vovvo
      Allocate(Y(v,o,v,o))            ! Y(b,m,e,j) = sum_f <bm||ef> t_f^j
      Y = Zero
      Call PZGEMM('N','N',v*ov,o,v,ZOne,Vvovv,v*ov,T1,v,ZOne,Y,v*ov)
      !$omp parallel do private(m,b,e,j)
      Do j = 1, o
      Do e = 1, v
      Do b = 1, v
      Do m = 1, o
        Wovvo(m,b,e,j) = Wovvo(m,b,e,j) - Y(b,m,e,j)       ! <mb||ef> = -<bm||ef>
      End Do
      End Do
      End Do
      End Do
      !$omp end parallel do
      Deallocate(Y)
      ! - t_b^n <mn||ej> = + t_b^n <mn||je> :  Wovvo(m,b,e,j) += sum_n Vooov(m,n,j,e) t1(b,n)
      !$omp parallel do private(e,j)
      Do j = 1, o
      Do e = 1, v
        Call ZGEMM('N','T',o,v,o,ZOne,Vooov(1,1,j,e),o,T1,v,ZOne,Wovvo(1,1,e,j),o)
      End Do
      End Do
      !$omp end parallel do
      ! - sum_{n,f} itm4(f,b,j,n) <mn||ef> : A(m,e | n,f) = Vmenf, B(n,f | b,j) = itm4(f,b,j,n)
      Allocate(P1(o,v,v,o), Y(o,v,v,o))   ! P1(n,f,b,j) ; Y(m,e,b,j)
      !$omp parallel do private(n,f,b,j)
      Do j = 1, o
      Do b = 1, v
      Do f = 1, v
      Do n = 1, o
        P1(n,f,b,j) = itm4(f,b,j,n)
      End Do
      End Do
      End Do
      End Do
      !$omp end parallel do
      Y = Zero
      Call PZGEMM('N','N',ov,ov,ov,ZOne,Vmenf,ov,P1,ov,ZOne,Y,ov)
      !$omp parallel do private(m,b,e,j)
      Do j = 1, o
      Do e = 1, v
      Do b = 1, v
      Do m = 1, o
        Wovvo(m,b,e,j) = Wovvo(m,b,e,j) - Y(m,e,b,j)
      End Do
      End Do
      End Do
      End Do
      !$omp end parallel do
      Deallocate(P1, Y, itm4)

!======================= res1 ==========================================!
!  r1(a,i) = f_ai + Fvv(a,e) t_e^i - t_a^m Foo(m,i) + Fov(m,e) t_im^ae
!          + t_e^m <am||ie> - 1/2 t_mn^af <mn||if> + 1/2 t_im^ef <am||ef>
      res1 = Fock(o+1:NSO,1:o) + matmul(Fvv,T1) - matmul(T1,Foo)
      Do m = 1, o
      Do e = 1, v
        Do i = 1, o
          res1(:,i) = res1(:,i) + Fov(m,e)*T2(:,e,i,m) + T1(e,m)*Vovvo(m,:,e,i)
        End Do
      End Do
      End Do
      ! -1/2 sum_{f,m,n} T2(a,f,m,n) <mn||if> : T2 as (a | f,m,n), Vfmni as (f,m,n | i)
      Call PZGEMM('N','N',v,o,v*no2,-F12*ZOne,T2,v,Vfmni,v*no2,ZOne,res1,v)
      ! +1/2 sum_{m,e,f} <am||ef> T2(e,f,i,m) : Vvovv as (a | m,e,f), P2(m,e,f,i) = T2(e,f,i,m)
      Allocate(P2(o,v,v,o))
      !$omp parallel do private(m,e,f,i)
      Do i = 1, o
      Do f = 1, v
      Do e = 1, v
      Do m = 1, o
        P2(m,e,f,i) = T2(e,f,i,m)
      End Do
      End Do
      End Do
      End Do
      !$omp end parallel do
      Call PZGEMM('N','N',v,o,o*nv2,F12*ZOne,Vvovv,v,P2,o*nv2,ZOne,res1,v)
      Deallocate(P2)

!======================= res2 ==========================================!
      res2 = Vvvoo
!---- + P(ab) T2(a,e,i,j) Xvv(b,e),  Xvv = Fvv - 1/2 t_b^m Fov(m,e)
      Xvv = Fvv - F12*matmul(T1,Fov)
      Allocate(G(v,v,o,o))
      G = Zero
      !$omp parallel do private(i,j)
      Do j = 1, o
      Do i = 1, o
        Call ZGEMM('N','T',v,v,v,ZOne,T2(1,1,i,j),v,Xvv,v,ZOne,G(1,1,i,j),v)
      End Do
      End Do
      !$omp end parallel do
      !$omp parallel do private(i,j)
      Do j = 1, o
      Do i = 1, o
        res2(:,:,i,j) = res2(:,:,i,j) + G(:,:,i,j) - transpose(G(:,:,i,j))
      End Do
      End Do
      !$omp end parallel do
!---- - P(ij) T2(a,b,i,m) Xoo(m,j),  Xoo = Foo + 1/2 t_e^j Fov(m,e)
      Xoo = Foo + F12*matmul(Fov,T1)
      G = Zero
      Call PZGEMM('N','N',nv2*o,o,o,ZOne,T2,nv2*o,Xoo,o,ZOne,G,nv2*o)
      !$omp parallel do private(i,j)
      Do j = 1, o
      Do i = 1, o
        res2(:,:,i,j) = res2(:,:,i,j) - G(:,:,i,j) + G(:,:,j,i)
      End Do
      End Do
      !$omp end parallel do
!---- + 1/2 tau_mn^ab Woooo(m,n,i,j)
      Call PZGEMM('N','N',nv2,no2,no2,F12*ZOne,tau,nv2,Woooo,no2,ZOne,res2,nv2)
!---- + 1/2 tau_ij^ef <ab||ef>  (ladder straight from H2; thread owns a b-slab of res2)
      If (.not. CH) Then
        !$omp parallel do schedule(dynamic) private(b,f)
        Do b = 1, v
          Do f = 1, v
            Call ZGEMM('N','N',v,no2,v,F12*ZOne,H2(o+1,o+b,o+1,o+f),NSO*NSO, &
                       tau(1,f,1,1),nv2,ZOne,res2(1,b,1,1),nv2)
          End Do
        End Do
        !$omp end parallel do
      Else
        ! conjg(<ab||ef>) = <ef||ab> : A(e,a) = H2(o+e,o+f,o+a,o+b), used transposed
        !$omp parallel do schedule(dynamic) private(b,f)
        Do b = 1, v
          Do f = 1, v
            Call ZGEMM('T','N',v,no2,v,F12*ZOne,H2(o+1,o+f,o+1,o+b),NSO*NSO, &
                       tau(1,f,1,1),nv2,ZOne,res2(1,b,1,1),nv2)
          End Do
        End Do
        !$omp end parallel do
      End If
!---- - 1/2 P(ab) t_b^m Z(a,m,i,j),  Z(a,m,i,j) = <am||ef> tau_ij^ef
      Allocate(Z(v,o,o,o))
      Z = Zero
      Call PZGEMM('N','N',ov,no2,nv2,ZOne,Vvovv,ov,tau,nv2,ZOne,Z,ov)
      G = Zero
      !$omp parallel do private(i,j)
      Do j = 1, o
      Do i = 1, o
        Call ZGEMM('N','T',v,v,o,ZOne,Z(1,1,i,j),v,T1,v,ZOne,G(1,1,i,j),v)   ! G(a,b,i,j) = Z(a,m,i,j) t_b^m
      End Do
      End Do
      !$omp end parallel do
      !$omp parallel do private(i,j)
      Do j = 1, o
      Do i = 1, o
        res2(:,:,i,j) = res2(:,:,i,j) - F12*(G(:,:,i,j) - transpose(G(:,:,i,j)))
      End Do
      End Do
      !$omp end parallel do
      Deallocate(Z)
!---- + P(ij)P(ab) [ T2(a,e,i,m) Wovvo(m,b,e,j) - t_e^i t_a^m <mb||ej> ]
      Allocate(P1(v,o,v,o), P2(v,o,v,o), Y(v,o,v,o))   ! P1(a,i,e,m)=T2(a,e,i,m); P2(e,m,b,j)=Wovvo(m,b,e,j); Y(a,i,b,j)
      !$omp parallel do private(a,e,i,m)
      Do m = 1, o
      Do e = 1, v
      Do i = 1, o
      Do a = 1, v
        P1(a,i,e,m) = T2(a,e,i,m)
      End Do
      End Do
      End Do
      End Do
      !$omp end parallel do
      !$omp parallel do private(m,b,e,j)
      Do j = 1, o
      Do b = 1, v
      Do m = 1, o
      Do e = 1, v
        P2(e,m,b,j) = Wovvo(m,b,e,j)
      End Do
      End Do
      End Do
      End Do
      !$omp end parallel do
      Y = Zero
      Call PZGEMM('N','N',ov,ov,ov,ZOne,P1,ov,P2,ov,ZOne,Y,ov)
      ! Q(m,b,i,j) = sum_e Vovvo(m,b,e,j) t_e^i ;  G(a,b,i,j) = sum_m t_a^m Q(m,b,i,j)
      Allocate(Q(o,v,o,o))
      Q = Zero
      !$omp parallel do private(j)
      Do j = 1, o
        Call ZGEMM('N','N',ov,o,v,ZOne,Vovvo(1,1,1,j),ov,T1,v,ZOne,Q(1,1,1,j),ov)
      End Do
      !$omp end parallel do
      G = Zero
      Call PZGEMM('N','N',v,v*no2,o,ZOne,T1,v,Q,o,ZOne,G,v)
      !$omp parallel do private(a,b,i,j)
      Do j = 1, o
      Do i = 1, o
      Do b = 1, v
      Do a = 1, v
        res2(a,b,i,j) = res2(a,b,i,j) + Y(a,i,b,j) - Y(a,j,b,i) - Y(b,i,a,j) + Y(b,j,a,i) &
                      - G(a,b,i,j) + G(a,b,j,i) + G(b,a,i,j) - G(b,a,j,i)
      End Do
      End Do
      End Do
      End Do
      !$omp end parallel do
      Deallocate(P1, P2, Y, Q)
!---- + P(ij) t_e^i <ab||ej>
      G = Zero
      !$omp parallel do private(j)
      Do j = 1, o
        Call ZGEMM('N','N',nv2,o,v,ZOne,Vvvvo(1,1,1,j),nv2,T1,v,ZOne,G(1,1,1,j),nv2)  ! G(a,b,i,j)
      End Do
      !$omp end parallel do
      !$omp parallel do private(i,j)
      Do j = 1, o
      Do i = 1, o
        res2(:,:,i,j) = res2(:,:,i,j) + G(:,:,i,j) - G(:,:,j,i)
      End Do
      End Do
      !$omp end parallel do
!---- - P(ab) t_a^m <mb||ij>
      G = Zero
      Call PZGEMM('N','N',v,v*no2,o,ZOne,T1,v,Vovoo,o,ZOne,G,v)              ! G(a,b,i,j) = t_a^m <mb||ij>
      !$omp parallel do private(i,j)
      Do j = 1, o
      Do i = 1, o
        res2(:,:,i,j) = res2(:,:,i,j) - G(:,:,i,j) + transpose(G(:,:,i,j))
      End Do
      End Do
      !$omp end parallel do
      Deallocate(G, tau, tb, Fvv, Foo, Fov, Xvv, Xoo, Woooo, Wovvo)
      End Subroutine CCRes12

    End Module CCRes
