    Module CCRes 
    Use Precision
    Use Constants
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
      ECorr = sum(Fock(1:NOcc,NOcc+1:NSO) * transpose(T1))
      Do A = NOcc+1,NSO
      Do I = 1,NOcc
        Do B = NOcc+1,NSO
        Do J = 1,NOcc
            ECorr = ECorr + ERI(I,J,A,B)*F14*(T2(A,B,I,J)           &
                + T1(A,I)*T1(B,J) - T1(B,I)*T1(A,J))
        End Do
        End Do
      End Do
      End Do
      End Subroutine CCEnergy

      Subroutine ASymm2(T2,NOcc,NSO)
      Implicit None
      Integer,           Intent(In)    :: NOcc, NSO
      Complex (Kind=pr), Intent(InOut) :: T2(NOcc+1:NSO,NOcc+1:NSO,NOcc,NOcc)
      Complex (Kind=pr), Allocatable   :: X2(:,:,:,:)
      Integer   :: I, J, A, B
      Integer   :: IAlloc
! Allocate necessary space
      Allocate(X2(NOcc+1:NSO,NOcc+1:NSO,NOcc,NOcc), Stat=IAlloc)
      If(IAlloc /= 0) Stop "Could not allocate in ASymm2"
! Antisymmetrize occ 
      !$omp parallel do
      Do I = 1, NOcc
      Do J = 1, NOcc 
      Do A = NOcc+1, NSO 
      Do B = NOcc+1, NSO 
        X2(A,B,I,J) = T2(A,B,I,J) - T2(A,B,J,I)
      EndDo
      EndDo
      EndDo
      EndDo
      !$omp end parallel do
! Antisymmetrize vir 
      !$omp parallel do
      Do I = 1, NOcc
      Do J = 1, NOcc 
      Do A = NOcc+1, NSO 
      Do B = NOcc+1, NSO 
        T2(A,B,I,J) = X2(A,B,I,J) - X2(B,A,I,J)
      EndDo
      EndDo
      EndDo
      EndDo
      !$omp end parallel do
! Deallocate necessary space
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
! Allocate necessary space
      Allocate(X3(NOcc+1:NSO,NOcc+1:NSO,NOcc+1:NSO,NOcc,NOcc,NOcc), Stat=IAlloc)
      If(IAlloc /= 0) Stop "Could not allocate in ASymm3"
! Antisymmetrize occ 
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
! Antisymmetrize vir 
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
! Deallocate necessary space
      Deallocate(X3, Stat=IAlloc)
      If(IAlloc /= 0) Stop "Could not deallocate in ASymm3"
      Return
      End Subroutine ASymm3

      Subroutine APerm(M,sort)
      Implicit None
      Integer,           Intent(In)    :: sort
      Complex (Kind=pr), Intent(InOut) :: M(:,:,:,:)
      Integer   :: I, J, NA, NB
      Complex (kind=pr), allocatable :: tmp(:,:,:,:)
      integer :: nblk, n, n1
      integer, allocatable :: load1(:), is1(:), ie1(:)
      integer, allocatable :: load2(:), is2(:), ie2(:)
      NA = size(M, dim=1)
      NB = size(M, dim=4)
      Allocate(tmp(NA,NA,NB,NB))
      nblk = nb / blksize
      Allocate(load1(nblk), is1(nblk), ie1(nblk))
      load1 = nb / nblk
      load1(1:mod(nb,nblk)) = load1(1:mod(nb,nblk)) + 1
      Do n = 1, nblk
        is1(n) = sum(load1(1:n-1))
        ie1(n) = is1(n) + load1(n)
      End Do
      tmp = M
      If (sort == 1) then
        !$omp parallel do
        Do I = 1, NA
        Do J = 1, NA
            M(I,J,:,:) = tmp(I,J,:,:) - tmp(J,I,:,:)
        End Do
        End Do
        !$omp end parallel do
      Else if (sort == 2) then
        !$omp parallel do
        Do n = 1, nblk
        Do n1 = 1, nblk
        Do I = is1(n)+1, ie1(n)
        Do J = is1(n1)+1, ie1(n1)
            M(:,:,I,J) = tmp(:,:,I,J) - tmp(:,:,J,I)
        End Do
        End Do
        End Do
        End Do
        !$omp end parallel do
      End If
      End Subroutine APerm

      Subroutine CCRes12(Fock,H2,T1,T2,res1,res2,NOcc,NVrt,NSO)
      Implicit None
      Integer,           Intent(In)    :: NOcc,NVrt,NSO
      Complex (Kind=pr), Intent(In)    :: Fock(NSO,NSO), H2(NSO,NSO,NSO,NSO)
      Complex (Kind=pr), Intent(In)    :: T1(NVrt,NOcc)
      Complex (Kind=pr), Intent(In)    :: T2(NVrt,NVrt,NOcc,NOcc)
      Complex (Kind=pr), Intent(Out)   :: res1(NOcc+1:NSO,NOcc)
      Complex (Kind=pr), Intent(Out)   :: res2(NOcc+1:NSO,NOcc+1:NSO,NOcc,NOcc)
      Integer :: a,b,c,d,i,j,k,l
      Complex (Kind=pr)                :: Foo(NOcc,NOcc), Fov(NOcc,NVrt)
      Complex (Kind=pr)                :: Fvv(NVrt,NVrt), tau(NVrt,NVrt,NOcc,NOcc)
      Complex (Kind=pr)                :: Woooo(NOcc,NOcc,NOcc,NOcc)
      Complex (Kind=pr)                :: Wvvvv(NVrt,NVrt,NVrt,NVrt)
      Complex (Kind=pr)                :: Wovvo(NOcc,NVrt,NVrt,NOcc)
      Complex (Kind=pr), allocatable :: r1(:,:), r2(:,:,:,:)
      Complex (Kind=pr), allocatable :: t2tmp(:,:,:,:), r2tmp(:,:,:,:)
      Complex (Kind=pr), allocatable :: Hvoov(:,:,:,:), Hooov(:,:,:,:)
      Complex (Kind=pr), allocatable :: Hvovv(:,:,:,:), Hovvo(:,:,:,:)
      Complex (Kind=pr), allocatable :: Hvvvo(:,:,:,:), Hovoo(:,:,:,:)
      Complex (Kind=pr), allocatable :: itm1(:,:),itm2(:,:,:,:)
      Complex (Kind=pr), allocatable :: itm3(:,:),itm4(:,:,:,:)
      integer :: nblkv, m, m1, no2, nv2, ov
      integer, allocatable :: loadv(:), isv(:), iev(:)
      Integer         :: cnt1, cnt2, clock_rate, clock_max
      nblkv = nvrt / blksize
      Allocate(loadv(nblkv), isv(nblkv), iev(nblkv))
      loadv = nvrt / nblkv
      loadv(1:mod(nvrt,nblkv)) = loadv(1:mod(nvrt,nblkv)) + 1
      Do m = 1, nblkv
        isv(m) = sum(loadv(1:m-1))
        iev(m) = isv(m) + loadv(m)
      End Do
      no2 = NOcc**2
      nv2 = NVrt**2
      ov = NOcc*NVrt
      Allocate(r1(NVrt,NOcc), r2(NVrt,NVrt,NOcc,NOcc), &
               t2tmp(NOcc,NOcc,NVrt,NVrt), r2tmp(NOcc,NOcc,NVrt,NVrt))
      ! build intermediates
      Call CCItms(Fock,H2,T1,T2,tau,Fvv,Foo,Fov, &
                  Woooo,Wvvvv,Wovvo,NOcc,NVrt,NSO)
      ! start
      r2tmp = Zero
      Call Transpose2(T2,t2tmp,nv2,no2)
      Allocate(Hvoov(NVrt,NOcc,NOcc,NVrt), Hooov(NOcc,NOcc,NOcc,NVrt), &
               Hvovv(NVrt,NOcc,NVrt,NVrt))
      Hvoov = H2(NOcc+1:NSO,:NOcc,:NOcc,NOcc+1:NSO)
      Hooov = H2(:NOcc,:NOcc,:NOcc,NOcc+1:NSO)
      Hvovv = H2(NOcc+1:NSO,:NOcc,NOcc+1:NSO,NOcc+1:NSO)
      ! build res1
      r1 = Fock(NOcc+1:NSO,:NOcc) + matmul(Fvv,T1) - matmul(t1,Foo)
      !$omp parallel do
      Do a = 1, NVrt
      Do i = 1, NOcc
        r1(a,i) = r1(a,i) + sum(Fov * t2tmp(i,:,a,:))
        r1(a,i) = r1(a,i) + sum(transpose(t1) * Hvoov(a,:,i,:))
      End Do
      End Do
      !$omp end parallel do
      !$omp parallel do reduction(+:r1)
      Do m1 = 1, nblkv
      Do b = isv(m1)+1, iev(m1)
        Call At_dot_B(t2tmp(:,:,:,b),-Hooov(:,:,:,b)/2,r1,no2,nvrt,nocc)
      End Do
      End Do
      !$omp end parallel do
      !$omp parallel do reduction(+:r1)
      Do m1 = 1, nblkv
      Do b = isv(m1)+1, iev(m1)
        Call A_dot_Bt(Hvovv(:,:,:,b),t2tmp(:,:,:,b)/2,r1,nvrt,ov,nocc)
      End Do
      End Do
      !$omp end parallel do
      ! build res2
      r2 = H2(NOcc+1:NSO,NOcc+1:NSO,:NOcc,:NOcc)
      Allocate(itm1(NVrt,NVrt))
      itm1 = Fvv - matmul(t1,Fov) / 2

      Allocate(itm2(NOcc,NOcc,NVrt,NVrt))
      itm2 = zero
      !$omp parallel do
      Do m = 1, nblkv
      Do a = isv(m)+1, iev(m)
        Call A_dot_Bt(t2tmp(:,:,a,:),itm1,itm2(:,:,a,:),no2,nvrt,nvrt)
      End Do
      End Do
      !$omp end parallel do
      Call APerm(itm2,2)
      r2tmp = r2tmp + itm2
      deallocate(itm1)
      Allocate(itm3(NOcc,NOcc))
      itm3 = Foo + matmul(Fov,t1) / 2
      itm2 = zero
      !$omp parallel do
      Do m = 1, nblkv
      Do b = isv(m)+1, iev(m)
      Do a = 1, NVrt
        itm2(:,:,a,b) = matmul(t2tmp(:,:,a,b),itm3)
      End Do
      End Do
      End Do
      !$omp end parallel do
      deallocate(itm3)
      Call APerm(itm2,1)
      r2tmp = r2tmp - itm2

      Allocate(itm4(NOcc,NOcc,NOcc,NVrt))
      itm4 = zero
      Call Transpose2(H2(:NOcc,NOcc+1:NSO,NOcc+1:NSO,:NOcc),Hvoov,ov,ov)
!      Hvoov = conjg(Hvoov)
      !$omp parallel do
      Do m = 1, nblkv
      Do b = isv(m)+1, iev(m)
        Call At_dot_B(t1,Hvoov(:,:,:,b),itm4(:,:,:,b),NVrt,NOcc,no2)
      End Do
      End Do
      !$omp end parallel do
      itm2 = zero
      !$omp parallel do
      Do m = 1, nblkv
      Do b = isv(m)+1, iev(m)
        Call A_dot_Bt(itm4(:,:,:,b),-t1,itm2(:,:,:,b),no2,nocc,nvrt)
      End Do
      End Do
      !$omp end parallel do
      Deallocate(itm4,Hvoov)
      !$omp parallel do
      Do m1 = 1, nblkv
      Do m = 1, nblkv
      Do i = 1, NOcc
      Do j = 1, NOcc
      Do a = isv(m)+1, iev(m)
      Do b = isv(m1)+1, iev(m1)
        itm2(i,j,a,b) = itm2(i,j,a,b) + sum(t2tmp(i,:,a,:) * Wovvo(:,b,:,j))
      End Do
      End Do
      End Do
      End Do
      End Do
      End Do
      !$omp end parallel do
      Call APerm(itm2,1)
      Call APerm(itm2,2)
      r2tmp = r2tmp + itm2

      Call Transpose1(Wvvvv,nv2)
      !$omp parallel do
      Do m = 1, nblkv
      Do b = isv(m)+1, iev(m)
        Call At_dot_B(tau/2,Wvvvv(:,:,:,b),r2tmp(:,:,:,b),nv2,no2,NVrt)
      End Do
      End Do
      !$omp end parallel do
      Call Transpose2(tau/2,t2tmp,nv2,no2)
      !$omp parallel do
      Do m = 1, nblkv
      Do b = isv(m)+1, iev(m)
        Call At_dot_B(Woooo,t2tmp(:,:,:,b),r2tmp(:,:,:,b),no2,no2,NVrt)
      End Do
      End Do
      !$omp end parallel do

      !Hvovv = conjg(Hvovv)
      !Hooov = conjg(Hooov)
      Call Transpose2(H2(NOcc+1:NSO,NOcc+1:NSO,NOcc+1:NSO,:NOcc),Hvovv,nv2,ov)
      Call Transpose2(H2(:NOcc,NOcc+1:NSO,:NOcc,:NOcc),Hooov,ov,no2)
      itm2 = Zero
      !$omp parallel do
      Do m = 1, nblkv
      Do b = isv(m)+1, iev(m)
        Call At_dot_B(t1,Hvovv(:,:,:,b),itm2(:,:,:,b),NVrt,NOcc,ov)
      End Do
      End Do
      !$omp end parallel do
      Call APerm(itm2,1)
      r2tmp = r2tmp + itm2
      itm2 = zero
      !$omp parallel do
      Do m = 1, nblkv
      Do b = isv(m)+1, iev(m)
        Call A_dot_Bt(Hooov(:,:,:,b),t1,itm2(:,:,:,b),no2,NOcc,NVrt)
      End Do
      End Do
      !$omp end parallel do
      Call APerm(itm2,2)
      r2tmp = r2tmp - itm2
      res1(NOcc+1:NSO,:NOcc) = r1
      Call Transpose2(r2tmp,res2(NOcc+1:NSO,NOcc+1:NSO,:NOcc,:NOcc),no2,nv2)
      res2(NOcc+1:NSO,NOcc+1:NSO,:NOcc,:NOcc) = res2(NOcc+1:NSO,NOcc+1:NSO,:NOcc,:NOcc) + r2
      deallocate(r1,r2,itm2)
      End Subroutine CCRes12

      Subroutine CCItms(Fock,H2,T1,T2,tau,Fvv,Foo,Fov, &
                        Woooo,Wvvvv,Wovvo,NOcc,NVrt,NSO)
      Implicit None
      Integer,           Intent(In)    :: NOcc,NVrt,NSO
      Complex (Kind=pr), Intent(In)    :: Fock(NSO,NSO), H2(NSO,NSO,NSO,NSO)
      Complex (Kind=pr), Intent(In)    :: T1(NVrt,NOcc)
      Complex (Kind=pr), Intent(In)    :: T2(NVrt,NVrt,NOcc,NOcc)
      Complex (Kind=pr), Intent(Out)   :: tau(NVrt,NVrt,NOcc,NOcc)
      Complex (Kind=pr), Intent(Out)   :: Fvv(NVrt,NVrt), Foo(NOcc,NOcc)
      Complex (Kind=pr), Intent(Out)   :: Fov(NOcc,NVrt)
      Complex (Kind=pr), Intent(Out)   :: Woooo(NOcc,NOcc,NOcc,NOcc)
      Complex (Kind=pr), Intent(Out)   :: Wvvvv(NVrt,NVrt,NVrt,NVrt)
      Complex (Kind=pr), Intent(Out)   :: Wovvo(NOcc,NVrt,NVrt,NOcc)
      Integer                          :: a,b,c,d,i,j,k,l
      Complex (Kind=pr), allocatable   :: tb(:,:,:,:), itm1(:,:,:,:)
      Complex (Kind=pr), allocatable   :: itm2(:,:,:,:), itm3(:,:,:,:)
      Complex (Kind=pr), allocatable   :: itm4(:,:,:,:), Hvovv(:,:,:,:)
      Complex (Kind=pr), allocatable   :: Hoovv(:,:,:,:), Hooov(:,:,:,:)
      Complex (Kind=pr), allocatable   :: Hoovo(:,:,:,:), Hovvv(:,:,:,:)
      Complex (Kind=pr), allocatable   :: t2tmp(:,:,:,:)
      integer :: nblkv, nblko, m, m1, no2, nv2, ov
      integer, allocatable :: loadv(:), isv(:), iev(:)
      integer, allocatable :: loado(:), iso(:), ieo(:)
      nblkv = nvrt / blksize
      Allocate(loadv(nblkv), isv(nblkv), iev(nblkv))
      loadv = nvrt / nblkv
      loadv(1:mod(nvrt,nblkv)) = loadv(1:mod(nvrt,nblkv)) + 1
      Do m = 1, nblkv
        isv(m) = sum(loadv(1:m-1))
        iev(m) = isv(m) + loadv(m)
      End Do
      no2 = NOcc**2
      nv2 = NVrt**2
      ov = NOcc*NVrt
      Allocate(tb(NVrt,NVrt,NOcc,NOcc), itm1(NVrt,NVrt,NOcc,NOcc), &
               itm4(NVrt,NVrt,NOcc,NOcc), Hoovv(NOcc,NOcc,NVrt,NVrt), &
               Hoovo(NOcc,NOcc,NVrt,NOcc),Hooov(NOcc,NOcc,NOcc,NVrt))
      Hoovv = H2(:NOcc,:NOcc,NOcc+1:NSO,NOcc+1:NSO)
      Hooov = H2(:NOcc,:NOcc,:NOcc,NOcc+1:NSO)
      Hoovo = H2(:NOcc,:NOcc,NOcc+1:NSO,:NOcc)
      ! tau
      !$omp parallel do
      Do a = 1, NVrt
      Do b = 1, NVrt
      Do i = 1, NOcc
      Do j = 1, NOcc
        itm1(a,b,i,j) = t1(a,i) * t1(b,j)
      End Do
      End Do
      End Do
      End Do
      !$omp end parallel do
      itm4 = T2 / 2 + itm1
      Call APerm(itm1,1)
      Call APerm(itm1,2)
      tau = T2 + itm1 / 2
      tb = T2 + itm1 / 4
      deallocate(itm1)
      Allocate(t2tmp(NOcc,NOcc,NVrt,NVrt))
      ! Fvv
      Call Transpose2(tb/2,t2tmp,nv2,no2)
      Allocate(Hvovv(NVrt,NOcc,NVrt,NVrt))
      Hvovv = H2(NOcc+1:NSO,:NOcc,NOcc+1:NSO,NOcc+1:NSO)
      Fvv = Fock(NOcc+1:NSO,NOcc+1:NSO) - matmul(T1,Fock(:NOcc,NOcc+1:NSO)) / 2
      !$omp parallel do
      Do m = 1, nblkv
      Do a = 1, NVrt
      Do b = isv(m)+1, iev(m)
        Fvv(a,b) = Fvv(a,b) + sum(transpose(T1) * Hvovv(a,:,b,:))
      End Do
      End Do
      End Do
      !$omp end parallel do
      !$omp parallel do reduction(+:Fvv)
      Do m = 1, nblkv
      Do b = isv(m)+1, iev(m)
        Call At_dot_B(-t2tmp(:,:,:,b),Hoovv(:,:,:,b),Fvv,no2,nvrt,nvrt)
      End Do
      End Do
      !$omp end parallel do
      ! Foo
      Foo = Fock(:NOcc,:NOcc) + matmul(Fock(:NOcc,NOcc+1:NSO), T1) / 2
      !$omp parallel do
      Do i = 1, NOcc
      Do j = 1, NOcc
        Foo(i,j) = Foo(i,j) + sum(transpose(T1) * Hooov(i,:,j,:))
      End Do
      End Do
      !$omp end parallel do
      !$omp parallel do reduction(+:Foo)
      Do m = 1, nblkv
      Do b = isv(m)+1, iev(m)
        Call A_dot_Bt(Hoovv(:,:,:,b),t2tmp(:,:,:,b),Foo,nocc,ov,nocc)
      End Do
      End Do
      !$omp end parallel do
      deallocate(tb)
      ! Fov
      Fov = Fock(:NOcc,NOcc+1:NSO)
      !$omp parallel do
      Do m = 1, nblkv
      Do i = 1, NOcc
      Do a = isv(m)+1, iev(m)
        Fov(i,a) = Fov(i,a) + sum(transpose(T1) * Hoovv(i,:,a,:))
      End Do
      End Do
      End Do
      !$omp end parallel do
      Call Transpose2(tau/4,t2tmp,nv2,no2)
      ! Woooo
      Woooo = H2(:NOcc,:NOcc,:NOcc,:NOcc)
      Call A_dot_B(Hoovv,tau/4,Woooo,no2,nv2,no2)
      Allocate(itm2(NOcc,NOcc,NOcc,NOcc))
      itm2 = zero
      Call A_dot_B(Hooov,T1,itm2,NOcc**3,NVrt,NOcc)
      Call APerm(itm2,2)
      Woooo = Woooo + itm2
      deallocate(itm2,Hooov)
      ! Wvvvv
      Wvvvv = H2(NOcc+1:NSO,NOcc+1:NSO,NOcc+1:NSO,NOcc+1:NSO)
      !$omp parallel do
      Do m = 1, nblkv
      Do b = isv(m)+1, iev(m)
        Call At_dot_B(t2tmp,Hoovv(:,:,:,b),Wvvvv(:,:,:,b),no2,nv2,NVrt)
      End Do
      End Do
      !$omp end parallel do
      Allocate(itm3(NVrt,NVrt,NVrt,NVrt))
      itm3 = zero
      !$omp parallel do
      Do m1 = 1, nblkv
      Do m = 1, nblkv
      Do c = isv(m)+1, iev(m)
      Do d = isv(m1)+1, iev(m1)
        itm3(:,:,c,d) = itm3(:,:,c,d) + matmul(Hvovv(:,:,c,d),transpose(t1))
      End Do
      End Do
      End Do
      End Do
      !$omp end parallel do
      Call APerm(itm3,1)
      Wvvvv = Wvvvv - itm3
      deallocate(itm3,Hvovv,t2tmp)

      ! Wovvo
      Allocate(Hovvv(NOcc,NVrt,NVrt,NVrt))
      Hovvv = H2(:NOcc,NOcc+1:NSO,NOcc+1:NSO,NOcc+1:NSO)
      Wovvo = H2(:NOcc,NOcc+1:NSO,NOcc+1:NSO,:NOcc)
      !$omp parallel do
      Do m = 1, nblkv
      Do b = isv(m)+1, iev(m)
        Call A_dot_B(Hovvv(:,:,b,:),T1,Wovvo(:,:,b,:),ov,NVrt,NOcc)
      End Do
      End Do
      !$omp end parallel do
      !$omp parallel do
      Do m = 1, nblkv
      Do j = 1, NOcc
      Do c = isv(m)+1, iev(m)
        Wovvo(:,:,c,j) = Wovvo(:,:,c,j) - matmul(Hoovo(:,:,c,j),transpose(t1))
      End Do
      End Do
      End Do
      !$omp end parallel do
      !$omp parallel do
      Do m = 1, nblkv
      Do m1 = 1, nblkv
      Do i = 1, NOcc
      Do j = 1, NOcc
      Do b = isv(m)+1, iev(m)
      Do c = isv(m1)+1, iev(m1)
        Wovvo(i,b,c,j) = Wovvo(i,b,c,j) - sum(transpose(itm4(:,b,j,:)) * Hoovv(i,:,c,:))
      End Do
      End Do
      End Do
      End Do
      End Do
      End Do
      !$omp end parallel do
      deallocate(itm4,Hovvv,Hoovo,Hoovv)

      End Subroutine CCItms


    End Module CCRes 
