
   Module ERIBlocks
!
!  Constant integral blocks of the (Dirac, antisymmetrised) Hamiltonian
!  H2(p,q,r,s) = <pq||rs>, extracted ONCE per run and reused by every grid
!  point and every Broyden iteration (the old code re-sliced them out of the
!  NSO^4 array at every kernel evaluation).  Index order of each block is the
!  one needed for contiguous ZGEMM operands.  The vvvv block is never stored:
!  the particle-particle ladder reads it straight from H2 with a strided
!  leading dimension (see CCRes12).
!
   Use Precision
   Use Constants
   Implicit None
   Logical :: ERIBlocksReady = .False.
   Logical :: ERIBlocksT     = .False.
   Logical :: ERIBlocksConj  = .False.
   Integer :: nOB = 0, nVB = 0, nSB = 0
   Complex (Kind=pr), Allocatable :: Voooo(:,:,:,:)   ! <mn||ij>  (o,o,o,o)
   Complex (Kind=pr), Allocatable :: Vooov(:,:,:,:)   ! <mn||ie>  (o,o,o,v)
   Complex (Kind=pr), Allocatable :: Voovv(:,:,:,:)   ! <mn||ef>  (o,o,v,v)
   Complex (Kind=pr), Allocatable :: Vovvo(:,:,:,:)   ! <mb||ej>  (o,v,v,o)
   Complex (Kind=pr), Allocatable :: Vovoo(:,:,:,:)   ! <mb||ij>  (o,v,o,o)
   Complex (Kind=pr), Allocatable :: Vvovv(:,:,:,:)   ! <am||ef>  (v,o,v,v)
   Complex (Kind=pr), Allocatable :: Vvvvo(:,:,:,:)   ! <ab||ej>  (v,v,v,o)
   Complex (Kind=pr), Allocatable :: Vvvoo(:,:,:,:)   ! <ab||ij>  (v,v,o,o)
! permuted copies used as ZGEMM operands
   Complex (Kind=pr), Allocatable :: Vfmne(:,:,:,:)   ! Vfmne(f,m,n,e) = <mn||ef>
   Complex (Kind=pr), Allocatable :: Vfmni(:,:,:,:)   ! Vfmni(f,m,n,i) = <mn||if>
   Complex (Kind=pr), Allocatable :: Vmenf(:,:,:,:)   ! Vmenf(m,e,n,f) = <mn||ef>
! (T) blocks, only when requested
   Complex (Kind=pr), Allocatable :: Tvvvo(:,:,:,:)   ! Tvvvo(d,k,b,c) = <bc||dk>
   Complex (Kind=pr), Allocatable :: Tovoo(:,:,:,:)   ! Tovoo(l,j,k,a) = <la||jk>

   Contains

      Subroutine SetUpERIBlocks(H2,NOcc,NSO,WithT,Conj)
!     Conj = .True. stores the blocks of conjg(H2) (element-wise, no NSO^4 temporary)
      Implicit None
      Integer,           Intent(In) :: NOcc, NSO
      Complex (Kind=pr), Intent(In) :: H2(NSO,NSO,NSO,NSO)
      Logical, Optional, Intent(In) :: WithT, Conj
      Integer :: o, v, i, j, a, b, m, e, f, k, IAlloc
      Logical :: DoT, DoC
      DoT = .False.
      DoC = .False.
      If (Present(WithT)) DoT = WithT
      If (Present(Conj))  DoC = Conj
      If (ERIBlocksReady) Then
        If (nOB == NOcc .and. nSB == NSO .and. (ERIBlocksT .or. .not. DoT) .and. (ERIBlocksConj .eqv. DoC)) Return
        DoT = DoT .or. ERIBlocksT
        Call ShutDownERIBlocks
      End If
      o = NOcc
      v = NSO - NOcc
      nOB = o; nVB = v; nSB = NSO
      Allocate(Voooo(o,o,o,o), Vooov(o,o,o,v), Voovv(o,o,v,v), Vovvo(o,v,v,o), &
               Vovoo(o,v,o,o), Vvovv(v,o,v,v), Vvvvo(v,v,v,o), Vvvoo(v,v,o,o), &
               Vfmne(v,o,o,v), Vfmni(v,o,o,o), Vmenf(o,v,o,v), Stat=IAlloc)
      If (IAlloc /= 0) Stop "Could not allocate in SetUpERIBlocks"
      Voooo = H2(1:o,1:o,1:o,1:o)
      Vooov = H2(1:o,1:o,1:o,o+1:NSO)
      Voovv = H2(1:o,1:o,o+1:NSO,o+1:NSO)
      Vovvo = H2(1:o,o+1:NSO,o+1:NSO,1:o)
      Vovoo = H2(1:o,o+1:NSO,1:o,1:o)
      Vvvoo = H2(o+1:NSO,o+1:NSO,1:o,1:o)
      !$omp parallel do schedule(static) private(e,m,a,f)
      Do f = 1, v
      Do e = 1, v
      Do m = 1, o
      Do a = 1, v
        Vvovv(a,m,e,f) = H2(o+a,m,o+e,o+f)
      End Do
      End Do
      End Do
      End Do
      !$omp end parallel do
      !$omp parallel do schedule(static) private(j,e,b,a)
      Do j = 1, o
      Do e = 1, v
      Do b = 1, v
      Do a = 1, v
        Vvvvo(a,b,e,j) = H2(o+a,o+b,o+e,j)
      End Do
      End Do
      End Do
      End Do
      !$omp end parallel do
      !$omp parallel do schedule(static) private(e,m,i,f)
      Do e = 1, v
      Do i = 1, o
      Do m = 1, o
      Do f = 1, v
        Vfmne(f,m,i,e) = Voovv(m,i,e,f)
        Vmenf(m,e,i,f) = Voovv(m,i,e,f)
      End Do
      End Do
      End Do
      End Do
      !$omp end parallel do
      Do i = 1, o
      Do m = 1, o
      Do k = 1, o
      Do f = 1, v
        Vfmni(f,m,k,i) = Vooov(m,k,i,f)
      End Do
      End Do
      End Do
      End Do
      If (DoT) Then
        Allocate(Tvvvo(v,o,v,v), Tovoo(o,o,o,v), Stat=IAlloc)
        If (IAlloc /= 0) Stop "Could not allocate (T) blocks in SetUpERIBlocks"
        !$omp parallel do schedule(static) private(b,k,a,e)
        Do e = 1, v
        Do b = 1, v
        Do k = 1, o
        Do a = 1, v
          Tvvvo(a,k,b,e) = H2(o+b,o+e,o+a,k)
        End Do
        End Do
        End Do
        End Do
        !$omp end parallel do
        Do a = 1, v
        Do k = 1, o
        Do j = 1, o
        Do m = 1, o
          Tovoo(m,j,k,a) = H2(m,o+a,j,k)
        End Do
        End Do
        End Do
        End Do
        ERIBlocksT = .True.
      End If
      If (DoC) Then
        Voooo = Conjg(Voooo); Vooov = Conjg(Vooov); Voovv = Conjg(Voovv); Vovvo = Conjg(Vovvo)
        Vovoo = Conjg(Vovoo); Vvovv = Conjg(Vvovv); Vvvvo = Conjg(Vvvvo); Vvvoo = Conjg(Vvvoo)
        Vfmne = Conjg(Vfmne); Vfmni = Conjg(Vfmni); Vmenf = Conjg(Vmenf)
        If (ERIBlocksT) Then
          Tvvvo = Conjg(Tvvvo); Tovoo = Conjg(Tovoo)
        End If
      End If
      ERIBlocksConj = DoC
      ERIBlocksReady = .True.
      End Subroutine SetUpERIBlocks

      Subroutine ShutDownERIBlocks
      Implicit None
      If (.not. ERIBlocksReady) Return
      Deallocate(Voooo, Vooov, Voovv, Vovvo, Vovoo, Vvovv, Vvvvo, Vvvoo, Vfmne, Vfmni, Vmenf)
      If (ERIBlocksT) Deallocate(Tvvvo, Tovoo)
      ERIBlocksReady = .False.
      ERIBlocksT = .False.
      ERIBlocksConj = .False.
      End Subroutine ShutDownERIBlocks

   End Module ERIBlocks
