
   Module CCSD_T
!
!  (T) correction, same formula as the original (W_abc, V_abc per virtual
!  triple a>=b>=c, three permuted calls of get_wv), but:
!    * the integral blocks come pre-permuted from ERIBlocks (Tvvvo(d,k,b,c) =
!      <bc||dk>, Tovoo(l,j,k,a) = <la||jk>, Vvvoo) and T2 is permuted once,
!      so every ZGEMM/ZGERU operand is contiguous or strided by a constant
!      increment - no array-section temporaries per triple,
!    * the triangular a>=b>=c loop is flattened and distributed dynamically,
!    * no nested OpenMP inside the per-triple helpers.
!
   Use Precision
   Use Constants
   Use ERIBlocks
   Use MPI
   Implicit None

   Contains

      Subroutine CCSD_T_Ene(F, MOE, T1, T2, NOcc, NVir, NSO, Ene, comm)
!     comm (optional): the a>=b>=c triples are dealt round-robin over its ranks and the energy
!     is summed with MPI_Allreduce, so every rank returns the full E(T).  Without comm the
!     routine is serial in MPI (used inside a grid loop that is already distributed).
!     MOE(NSO): orbital energies for the triples denominators (semicanonical
!     Fock diagonal, or the PHF effective-Fock eigenvalues).  F supplies f_ov.
      Implicit None
      Integer, Intent(In) :: NSO, NOcc, NVir
      Real (Kind=pr), Intent(In) :: MOE(NSO)
      Complex (Kind=pr), Intent(In) :: T1(nvir,nocc), T2(nvir,nvir,nocc,nocc)
      Complex (Kind=pr), Intent(In) :: F(NSO,NSO)
      Complex (Kind=pr), Intent(Out) :: Ene
      Integer, Intent(In), Optional :: comm
      Integer :: myrank, nrank, ierr, mpi_dp
      Complex (Kind=pr) :: EneLoc
      Complex (Kind=pr), Allocatable :: T2p(:,:,:,:), T2q(:,:,:,:), fvo(:,:), eijk(:,:,:)
      Complex (Kind=pr), Allocatable :: wabc(:,:,:), vabc(:,:,:), wcab(:,:,:), vcab(:,:,:)
      Complex (Kind=pr), Allocatable :: wbac(:,:,:), vbac(:,:,:), w(:,:,:), v(:,:,:)
      Complex (Kind=pr) :: eabc, Et
      Integer :: a, b, c, d, i, j, k, l, it, ntrip, IAlloc
      Integer, Allocatable :: ta(:), tb(:), tc(:)
      If (.not. ERIBlocksReady .or. .not. ERIBlocksT) Stop "CCSD_T_Ene: (T) integral blocks not set up"
      Allocate(T2p(nvir,nocc,nocc,nvir), T2q(nocc,nocc,nvir,nvir), fvo(nvir,nocc), &
               eijk(nocc,nocc,nocc), Stat=IAlloc)
      If (IAlloc /= 0) Stop "Could not allocate in CCSD_T_Ene"
      fvo = F(nocc+1:nso, 1:nocc)
      !$omp parallel do private(a,d,i,j)
      Do a = 1, nvir
      Do j = 1, nocc
      Do i = 1, nocc
      Do d = 1, nvir
        T2p(d,i,j,a) = T2(a,d,i,j)       ! T2p(:,:,:,a) : (d | i,j) contiguous
      End Do
      End Do
      End Do
      End Do
      !$omp end parallel do
      !$omp parallel do private(b,c,i,l)
      Do c = 1, nvir
      Do b = 1, nvir
      Do l = 1, nocc
      Do i = 1, nocc
        T2q(i,l,b,c) = T2(b,c,i,l)       ! T2q(:,:,b,c) : (i | l) contiguous
      End Do
      End Do
      End Do
      End Do
      !$omp end parallel do
      Do k = 1, nocc
      Do j = 1, nocc
      Do i = 1, nocc
        eijk(i,j,k) = MOE(i) + MOE(j) + MOE(k)
      End Do
      End Do
      End Do
! flatten the triangle a >= b >= c
      ntrip = nvir*(nvir+1)*(nvir+2)/6
      Allocate(ta(ntrip), tb(ntrip), tc(ntrip))
      it = 0
      Do a = 1, nvir
      Do b = 1, a
      Do c = 1, b
        it = it + 1
        ta(it) = a; tb(it) = b; tc(it) = c
      End Do
      End Do
      End Do
      myrank = 0; nrank = 1
      If (Present(comm)) Then
        Call MPI_Comm_rank(comm, myrank, ierr); Call MPI_Comm_size(comm, nrank, ierr)
      End If
      Ene = Zero
      !$omp parallel default(shared) private(a,b,c,it,eabc,Et,w,wabc,wcab,wbac,v,vabc,vcab,vbac)
      Allocate(wabc(nocc,nocc,nocc), vabc(nocc,nocc,nocc), wcab(nocc,nocc,nocc), vcab(nocc,nocc,nocc), &
               wbac(nocc,nocc,nocc), vbac(nocc,nocc,nocc), w(nocc,nocc,nocc), v(nocc,nocc,nocc))
      Et = Zero
      !$omp do schedule(dynamic,4)
      Do it = 1 + myrank, ntrip, nrank
        a = ta(it); b = tb(it); c = tc(it)
        Call get_wv(wabc,vabc,T1,T2p,T2q,a,b,c,fvo,F,nocc,nvir,nso)
        Call get_wv(wcab,vcab,T1,T2p,T2q,c,a,b,fvo,F,nocc,nvir,nso)
        Call get_wv(wbac,vbac,T1,T2p,T2q,b,a,c,fvo,F,nocc,nvir,nso)
        w = wabc + wcab - wbac
        v = vabc + vcab - vbac
        eabc = MOE(nocc+a) + MOE(nocc+b) + MOE(nocc+c)
        w = w / (eijk - eabc)
        Et = Et + sum(w * conjg(v))
      End Do
      !$omp end do
      !$omp critical
      Ene = Ene + Et
      !$omp end critical
      Deallocate(wabc, vabc, wcab, vcab, wbac, vbac, w, v)
      !$omp end parallel
      If (nrank > 1) Then
        Call MPI_TYPE_CREATE_F90_COMPLEX(15, 307, mpi_dp, ierr)
        EneLoc = Ene
        Call MPI_Allreduce(EneLoc, Ene, 1, mpi_dp, MPI_Sum, comm, ierr)
      End If
      Ene = Ene / Two
      Deallocate(T2p, T2q, fvo, eijk, ta, tb, tc)
      End Subroutine CCSD_T_Ene


      Subroutine get_wv(w,v,T1,T2p,T2q,a,b,c,fvo,F,nocc,nvir,nso)
!     w(i,j,k) = P(ijk)[ sum_d <bc||dk> t_ij^ad - sum_l t_il^bc <la||jk> ]
!     v(i,j,k) = w(before P) + t_a^i <bc||jk> + f_ai t_jk^bc
      Implicit None
	  Integer, Intent(In) :: nocc, nvir, nso, a, b, c
      Complex (Kind=pr), Intent(In) :: T1(nvir,nocc), T2p(nvir,nocc,nocc,nvir), T2q(nocc,nocc,nvir,nvir)
	  Complex (Kind=pr), Intent(In) :: fvo(nvir,nocc), F(nso,nso)
	  Complex (Kind=pr), Intent(Out) :: w(nocc,nocc,nocc), v(nocc,nocc,nocc)
      Integer :: no2
      no2 = nocc*nocc
	  w = Zero
      ! w(k,ij) += sum_d Tvvvo(d,k,b,c) T2p(d,i,j,a)      [A^T B,  A = Tvvvo(:,:,b,c) (d x k)]
      Call ZGEMM('T','N',nocc,no2,nvir,ZOne,Tvvvo(1,1,b,c),nvir,T2p(1,1,1,a),nvir,ZOne,w,nocc)
      ! w(i,jk) -= sum_l T2q(i,l,b,c) Tovoo(l,j,k,a)
      Call ZGEMM('N','N',nocc,no2,nocc,-ZOne,T2q(1,1,b,c),nocc,Tovoo(1,1,1,a),nocc,ZOne,w,nocc)
	  v = Zero
      ! v(i,jk) += t1(a,i) Vvvoo(b,c,j,k)   (strided vectors, no copies)
      Call ZGERU(nocc,no2,ZOne,T1(a,1),nvir,Vvvoo(b,c,1,1),nvir*nvir,v,nocc)
      ! v(i,jk) += fvo(a,i) T2q(j,k,b,c)
      Call ZGERU(nocc,no2,ZOne,fvo(a,1),nvir,T2q(1,1,b,c),1,v,nocc)
	  v = v + w
	  Call Perm3(w, nocc)
      End Subroutine get_wv

	  Subroutine Perm3(w, n)
	  Implicit None
	  Integer, Intent(In) :: n
	  Complex (Kind=pr), Intent(InOut) :: w(n,n,n)
	  Complex (Kind=pr) :: tmp(n,n,n)
	  Integer :: i,j,k
	  tmp = w
	  Do k = 1, n
	  Do j = 1, n
	  Do i = 1, n
	  	  w(i,j,k) = tmp(i,j,k) + tmp(k,i,j) + tmp(j,k,i)
	  End Do
	  End Do
	  End Do
	  End Subroutine Perm3

   End Module CCSD_T
