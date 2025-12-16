    Module PHFTools
    Use Precision
    Use Constants
    Use BuildSC_dr
	Use mpi
    Implicit None

    Contains
        
        Subroutine BuildHSG(HOne,HTwo,Ref,Z,Ngrid,Npg,nci,CmplxConj,ncisp,ncipg,NOcc,NSO, &
                            R1,R2,Rpg,weightsp,weightpg,roota,rootb,rooty,QJ,SP, &
                            Hmat,Smat,Gradmat,Rdmmat,comm)
        Implicit None
! PHF 
        Integer,            Intent(in)  :: NOcc,NSO,Ngrid(2),Npg,QJ,SP
        Integer,            Intent(in)  :: nci,CmplxConj,ncisp,ncipg
        Complex(kind=pr),   Intent(in)  :: HOne(NSO,NSO), HTwo(NSO,NSO,NSO,NSO)
        Complex(kind=pr),   Intent(in)  :: Ref(NSO,NSO),Z(NSO-NOcc,NOcc)
        Complex(kind=pr),   Intent(in)  :: R1(:,:,:), R2(:,:,:)
         Complex(kind=pr),  Intent(in)  :: Rpg(npg,NSO,NSO)
        Real    (Kind=pr),  Intent(In)  :: Weightsp(:,:),Weightpg(Npg,ncipg,ncipg)
        Real    (Kind=pr),  Intent(In)  :: Roota(:), Rootb(:), Rooty(:)
        Integer :: I, J, ni, nj, p ,q, m, n, d, il, iy, ipg
        Complex(kind=pr), Intent(Out) :: Hmat(nci,nci), Smat(nci,nci)
        Complex(kind=pr), Intent(Out) :: Gradmat(nci,nci,NSO,NSO)
        Complex(kind=pr), Intent(Out) :: Rdmmat(nci,nci,NSO,NSO)
        Complex(kind=pr) :: newMO(NSO,NOcc), R(NSO,NSO), Rtmp(NSO,NSO)
        Complex(kind=pr) :: rho0(NSO,NSO), s0, h0, g0(NSO,NSO)
        Complex(kind=pr) :: rho1(NSO,NSO), s1, h1, g1(NSO,NSO)
        Complex(kind=pr) :: rho0c(NSO,NSO), s0c, g0c(NSO,NSO)
        Complex(kind=pr) :: rho1c(NSO,NSO), s1c, g1c(NSO,NSO)
        Complex(kind=pr) :: gama0(NSO,NSO), gama1(NSO,NSO)
        Complex(kind=pr) :: gama0c(NSO,NSO), gama1c(NSO,NSO)
		Complex(kind=pr) :: Ztmp(NSO-NOcc,NOcc)
        Real(kind=pr) :: w1, w2, Tol = 1e-5
        Complex(kind=pr) :: norm, wig, w
        Integer         :: cnt1, cnt2, clock_rate, clock_max
! openmp
        INTEGER, EXTERNAL :: OMP_GET_THREAD_NUM, OMP_GET_NUM_THREADS
! mpi
		integer, intent(in) :: comm
		integer :: rank, ntask, ierr, stat, n1, n2, mpi_dp
		integer, allocatable :: workload(:)
		integer :: istart, iend
        integer, allocatable :: NIter(:)
		Complex(kind=pr) :: iHmat(nci,nci), iSmat(nci,nci)
		Complex(kind=pr) :: iGradmat(nci,nci,NSO,NSO)
		Complex(kind=pr) :: iRdmmat(nci,nci,NSO,NSO)
!		Call MPI_INIT ( ierr )
		Call MPI_Comm_size(comm, ntask, ierr)
		Call MPI_Comm_rank(comm, rank, ierr)
		Call MPI_TYPE_CREATE_F90_COMPLEX(15, 307, mpi_dp, ierr)
		Allocate(workload(ntask))
		workload = Ngrid(2) / ntask
		workload(1:mod(Ngrid(2),ntask)) = workload(1:mod(Ngrid(2),ntask)) + 1
		istart = sum(workload(1:rank))
		iend = istart + workload(rank+1)
		iHmat = Zero
		iSmat = Zero
		iGradmat = Zero
		iRdmmat = Zero
		n1 = NSO * NSO
		n2 = nci * nci
! initialize PHF
        Hmat = Zero
        Smat = Zero
        Gradmat = Zero
        Rdmmat = Zero
        newMO = Zero
		Ztmp = Z
        If (CmplxConj == 0) Then
            Ztmp = Cmplx(Real(Z,kind=pr),Zero,kind=pr)
        End If
        Do I= 1, NOcc
            newMO(I,I) = One
        End Do
        newMO(NOcc+1:NSO,1:NOcc) = Ztmp
! start PHF
		!$omp parallel default(shared)

		!$omp do schedule(static) reduction(+:iHmat,iSmat,iGradmat,iRdmmat), &
		!$omp& private(Rtmp,R,rho0,s0,h0,g0,rho1,s1,h1,g1,rho0c,s0c,g0c,rho1c,s1c,g1c), &
		!$omp& private(gama0,gama1,gama0c,gama1c,wig,w,ni,nj,p,q,m,n,il,iy,ipg)
		Do il = 1, Ngrid(1)
		Do iy = 1, Ngrid(2)
!        Do iy = 1+istart, iend
		Do ipg = 1, Npg
            Rtmp = MatMul(R1(il,:,:), R2(iy,:,:))
            Rtmp = MatMul(Rtmp, Rpg(ipg,:,:))
            Call ao2mo2(R,Rtmp,Ref,NSO)
            Call RDM(R,newMO,NOcc,NSO,rho0,s0)
            Call MakeGamma(HTwo,rho0,NSO,gama0)
            Call Kernels(HOne,HTwo,gama0,rho0,s0,NSO,h0)
            Call Gradient(HOne,HTwo,gama0,rho0,s0,NSO,g0)
            If (CmplxConj == 2) Then
                Call RDMConj(R,newMO,NOcc,NSO,rho1,s1)
                Call MakeGamma(HTwo,rho1,NSO,gama1)
                Call Kernels(HOne,HTwo,gama1,rho1,s1,NSO,h1)
                Call Gradient(HOne,HTwo,gama1,rho1,s1,NSO,g1)

                Call RDM(R,Conjg(newMO),NOcc,NSO,rho0c,s0c)
                Call MakeGamma(HTwo,rho0c,NSO,gama0c)
                Call Gradient(HOne,HTwo,gama0c,rho0c,s0c,NSO,g0c)

                Call RDMConj(R,Conjg(newMO),NOcc,NSO,rho1c,s1c)
                Call MakeGamma(HTwo,rho1c,NSO,gama1c)
                Call Gradient(HOne,HTwo,gama1c,rho1c,s1c,NSO,g1c)
            Else
                h1 = h0
                g1 = g0
                s1 = s0
                rho1 = rho0
                g0c = g0
                s0c = s0
                rho0c = rho0
                g1c = g0
                s1c = s0
                rho1c = rho0
            End If
            rho0 = rho0 * weightsp(il,iy) * s0
            rho1 = rho1 * weightsp(il,iy) * s1
            rho0c = rho0c * weightsp(il,iy) * s0c
            rho1c = rho1c * weightsp(il,iy) * s1c

            h0 = weightsp(il,iy) * h0
            h1 = weightsp(il,iy) * h1

            s0 = weightsp(il,iy) * s0
            s1 = weightsp(il,iy) * s1
            s0c = s0c * weightsp(il,iy)
            s1c = s1c * weightsp(il,iy)

            g0 = g0 * weightsp(il,iy)
            g1 = g1 * weightsp(il,iy)
            g0c = g0c * weightsp(il,iy)
            g1c = g1c * weightsp(il,iy)
            Do ni = -QJ, QJ
            Do nj = -QJ, QJ
                if (SP == 2) then
                    Call Wignerfac(QJ,ni,nj,roota(il),rootb(iy),rooty(il),wig)
                else
                    Call Wignerfac(QJ,ni,nj,roota(il),rootb(il),rooty(iy),wig)
                end if
                Do p = 1, ncipg
                Do q = 1, ncipg
                    w = weightpg(ipg,p,q) * wig
                    m = ni + 1 + QJ + (p-1) * ncisp
                    n = nj + 1 + QJ + (q-1) * ncisp
                    iHmat(m,n) = iHmat(m,n) + h0 * w
                    iSmat(m,n) = iSmat(m,n) + s0 * w
                    iGradmat(m,n,:,:) = iGradmat(m,n,:,:) + g0 * w
                    iRdmmat(m,n,:,:) = iRdmmat(m,n,:,:) + rho0 * w
                    If (CmplxConj == 2) Then
                        d = ncisp * ncipg
                        iHmat(m,n+d) = iHmat(m,n+d) + h1 * w
                        iSmat(m,n+d) = iSmat(m,n+d) + s1 * w
                        iGradmat(m,n+d,:,:) = iGradmat(m,n+d,:,:) + g1 * w
                        iRdmmat(m,n+d,:,:) = iRdmmat(m,n+d,:,:) + rho1 * w
                        iGradmat(m+d,n,:,:) = iGradmat(m+d,n,:,:) + g1c * w
                        iRdmmat(m+d,n,:,:) = iRdmmat(m+d,n,:,:) + rho1c * w
                        iGradmat(m+d,n+d,:,:) = iGradmat(m+d,n+d,:,:) + g0c * w
                        iRdmmat(m+d,n+d,:,:) = iRdmmat(m+d,n+d,:,:) + rho0c * w
                    End If
                End Do
                End Do
            End Do
            End Do
        End Do
        End Do
        End Do
		!$omp end do

		!$omp end parallel
        If (CmplxConj == 2) Then
            iHmat(d+1:,1:d) = Conjg(Transpose(iHmat(1:d,d+1:)))
            iSmat(d+1:,1:d) = Conjg(Transpose(iSmat(1:d,d+1:)))
            iHmat(d+1:,d+1:) = iHmat(1:d,1:d)
            iSmat(d+1:,d+1:) = iSmat(1:d,1:d)
        End If
! collect result
		Call MPI_reduce(iHmat,Hmat,n2,mpi_dp,MPI_Sum,0,comm,ierr)
		Call MPI_reduce(iSmat,Smat,n2,mpi_dp,MPI_Sum,0,comm,ierr)
		Call MPI_reduce(iGradmat,Gradmat,n1*n2,mpi_dp,MPI_Sum,0,comm,ierr)
		Call MPI_reduce(iRDMmat,RDMmat,n1*n2,mpi_dp,MPI_Sum,0,comm,ierr)
		Call MPI_Bcast(Hmat,n2,mpi_dp,0,comm,ierr)
		Call MPI_Bcast(Smat,n2,mpi_dp,0,comm,ierr)
		Call MPI_Bcast(Gradmat,n1*n2,mpi_dp,0,comm, ierr)
		Call MPI_Bcast(RDMmat,n1*n2,mpi_dp,0,comm, ierr)
!		Call MPI_FINALIZE ( ierr )
!		Hmat = iHmat
!		Smat = iSmat
!		Gradmat = iGradmat
!		RDMmat = iRDMmat
!		Print *, "HImag",maxval(abs(aimag(Hmat))), maxval(abs(aimag(Smat))),maxval(abs(real(Ztmp)))
        w1 = maxval( abs( DImag((Conjg(Transpose(Hmat)) - Hmat))) )
        w2 = maxval( abs( DImag((Conjg(Transpose(Smat)) - Smat))) )
        If ((w1 > maxval(abs(Hmat)) * tol) .or. (w2 > maxval(abs(Smat)) * tol)) then
            Print *, maxval(abs(Ztmp))
            Stop "Not Hermitian"
        End If
        End Subroutine BuildHSG

        Subroutine LocalG(Gradmat,Rdmmat,Smat,E0,NSO,fsp,fpg,fk,nci,ncisp,ncipg,CmplxConj,G)
        Implicit None
        Integer,            Intent(in)  :: NSO, nci, ncisp, ncipg, CmplxConj
        Complex(kind=pr),   Intent(in)  :: Gradmat(nci,nci,NSO,NSO)
        Complex(kind=pr),   Intent(in)  :: Rdmmat(nci,nci,NSO,NSO), Smat(nci,nci)
        Complex(kind=pr),   Intent(in)  :: fsp(ncisp),fpg(ncipg),fk(:),E0
        Complex(kind=pr),   Intent(out) :: G(NSO,NSO)
        Integer :: ni, nj, p ,q, m, n, d
        Complex(kind=pr) :: Rdm(NSO,NSO), Ovlp, fac
        G = Zero
        Rdm = Zero
        Ovlp = Zero
        Do ni = 1, ncisp
        Do nj = 1, ncisp
        Do p = 1, ncipg
        Do q = 1, ncipg
            m = ni + (p-1) * ncisp
            n = nj + (q-1) * ncisp
            G = G + Conjg(fsp(ni)) * fsp(nj) * Conjg(fpg(p)) * fpg(q) * &
            Conjg(fk(1)) * fk(1) * Gradmat(m,n,:,:)
            Rdm = Rdm + Conjg(fsp(ni)) * fsp(nj) * Conjg(fpg(p)) * fpg(q) * &
            Conjg(fk(1)) * fk(1) * Rdmmat(m,n,:,:)
            Ovlp = Ovlp + Conjg(fsp(ni)) * fsp(nj) * Conjg(fpg(p)) * fpg(q) * &
            Conjg(fk(1)) * fk(1) * Smat(m,n)
            If (CmplxConj == 2) Then
                d = ncisp * ncipg
                fac = Conjg(fsp(ni)) * fsp(nj) * Conjg(fpg(p)) * fpg(q)
                G = G + fac * Conjg(fk(2)) * fk(2) * Conjg(Gradmat(m+d,n+d,:,:))
                Rdm = Rdm + fac * Conjg(fk(2)) * fk(2) * Conjg(Rdmmat(m+d,n+d,:,:))
                Ovlp = Ovlp + fac * Conjg(fk(2)) * fk(2) * Smat(m+d,n+d)
                                                                                    
                G = G + fac * Conjg(fk(1)) * fk(2) * Gradmat(m,n+d,:,:)
                Rdm = Rdm + fac * Conjg(fk(1)) * fk(2) * Rdmmat(m,n+d,:,:)
                Ovlp = Ovlp + fac * Conjg(fk(1)) * fk(2) * Smat(m,n+d)
                                                                                    
                G = G + fac * Conjg(fk(1)) * fk(2) * Conjg(Gradmat(m+d,n,:,:))
                Rdm = Rdm + fac * Conjg(fk(1)) * fk(2) * Conjg(Rdmmat(m+d,n,:,:))
                Ovlp = Ovlp + fac * Conjg(fk(2)) * fk(1) * Smat(m+d,n)
            End If
        End Do
        End Do
        End Do
        End Do
        G = (G - E0 * Rdm) / Ovlp
        End Subroutine LocalG

        Subroutine RDM(R,MO,NOcc,NSO,rho,det)
        Implicit None
        Integer,            Intent(in)  :: NOcc, NSO
        Complex(kind=pr),   Intent(in)  :: R(NSO,NSO)
        Complex(kind=pr),   Intent(in)  :: MO(NSO,NOcc)
        Complex(kind=pr),   Intent(Out) :: det, rho(NSO,NSO)
        Complex(kind=pr) :: M(NOcc,NOcc), Minv(NOcc,NOcc)
		 Complex(kind=pr) :: tmp(NSO,NOcc)
		tmp = MatMul(R, MO)
        M = MatMul(Conjg(Transpose(MO)), tmp)
        Call InvertC(M,Minv,NOcc)
        Call determinantC(det,M,NOcc)
        tmp = MatMul(tmp, Minv)
        rho = MatMul(tmp, Conjg(Transpose(MO)))
        End Subroutine RDM

        Subroutine RDMConj(R,MO,NOcc,NSO,rho,det)
        Implicit None
        Integer,            Intent(in)  :: NOcc, NSO
        Complex(kind=pr),   Intent(in)  :: R(NSO,NSO)
        Complex(kind=pr),   Intent(in)  :: MO(NSO,NOcc)
        Complex(kind=pr),   Intent(Out) :: det, rho(NSO,NSO)
        Complex(kind=pr) :: M(NOcc,NOcc), Minv(NOcc,NOcc)
	    Complex(kind=pr) :: tmp(NSO,NOcc)
		tmp = MatMul(R, Conjg(MO))
        M = MatMul(Conjg(Transpose(MO)), tmp)
        Call InvertC(M,Minv,NOcc)
        Call determinantC(det,M,NOcc)
        tmp = MatMul(tmp, Minv)
        rho = MatMul(tmp, Conjg(Transpose(MO)))
        End Subroutine RDMConj

        Subroutine Kernels(H1,H2,gama,rho,det,NSO,h)
        Implicit None
        Integer,            Intent(in)  :: NSO
        Complex(kind=pr),   Intent(in)  :: H1(NSO,NSO), H2(NSO,NSO,NSO,NSO)
        Complex(kind=pr),   Intent(in)  :: gama(NSO,NSO), rho(NSO,NSO), det
        Complex(kind=pr),   Intent(out) :: h
        Call Makeh(H1,H2,gama,rho,NSO,h)
        h = h * det
        End Subroutine Kernels

        Subroutine Gradient(H1,H2,gama,rho,det,NSO,G)
        Implicit None
        Integer,            Intent(in)  :: NSO
        Complex(kind=pr),   Intent(in)  :: H1(NSO,NSO), H2(NSO,NSO,NSO,NSO)
        Complex(kind=pr),   Intent(in)  :: gama(NSO,NSO), rho(NSO,NSO), det
        Complex(kind=pr),   Intent(out) :: G(NSO,NSO)
        Complex(kind=pr) :: h, Id(NSO,NSO)
        Integer :: I
        Id = Zero
        Do I = 1, NSO
            Id(I,I) = One
        End Do
        Call Makeh(H1,H2,gama,rho,NSO,h)
        G = MatMul(H1+gama,rho)
        G = MatMul(Id-rho,G)
        G = G + h * rho
        G = G * det
        End Subroutine Gradient

        Subroutine MakeGamma(H2,rho,NSO,gama)
        Implicit None
        Integer,            Intent(in)  :: NSO
        Complex(kind=pr),   Intent(in)  :: H2(NSO,NSO,NSO,NSO)
        Complex(kind=pr),   Intent(in)  :: rho(NSO,NSO)
        Complex(kind=pr),   Intent(out) :: gama(NSO,NSO)
        Integer :: i,j,k,l
        integer, allocatable :: workload(:), istart(:), iend(:)
        integer :: m, nblk
        Integer         :: cnt1, cnt2, clock_rate, clock_max
        nblk = NSO / 4
        Allocate(workload(nblk), istart(nblk), iend(nblk))
        workload = nso / nblk
        workload(1:mod(nso,nblk)) = workload(1:mod(nso,nblk)) + 1
        do m = 1, nblk
            istart(m) = sum(workload(1:m-1))
            iend(m) = istart(m) + workload(m)
        end do
        !$omp parallel do
        do m = 1, nblk
        Do i = 1,NSO
        Do k = istart(m)+ 1, iend(m)
            gama(i,k) =  sum(H2(i,:,k,:) * transpose(rho))
        End Do
        End Do
        end do
        !$omp end parallel do
        End Subroutine MakeGamma

        Subroutine Makeh(H1,H2,gama,rho,NSO,h)
        Implicit None
        Integer,            Intent(in)  :: NSO
        Complex(kind=pr),   Intent(in)  :: H1(NSO,NSO), H2(NSO,NSO,NSO,NSO)
        Complex(kind=pr),   Intent(in)  :: gama(NSO,NSO), rho(NSO,NSO)
        Complex(kind=pr),   Intent(out) :: h
        !h = Trace(MatMul(H1,rho),NSO)
        h = sum(transpose(H1) * rho)
        h = h + sum(gama * transpose(rho)) / 2
        End Subroutine Makeh

        Subroutine Wignerfac(J,N,M,alpha,beta,gama,wig)
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
        Djnm = Zero
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
        End Subroutine Wignerfac

        Function Trace(M,N)
        Implicit None
        Integer,            Intent(in)  :: N
        Complex(kind=pr),   Intent(in)  :: M(N,N)
        Complex(kind=pr) :: Trace
        Integer :: I
        Trace = Zero
        Do I = 1, N
            Trace = Trace + M(I,I)
        End Do
        End Function Trace

        Subroutine ao2mo2(newH,H,M,N)
        Implicit None
        Integer,            Intent(in)  :: N
        Complex(kind=pr),   Intent(in)  :: M(N,N)
        Complex(kind=pr),   Intent(in)  :: H(N,N)
        Complex(kind=pr),   Intent(out) :: newH(N,N)
        newH = MatMul(H,M)
        newH = MatMul(Conjg(Transpose(M)), newH)
        End Subroutine ao2mo2

        Subroutine ao2mo4(newH,H,M,N)
        Implicit None
        Integer,            Intent(in)  :: N
        Complex(kind=pr),   Intent(in)  :: M(N,N)
        Complex(kind=pr),   Intent(in)  :: H(N,N,N,N)
        Complex(kind=pr),   Intent(out) :: newH(N,N,N,N)
        Integer :: p,q,r,s, i,j,k,l
        Complex(kind=pr) :: tmp(N,N,N,N)
		newH = H
		
		!$omp parallel default(shared)

		!$omp do schedule(static)
		Do p = 1,N
		Do q = 1,N
		Do r = 1,N
		Do l = 1,N
			tmp(p,q,r,l) = Dot_Product(Conjg(M(:,l)), H(p,q,r,:))
		End Do
		End Do
		End Do
		End Do
		!$omp end do

		!$omp do schedule(static)
		Do p = 1,N
		Do q = 1,N
		Do k = 1,N
		Do l = 1,N
			newH(p,q,k,l) = Dot_Product(Conjg(M(:,k)), tmp(p,q,:,l))
		End Do
		End Do
		End Do
		End Do
		!$omp end do

		!$omp do schedule(static)
		Do p = 1,N
		Do j = 1,N
		Do k = 1,N
		Do l = 1,N
			tmp(p,j,k,l) = Dot_Product(M(:,j), newH(p,:,k,l))
		End Do
		End Do
		End Do
		End Do
		!$omp end do

		!$omp do schedule(static)
		Do i = 1,N
		Do j = 1,N
		Do k = 1,N
		Do l = 1,N
			newH(i,j,k,l) = Dot_Product(M(:,i), tmp(:,j,k,l))
		End Do
		End Do
		End Do
		End Do
		!$omp end do

		!$omp end parallel

        End Subroutine ao2mo4

        Subroutine BuildSC(t0,t1,t2,c0,c1,c2,NOcc,NSO,Ngrid,Npg,QJ, &
                            ncisp,ncipg,ncik,R1, R2, Rpg, Rk, &
                            Roota,Rootb,Rooty,Weightsp,Weightpg,fsp,fpg,fk)
        Implicit None
        Integer,            Intent(in)  :: NOcc,NSO,Ngrid(2),Npg,QJ
        Integer,            Intent(in)  :: ncisp,ncipg,ncik
        Real    (Kind=pr),  Intent(In)  :: Weightsp(:,:),Weightpg(Npg,ncipg,ncipg)
        Real    (Kind=pr),  Intent(In)  :: Roota(:), Rootb(:), Rooty(:)
        Complex(kind=pr),  Intent(in)  :: R1(:,:,:), R2(:,:,:)
        Complex(kind=pr),  Intent(in)  :: Rpg(npg,NSO,NSO), Rk(ncik,NSO,NSO)
        Complex(kind=pr),  Intent(in)  :: fsp(ncisp),fpg(ncipg),fk(ncik)
        Complex(kind=pr),  Intent(in)  :: c0, c1(NSO,NSO)
        Complex(kind=pr),  Intent(in)  :: c2(NSO,NSO,NSO,NSO)
        Complex(kind=pr),  Intent(out) :: t0, t1(NSO,NSO)
        Complex(kind=pr),  Intent(out) :: t2(NSO,NSO,NSO,NSO)
        Integer :: I, J, ni, nj, p ,q, m, n, d, il, iy, ipg, ik
        Complex(kind=pr) :: R(NSO,NSO), rho(NSO,NSO), ovlp, s, wig
        Complex(kind=pr) :: newMO(NSO,NOcc)
        Complex(kind=pr) :: a0, a1(NSO,NSO), a2(NSO,NSO,NSO,NSO)
        Complex(kind=pr) :: c0tmp, c1tmp(NSO,NSO), c2tmp(NSO,NSO,NSO,NSO)
        Real(kind=pr) :: w, w1, Tol = 1e-6
        t0 = Zero
        t1 = Zero
        t2 = Zero
		ovlp = Zero
        newMO = Zero
        Do I= 1, NOcc
            newMO(I,I) = One
        End Do

        !$omp parallel default(shared)
        !$omp do schedule(static) collapse(2) reduction(+:t0,t1,t2,ovlp), &
        !$omp& private(R,rho,s,c0tmp,c1tmp,c2tmp,a0,a1,a2), &
        !$omp& private(wig,w,m,n,p,q)
        Do il = 1, Ngrid(1)
        Do iy = 1, Ngrid(2)
        Do ipg = 1, Npg
            R = MatMul(R1(il,:,:), R2(iy,:,:))
            R = MatMul(R, Rpg(ipg,:,:))
            Call RDM(R,newMO,NOcc,NSO,rho,s)
            c0tmp = c0
            Call ao2mo2(c1tmp,c1,Conjg(Transpose(R)),NSO)
            Call ao2mo4(c2tmp,c2,Conjg(Transpose(R)),NSO)
            Call BuildSC_Vec(c0tmp,c1tmp,c2tmp,a0,a1,a2,rho,NSO)
            s = weightsp(il,iy) * s
            a0 = a0 * s 
            a1 = a1 * s 
            a2 = a2 * s
            Do m = -QJ, QJ
            Do n = -QJ, QJ
                Call Wignerfac(QJ,m,n,roota(il),rootb(il),rooty(iy),wig)
				wig = One
                Do p = 1, ncipg
                Do q = 1, ncipg
                    w = weightpg(ipg,p,q) * wig
                    w = w * conjg(fsp(m+1+QJ)) * fsp(n+1+QJ)
                    w = w * conjg(fpg(p)) * fpg(q)
                    t0 = t0 + a0 * w
                    t1 = t1 + a1 * w
                    t2 = t2 + a2 * w
                    ovlp = ovlp + s * w
                End Do
                End Do
            End Do
            End Do
        End Do
        End Do
        End Do
        !$omp end do
        !$omp end parallel
        t0 = t0 / ovlp
        t1 = t1 / ovlp
        t2 = t2 / ovlp
        End Subroutine BuildSC

		Subroutine ExtraGradient(H1,H2,gama,rho,det,NSO,G)
		Implicit None
		Integer,            Intent(in)  :: NSO
		Complex(kind=pr),   Intent(in)  :: H1(NSO,NSO), H2(NSO,NSO,NSO,NSO)
		Complex(kind=pr),   Intent(in)  :: gama(NSO,NSO), rho(NSO,NSO), det
		Complex(kind=pr),   Intent(out) :: G(NSO,NSO)
		Complex(kind=pr) :: h
		Call Makeh(H1,H2,gama,rho,NSO,h)
		G = MatMul(rho, H1+gama)
		G = MatMul(G, rho)
		G = -G + h * rho
		G = G * det
		End Subroutine ExtraGradient

		Subroutine EffFock(HOne,HTwo,Ref,Ngrid,Npg,nci,CmplxConj,ncisp,ncipg,NOcc,NSO, &
						   R1,R2,Rpg,weightsp,weightpg,roota,rootb,rooty,QJ,SP, &
						   Hmat,Smat,Gradmat,Rdmmat,comm)
		Implicit None
		Integer,            Intent(in)  :: NOcc,NSO,Ngrid(2),Npg,QJ,SP
		Integer,            Intent(in)  :: nci,CmplxConj,ncisp,ncipg
		Complex(kind=pr),   Intent(in)  :: HOne(NSO,NSO), HTwo(NSO,NSO,NSO,NSO)
		Complex(kind=pr),   Intent(in)  :: Ref(NSO,NSO)
		Complex(kind=pr),   Intent(in)  :: R1(:,:,:), R2(:,:,:)
		Complex(kind=pr),   Intent(in)  :: Rpg(npg,NSO,NSO)
		Real    (Kind=pr),  Intent(In)  :: Weightsp(:,:),Weightpg(Npg,ncipg,ncipg)
		Real    (Kind=pr),  Intent(In)  :: Roota(:), Rootb(:), Rooty(:)
		Integer :: I, J, ni, nj, p ,q, m, n, d, il, iy, ipg
		integer, intent(in) :: comm
		integer :: rank, ntask, ierr, stat, n1, n2, mpi_dp
		Fock = Zero
		Do il = 1, Ngrid(1)
		Do iy = 1, Ngrid(2)
		Do ipg = 1, Npg
			Rtmp = MatMul(R1(il,:,:), R2(iy,:,:))
			Rtmp = MatMul(Rtmp, Rpg(ipg,:,:))
			Call ao2mo2(R,Rtmp,Ref,NSO)
			Call RDM(R,newMO,NOcc,NSO,rho0,s0)
			Call MakeGamma(HTwo,rho0,NSO,gama0)
			Call Kernels(HOne,HTwo,gama0,rho0,s0,NSO,h0)
			Call Gradient(HOne,HTwo,gama0,rho0,s0,NSO,g0)
			Call ExtraGradient(HOne,HTwo,gama0,rho0,s0,NSO,g1)

			rho0 = rho0 * weightsp(il,iy) * s0
			h0 = weightsp(il,iy) * h0
			s0 = weightsp(il,iy) * s0
			g0 = g0 * weightsp(il,iy)
			g1 = g1 * weightsp(il,iy)
			Do ni = -QJ, QJ
			Do nj = -QJ, QJ
				if (SP == 2) then
					Call Wignerfac(QJ,ni,nj,roota(il),rootb(iy),rooty(il),wig)
				else
					Call Wignerfac(QJ,ni,nj,roota(il),rootb(il),rooty(iy),wig)
				end if
				Do p = 1, ncipg
				Do q = 1, ncipg
					w = weightpg(ipg,p,q) * wig
					m = ni + 1 + QJ + (p-1) * ncisp
					n = nj + 1 + QJ + (q-1) * ncisp
				End Do
				End Do
			End Do
			End Do
		End Do
		End Do
		End Do

		End Subroutine EffFock

		Subroutine FockSolver()
		Implicit None
		End Subroutine FockSolver


    End Module PHFTools

