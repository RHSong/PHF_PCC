    Module FockTools
!
!   Fock-matrix formulation of PHF (PHFB_formula.pdf, Eqs. 29-33, PHF special case).
!
!   All quantities are expressed in the (orthonormal) basis of HOne/HTwo/R1/R2/Rpg;
!   the current orbitals MO are given in that same basis.  Nothing is sliced by
!   occupied/virtual index range: the projectors rho, 1-rho, rho(g), 1-rho(g)
!   enforce the block structure of Eqs. 30/32 in any basis.
!
!   Per grid point g the un-Hermitized integrand is
!       X(g) = H(g) rho(g) + (1-rho(g)) Gamma(g) rho(g)          [Eq. 30, before -E]
!            + 1/2 rho Gamma(g) rho(g)                            [Eq. 32b, F_-]
!            + 1/2 (1-rho(g)) Gamma(g) (1-rho)                    [Eq. 32a, F_+]
!   with rho = C C^dag, rho(g) = R C M^-1 C^dag, Gamma(g) = h + vbar.rho(g),
!   H(g) = Tr[(h + Gamma/2) rho(g)] and -sigma^dag(g) = 1 - rho(g).
!   BuildHSG accumulates  det(g) w(g) X(g)  and  det(g) w(g) rho(g)  over the
!   projection grid into CI-indexed matrices; BuildFock contracts with the CI
!   coefficients, subtracts E, adds the Hermitian conjugate and normalises.
!
    Use Precision
    Use Constants
    Implicit None

    Contains

        Subroutine BuildHSG(HOne,HTwo,MO,Ngrid,Npg,nci,CmplxConj,ncisp,ncipg,NOcc,NSO, &
                            R1,R2,Rpg,weightsp,weightpg,roota,rootb,rooty,QJ,SP, &
                            Hmat,Smat,Fockmat,Rdmmat)
        Implicit None
! PHF
        Integer,            Intent(in)  :: NOcc,NSO,Ngrid(2),Npg,QJ,SP
        Integer,            Intent(in)  :: nci,CmplxConj,ncisp,ncipg
        Complex(kind=pr),   Intent(in)  :: HOne(NSO,NSO), HTwo(NSO,NSO,NSO,NSO)
        Complex(kind=pr),   Intent(in)  :: MO(NSO,NSO)
        Complex(kind=pr),   Intent(in)  :: R1(:,:,:), R2(:,:,:)
        Complex(kind=pr),   Intent(in)  :: Rpg(npg,NSO,NSO)
        Real    (Kind=pr),  Intent(In)  :: Weightsp(:,:),Weightpg(Npg,ncipg,ncipg)
        Real    (Kind=pr),  Intent(In)  :: Roota(:), Rootb(:), Rooty(:)
        Complex(kind=pr), Intent(Out) :: Hmat(nci,nci), Smat(nci,nci)
        Complex(kind=pr), Intent(Out) :: Fockmat(nci,nci,NSO,NSO)
        Complex(kind=pr), Intent(Out) :: Rdmmat(nci,nci,NSO,NSO)
        Integer :: ni, nj, p ,q, m, n, d, il, iy, ipg
        Complex(kind=pr) :: newMO(NSO,NOcc), R(NSO,NSO)
        Complex(kind=pr) :: rho0(NSO,NSO), s0, h0, g0(NSO,NSO)
        Complex(kind=pr) :: rho1(NSO,NSO), s1, h1, g1(NSO,NSO)
        Complex(kind=pr) :: rho0c(NSO,NSO), s0c, g0c(NSO,NSO)
        Complex(kind=pr) :: rho1c(NSO,NSO), s1c, g1c(NSO,NSO)
        Complex(kind=pr) :: gama0(NSO,NSO), gama1(NSO,NSO)
        Complex(kind=pr) :: gama0c(NSO,NSO), gama1c(NSO,NSO)
        Complex(kind=pr) :: dm0(NSO,NSO)
        Real(kind=pr) :: w1, w2, Tol = 1e-5
        Complex(kind=pr) :: wig, w
! initialize PHF
        Hmat = Zero
        Smat = Zero
        Fockmat = Zero
        Rdmmat = Zero
        d = ncisp * ncipg
        newMO = MO(:, 1:NOcc)
        If (CmplxConj == 0) Then
            newMO = Cmplx(Real(newMO,kind=pr),Zero,kind=pr)
        End If
        dm0 = MatMul(newMO, Conjg(Transpose(newMO)))
! start PHF
        !$omp parallel default(shared)

        ! collapse(2): distribute all (Lebedev, gamma) grid points over the threads, not only the
        ! Lebedev index (which capped the useful thread count at Ngrid(1) and unbalanced the load)
        !$omp do schedule(dynamic) collapse(2) reduction(+:Hmat,Smat,Fockmat,Rdmmat), &
        !$omp& private(R,rho0,s0,h0,g0,rho1,s1,h1,g1,rho0c,s0c,g0c,rho1c,s1c,g1c), &
        !$omp& private(gama0,gama1,gama0c,gama1c,wig,w,ni,nj,p,q,m,n,il,iy,ipg)
        Do il = 1, Ngrid(1)
        Do iy = 1, Ngrid(2)
        Do ipg = 1, Npg
            R = MatMul(R1(il,:,:), R2(iy,:,:))
            R = MatMul(R, Rpg(ipg,:,:))
            Call RDM(R,newMO,NOcc,NSO,rho0,s0)
            Call MakeGamma(HTwo,rho0,NSO,gama0)
            Call Kernels(HOne,HTwo,gama0,rho0,s0,NSO,h0)
            Call FockKernel(HOne,gama0,rho0,dm0,s0,NSO,g0)
            If (CmplxConj == 2) Then
! NOTE: the complex-conjugation blocks are carried over from the gradient code
! and have NOT been re-derived for the Fock formulation.  Do not rely on them.
                Call RDMConj(R,newMO,NOcc,NSO,rho1,s1)
                Call MakeGamma(HTwo,rho1,NSO,gama1)
                Call Kernels(HOne,HTwo,gama1,rho1,s1,NSO,h1)
                Call FockKernel(HOne,gama1,rho1,dm0,s1,NSO,g1)

                Call RDM(R,Conjg(newMO),NOcc,NSO,rho0c,s0c)
                Call MakeGamma(HTwo,rho0c,NSO,gama0c)
                Call FockKernel(HOne,gama0c,rho0c,Conjg(dm0),s0c,NSO,g0c)

                Call RDMConj(R,Conjg(newMO),NOcc,NSO,rho1c,s1c)
                Call MakeGamma(HTwo,rho1c,NSO,gama1c)
                Call FockKernel(HOne,gama1c,rho1c,Conjg(dm0),s1c,NSO,g1c)
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
                    Hmat(m,n) = Hmat(m,n) + h0 * w
                    Smat(m,n) = Smat(m,n) + s0 * w
                    Fockmat(m,n,:,:) = Fockmat(m,n,:,:) + g0 * w
                    Rdmmat(m,n,:,:) = Rdmmat(m,n,:,:) + rho0 * w
                    If (CmplxConj == 2) Then
                        Hmat(m,n+d) = Hmat(m,n+d) + h1 * w
                        Smat(m,n+d) = Smat(m,n+d) + s1 * w
                        Fockmat(m,n+d,:,:) = Fockmat(m,n+d,:,:) + g1 * w
                        Rdmmat(m,n+d,:,:) = Rdmmat(m,n+d,:,:) + rho1 * w
                        Fockmat(m+d,n,:,:) = Fockmat(m+d,n,:,:) + g1c * w
                        Rdmmat(m+d,n,:,:) = Rdmmat(m+d,n,:,:) + rho1c * w
                        Fockmat(m+d,n+d,:,:) = Fockmat(m+d,n+d,:,:) + g0c * w
                        Rdmmat(m+d,n+d,:,:) = Rdmmat(m+d,n+d,:,:) + rho0c * w
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
            Hmat(d+1:,1:d) = Conjg(Transpose(Hmat(1:d,d+1:)))
            Smat(d+1:,1:d) = Conjg(Transpose(Smat(1:d,d+1:)))
            Hmat(d+1:,d+1:) = Hmat(1:d,1:d)
            Smat(d+1:,d+1:) = Smat(1:d,1:d)
        End If
! sanity check
        w1 = maxval( abs( DImag((Conjg(Transpose(Hmat)) - Hmat))) )
        w2 = maxval( abs( DImag((Conjg(Transpose(Smat)) - Smat))) )
        If ((w1 > maxval(abs(Hmat)) * tol) .or. (w2 > maxval(abs(Smat)) * tol)) then
            Stop "Not Hermitian"
        End If
        End Subroutine BuildHSG

        Subroutine BuildFock(Fockmat,Rdmmat,Smat,E0,NSO,fsp,fpg,fk,nci,ncisp,ncipg,CmplxConj,Fock)
!
!   Contract the CI-indexed integrands with the CI coefficients, subtract E0,
!   add the Hermitian conjugate (the "+h.c." of Eqs. 30/32) and normalise by
!   <Phi|P|Phi>.  Returns the PHF Fock matrix F = F0 + F+ + F- of Eq. 29.
!
        Implicit None
        Integer,            Intent(in)  :: NSO, nci, ncisp, ncipg, CmplxConj
        Complex(kind=pr),   Intent(in)  :: Fockmat(nci,nci,NSO,NSO)
        Complex(kind=pr),   Intent(in)  :: Rdmmat(nci,nci,NSO,NSO), Smat(nci,nci)
        Complex(kind=pr),   Intent(in)  :: fsp(ncisp),fpg(ncipg),fk(:),E0
        Complex(kind=pr),   Intent(out) :: Fock(NSO,NSO)
        Integer :: ni, nj, p ,q, m, n, d
        Complex(kind=pr) :: Rdm(NSO,NSO), Ovlp, fac
        Fock = Zero
        Rdm = Zero
        Ovlp = Zero
        d = ncisp * ncipg
        Do ni = 1, ncisp
        Do nj = 1, ncisp
        Do p = 1, ncipg
        Do q = 1, ncipg
            m = ni + (p-1) * ncisp
            n = nj + (q-1) * ncisp
            fac = Conjg(fsp(ni)) * fsp(nj) * Conjg(fpg(p)) * fpg(q)
            Fock = Fock + fac * Conjg(fk(1)) * fk(1) * Fockmat(m,n,:,:)
            Rdm = Rdm + fac * Conjg(fk(1)) * fk(1) * Rdmmat(m,n,:,:)
            Ovlp = Ovlp + fac * Conjg(fk(1)) * fk(1) * Smat(m,n)
            If (CmplxConj == 2) Then
! NOTE: carried over from LocalG, not re-derived for the Fock formulation.
                Fock = Fock + fac * Conjg(fk(2)) * fk(2) * Conjg(Fockmat(m+d,n+d,:,:))
                Rdm = Rdm + fac * Conjg(fk(2)) * fk(2) * Conjg(Rdmmat(m+d,n+d,:,:))
                Ovlp = Ovlp + fac * Conjg(fk(2)) * fk(2) * Smat(m+d,n+d)

                Fock = Fock + fac * Conjg(fk(1)) * fk(2) * Fockmat(m,n+d,:,:)
                Rdm = Rdm + fac * Conjg(fk(1)) * fk(2) * Rdmmat(m,n+d,:,:)
                Ovlp = Ovlp + fac * Conjg(fk(1)) * fk(2) * Smat(m,n+d)

                Fock = Fock + fac * Conjg(fk(1)) * fk(2) * Conjg(Fockmat(m+d,n,:,:))
                Rdm = Rdm + fac * Conjg(fk(1)) * fk(2) * Conjg(Rdmmat(m+d,n,:,:))
                Ovlp = Ovlp + fac * Conjg(fk(2)) * fk(1) * Smat(m+d,n)
            End If
        End Do
        End Do
        End Do
        End Do
        Fock = Fock - E0 * Rdm
        Fock = (Fock + Conjg(Transpose(Fock))) / Ovlp
        End Subroutine BuildFock

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

        Subroutine FockKernel(H1,gama,rho,dm0,det,NSO,G)
!
!   Un-Hermitized Fock integrand X(g) (see module header), times det(g).
!   rho = rho(g), dm0 = reference density rho = C C^dag, H1+gama = Gamma(g).
!
        Implicit None
        Integer,            Intent(in)  :: NSO
        Complex(kind=pr),   Intent(in)  :: H1(NSO,NSO)
        Complex(kind=pr),   Intent(in)  :: gama(NSO,NSO), rho(NSO,NSO)
        Complex(kind=pr),   Intent(in)  :: dm0(NSO,NSO), det
        Complex(kind=pr),   Intent(out) :: G(NSO,NSO)
        Complex(kind=pr) :: h, Id(NSO,NSO), Gam(NSO,NSO), tmp(NSO,NSO)
        Integer :: I
        Id = Zero
        Do I = 1, NSO
            Id(I,I) = One
        End Do
        Gam = H1 + gama
! scalar kernel H(g) = Tr[(h + Gamma/2) rho(g)]
        h = sum(transpose(H1) * rho) + sum(gama * transpose(rho)) / 2
! Eq. 30 integrand (before -E):  H(g) rho(g) + (1-rho(g)) Gamma(g) rho(g)
        tmp = MatMul(Gam, rho)
        G = h * rho + MatMul(Id - rho, tmp)
! Eq. 32b, F_-:  1/2 rho Gamma(g) rho(g)
        G = G + F12 * MatMul(dm0, tmp)
! Eq. 32a, F_+:  1/2 (1-rho(g)) Gamma(g) (1-rho)
        tmp = MatMul(Gam, Id - dm0)
        G = G + F12 * MatMul(Id - rho, tmp)
        G = G * det
        End Subroutine FockKernel

        Subroutine MakeGamma(H2,rho,NSO,gama)
        Implicit None
        Integer,            Intent(in)  :: NSO
        Complex(kind=pr),   Intent(in)  :: H2(NSO,NSO,NSO,NSO)
        Complex(kind=pr),   Intent(in)  :: rho(NSO,NSO)
        Complex(kind=pr),   Intent(out) :: gama(NSO,NSO)
        Integer :: i,k
        !$omp parallel do private(i,k)
        Do i = 1,NSO
        Do k = 1,NSO
            gama(i,k) =  sum(H2(i,:,k,:) * transpose(rho))
        End Do
        End Do
        !$omp end parallel do
        End Subroutine MakeGamma

        Subroutine Makeh(H1,H2,gama,rho,NSO,h)
        Implicit None
        Integer,            Intent(in)  :: NSO
        Complex(kind=pr),   Intent(in)  :: H1(NSO,NSO), H2(NSO,NSO,NSO,NSO)
        Complex(kind=pr),   Intent(in)  :: gama(NSO,NSO), rho(NSO,NSO)
        Complex(kind=pr),   Intent(out) :: h
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

    End Module FockTools
