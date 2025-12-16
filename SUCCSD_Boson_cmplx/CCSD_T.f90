
   Module CCSD_T

   Use Precision
   Use Constants
   Implicit None

   Contains

      Subroutine CCSD_T_Ene(F, H2, T1, T2, NOcc, NVir, NSO, Ene)
      Implicit None
      Integer, Intent(In) :: NSO, NOcc, NVir
	  Complex (Kind=pr), Intent(In) :: T1(nvir,nocc), T2(nvir,nvir,nocc,nocc)
	  Complex (Kind=pr), Intent(In) :: F(NSO,NSO), H2(NSO,NSO,NSO,NSO)
	  Complex (Kind=pr), Intent(Out) :: Ene
	  Complex (Kind=pr) :: fvo(nvir,nocc), fd(nso), vvvo(nvir,nvir,nvir,nocc)
	  Complex (Kind=pr) :: ovoo(nocc,nvir,nocc,nocc), vvoo(nvir,nvir,nocc,nocc)
	  Complex (Kind=pr) :: eijk(nocc,nocc,nocc), eabc(nvir,nvir,nvir)
	  Complex (Kind=pr) :: wabc(nocc,nocc,nocc), vabc(nocc,nocc,nocc)
	  Complex (Kind=pr) :: wcab(nocc,nocc,nocc), vcab(nocc,nocc,nocc)
	  Complex (Kind=pr) :: wbac(nocc,nocc,nocc), vbac(nocc,nocc,nocc)
	  Complex (Kind=pr) :: w(nocc,nocc,nocc), v(nocc,nocc,nocc)
	  Integer :: a, b, c
	  fvo = F(nocc+1:nso, 1:nocc)
	  vvvo = H2(nocc+1:nso, nocc+1:nso, nocc+1:nso, 1:nocc)
	  ovoo = H2(1:nocc, nocc+1:nso, 1:nocc, 1:nocc)
	  vvoo = H2(nocc+1:nso, nocc+1:nso, 1:nocc, 1:nocc)
	  Do a = 1, nso
	  	  fd(a) = F(a,a)
	  End Do
	  Call DirectSum3(eijk,fd(1:nocc),fd(1:nocc),fd(1:nocc),nocc)
	  Call DirectSum3(eabc,fd(nocc+1:nso),fd(nocc+1:nso),fd(nocc+1:nso),nvir)
	  Ene = Zero
	  !$omp parallel do reduction(+:Ene), &
	  !$omp& private(w, wabc, wcab, wbac, v, vabc, vcab, vbac)
	  Do a = 1, nvir
	  Do b = 1, a
	  Do c = 1, b
	  	  Call get_wv(wabc,vabc,T1,T2,a,b,c,fvo,vvvo,ovoo,vvoo,nocc,nvir)
		  Call get_wv(wcab,vcab,T1,T2,c,a,b,fvo,vvvo,ovoo,vvoo,nocc,nvir)
		  Call get_wv(wbac,vbac,T1,T2,b,a,c,fvo,vvvo,ovoo,vvoo,nocc,nvir)
		  w = wabc + wcab - wbac
		  v = vabc + vcab - vbac
		  w = w / (eijk - eabc(a,b,c))
		  Ene = Ene + sum(w * conjg(v))
	  End Do
	  End Do
	  End Do
	  !$omp end parallel do
	  Ene = Ene / Two
      End Subroutine CCSD_T_Ene


      Subroutine get_wv(w,v,T1,T2,a,b,c,fvo,vvvo,ovoo,vvoo,nocc,nvir)
      Implicit None
	  Integer, Intent(In) :: nocc, nvir, a, b, c
      Complex (Kind=pr), Intent(In) :: T1(nvir,nocc), T2(nvir,nvir,nocc,nocc)
	  Complex (Kind=pr), Intent(In) :: fvo(nvir,nocc), vvvo(nvir,nvir,nvir,nocc)
	  Complex (Kind=pr), Intent(In) :: ovoo(nocc,nvir,nocc,nocc), vvoo(nvir,nvir,nocc,nocc)
	  Complex (Kind=pr), Intent(Out) :: w(nocc,nocc,nocc), v(nocc,nocc,nocc)
      Integer :: I, J, K, L
	  w = Zero
	  Call At_dot_B(vvvo(b,c,:,:), T2(a,:,:,:), w, nvir, nocc, nocc*nocc)
	  Call A_dot_B(-T2(b,c,:,:), ovoo(:,a,:,:), w, nocc, nocc, nocc*nocc)
	  v = Zero
	  Call A_outer_B(T1(a,:), vvoo(b,c,:,:), v, nocc, nocc*nocc)
	  Call A_outer_B(fvo(a,:), T2(b,c,:,:), v, nocc, nocc*nocc)
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
	  !$omp parallel do
	  Do i = 1, n
	  Do j = 1, n
	  Do k = 1, n
	  	  w(i,j,k) = tmp(i,j,k) + tmp(k,i,j) + tmp(j,k,i)
	  End Do
	  End Do
	  End Do
	  !$omp end parallel do
	  End Subroutine Perm3

	  Subroutine DirectSum3(w,a,b,c,n)
	  Implicit None
	  Integer, Intent(In) :: n
	  Complex (Kind=pr), Intent(In) :: a(n), b(n), c(n)
	  Complex (Kind=pr), Intent(Out) :: w(n,n,n)
	  Integer :: i,j,k
	  w = zero
	  !$omp parallel do
	  Do i = 1, n 
	  Do j = 1, n
	  Do k = 1, n
	  	  w(i,j,k) = a(i) + b(j) + c(k)
	  End Do
	  End Do
	  End Do
	  !$omp end parallel do
	  End Subroutine DirectSum3


   End Module CCSD_T

