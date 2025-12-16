
   Module Spin
   Use Precision
   Use Constants
   Use IntTrans
   Use FPUCC
   Implicit None
   Contains

      Subroutine BuildSz(SOne,NAO,NSO)
      Implicit None
      Integer,        Intent(In)    :: NAO, NSO
      Complex (Kind=pr), Intent(Out)   :: SOne(NSO,NSO)
      Integer :: i, j
      Do i = 1, NAO
        SOne(i,i) = 0.5
        SOne(i+NAO,i+NAO) = -0.5
      End Do
      End Subroutine BuildSz

      Subroutine BuildS2(SOne,STwo,NAO,NSO)
      Implicit None
      Integer,        Intent(In)    :: NAO, NSO
      Complex (Kind=pr), Intent(Out)   :: SOne(NSO,NSO),STwo(NSO,NSO,NSO,NSO)
      Integer :: i, j
      Do i = 1, NAO
         SOne(i,i) = 0.75
         SOne(i+NAO,i+NAO) = 0.75
      End Do

      Do i = 1, NAO
      Do j = 1, NAO
         STwo(i+NAO,i+NAO,j+NAO,j+NAO) = 0.25
         STwo(i,i,j,j) = 0.25
         STwo(i,i+NAO,j+NAO,j) = 1
         STwo(i,i,j+NAO,j+NAO) = -0.5
      End Do
      End Do
      Return
      End Subroutine BuildS2

      Subroutine CalcS2(Ref,S2,Sz,T1,T2,NOcc,NAO,NSO,QJ,SP,&
	  					NPoints, Npg, ncisp, ncipg, ncik, R1, R2, Rpg, Rk, &
						Roota,Rootb,Rooty,Weightsp,Weightpg,fsp,fpg,fk,TrunODE,comm)
      Implicit None
      Integer,           Intent(In)    :: NOcc, NAO, NSO
      Integer,           Intent(In)    :: NPoints(2), TrunODE, QJ,SP, Npg,comm
	  Integer,           Intent(In)    :: ncisp, ncipg, ncik
      Real    (Kind=pr), Intent(In)    :: Roota(:), Rootb(:), Rooty(:)
      Real    (Kind=pr), Intent(In)    :: Weightsp(:,:), Weightpg(Npg,ncipg,ncipg)
	  Complex (Kind=pr), Intent(In)    :: fsp(ncisp),fpg(ncipg),fk(ncik)
	  Complex (Kind=pr), Intent(In)    :: R1(:,:,:), R2(:,:,:)
	  Complex (Kind=pr), Intent(In)    :: Rpg(npg,NSO,NSO), Rk(ncik,NSO,NSO)
      Complex (Kind=pr), Intent(In)    :: Ref(NSO,NSO)
      Complex (Kind=pr), Intent(In)    :: T1(NOcc+1:NSO,NOcc)
      Complex (Kind=pr), Intent(In)    :: T2(NOcc+1:NSO,NOcc+1:NSO,NOcc,NOcc)
      Complex (Kind=pr), Intent(Out)   :: S2, Sz
      Complex (Kind=pr) :: SOne(NSO,NSO),STwo(NSO,NSO,NSO,NSO)
      Complex (Kind=pr) :: Olap0, EOlap0
      Complex (Kind=pr) :: OlapEx1(NOcc+1:NSO,NOcc)
      Complex (Kind=pr) :: OlapEx2(NOcc+1:NSO,NOcc+1:NSO,NOcc,NOcc)
      Complex (Kind=pr) :: EOlapEx1(NOcc+1:NSO,NOcc)
      Complex (Kind=pr) :: EOlapEx2(NOcc+1:NSO,NOcc+1:NSO,NOcc,NOcc)
      Real (Kind=pr) :: ENuc = 0
      
      Call BuildS2(SOne,STwo,NAO,NSO)
      STwo = 2*STwo
      Call Mulliken2Dirac(STwo,NSO)
      Call ASymmH2(STwo,NSO)
      Call IntTran2(SOne,Ref,NSO)
      Call IntTran4(STwo,Ref,NSO)
      Call PbarHbarOlap(Olap0,OlapEx1,OlapEx2,EOlap0,EOlapEx1,        &
                        EOlapEx2,SOne,STwo,Ref,T1,T2,NOcc,NAO,NSO,ENuc,QJ,SP,&
                        NPoints, Npg, ncisp, ncipg, ncik, R1, R2, Rpg, Rk, &
                        Roota,Rootb,Rooty,Weightsp,Weightpg,fsp,fpg,fk,TrunODE,comm)
      S2 = EOlap0 / Olap0
      Call BuildSz(SOne,NAO,NSO)
      Call IntTran2(SOne,Ref,NSO)
      Call PbarHbarOlap(Olap0,OlapEx1,OlapEx2,EOlap0,EOlapEx1,        &
                        EOlapEx2,SOne,0*STwo,Ref,T1,T2,NOcc,NAO,NSO,ENuc,QJ,SP,&
                        NPoints, Npg, ncisp, ncipg, ncik, R1, R2, Rpg, Rk, &
                        Roota,Rootb,Rooty,Weightsp,Weightpg,fsp,fpg,fk,TrunODE,comm)
      Sz = EOlap0 / Olap0
      End Subroutine CalcS2

   End Module Spin

