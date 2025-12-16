import numpy as np
from pyscf import symm
# get Ops matrix in Symmetry-adaptive basis, only for abelian group
#  C2v: A1 A2 B1 B2 -> 0 1 2 3 
#  D2h: Ag B1g B2g B3g Au B1u B2u B3u -> 0 1 2 3 4 5 6 7
def PG_Proj_SO(mol, RHF, PG, irrep):
	nci = 1
	NAO = mol.nao
	orbsym = symm.label_orb_symm(mol, mol.irrep_id, mol.symm_orb, RHF.mo_coeff)
	if (PG == 'C2v'):
		ID = np.ones([4,nci,nci])
		if (irrep == 'A1'):
			ID = ID
		if (irrep == 'A2'):
			ID[2,:,:] = -1
			ID[3,:,:] = -1
		if (irrep == 'B1'):
			ID[1,:,:] = -1
			ID[3,:,:] = -1
		if (irrep == 'B2'):
			ID[1,:,:] = -1
			ID[2,:,:] = -1
		Ops = np.zeros([4,NAO,NAO])
		Ops[0,:,:] = np.eye(NAO) # E
		tmp = orbsym // 2
		Ops[1,:,:] = np.diag((-1)**tmp) # C2
		Ops[2,:,:] = np.diag((-1)**orbsym) # sigma xz
		Ops[3,:,:] = Ops[1,:,:] @ Ops[2,:,:] # sigma yz
	if (PG == 'D2h'):
		ID = np.ones([8,nci,nci])
		if (irrep == 'Ag'):
			ID = ID
		if (irrep == 'B1g'):
			ID[2:4,:,:] = -1
			ID[6:,:,:] = -1
		if (irrep == 'B2g'):
			ID[1,:,:] = -1
			ID[3,:,:] = -1
			ID[5,:,:] = -1
			ID[7,:,:] = -1
		if (irrep == 'B3g'):
			ID[1:3,:,:] = -1
			ID[5:7,:,:] = -1
		if (irrep == 'Au'):
			ID[4:,:,:] = -1
		if (irrep == 'B1u'):
			ID[2:6,:,:] = -1
		if (irrep == 'B2u'):
			ID[1,:,:] = -1
			ID[3:5,:,:] = -1
			ID[6,:,:] = -1
		if (irrep == 'B3u'):
			ID[1:3,:,:] = -1
			ID[4,:,:] = -1
			ID[7,:,:] = -1
		Ops = np.zeros([8,NAO,NAO])
		Ops[0,:,:] = np.eye(NAO) # E
		tmp = orbsym // 2
		Ops[1,:,:] = np.diag((-1)**tmp) # C2z
		Ops[2,:,:] = np.diag((-1)**orbsym) # C2y
		Ops[3,:,:] = Ops[1,:,:] @ Ops[2,:,:]
		tmp = orbsym // 4
		Ops[4,:,:] = np.diag((-1)**tmp) # i
		Ops[5,:,:] = Ops[4,:,:] @ Ops[1,:,:] # sigma xy
		Ops[6,:,:] = Ops[4,:,:] @ Ops[2,:,:] # sigma xz
		Ops[7,:,:] = Ops[4,:,:] @ Ops[3,:,:] # sigma yz
	return nci, ID, Ops 


