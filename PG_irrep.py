import numpy as np
from pyscf import symm
# get Ops matrix in Symmetry-adaptive basis, only for abelian group
#  C2v: A1 A2 B1 B2 -> 0 1 2 3 
#  D2h: Ag B1g B2g B3g Au B1u B2u B3u -> 0 1 2 3 4 5 6 7
# C4v / D4h characters for linear molecules.  sigma_v = xz, yz; C2' along x, y.  Correlation with Coov / Dooh:
# Sigma+ -> A1 (C4v), Sigma_g+ -> A1g, Sigma_u+ -> A2u, Sigma_u- -> A1u, Pi -> E, Delta -> B1 + B2 (project B1 or B2)
_D4H_CLASSES = ['E', 'C4', 'C2', "C2'", "C2''", 'i', 'S4', 'sh', 'sv', 'sd']
_D4H_CHAR = {
	'A1g': [1, 1, 1, 1, 1, 1, 1, 1, 1, 1],    'A2g': [1, 1, 1, -1, -1, 1, 1, 1, -1, -1],
	'B1g': [1, -1, 1, 1, -1, 1, -1, 1, 1, -1], 'B2g': [1, -1, 1, -1, 1, 1, -1, 1, -1, 1],
	'Eg':  [2, 0, -2, 0, 0, 2, 0, -2, 0, 0],
	'A1u': [1, 1, 1, 1, 1, -1, -1, -1, -1, -1], 'A2u': [1, 1, 1, -1, -1, -1, -1, -1, 1, 1],
	'B1u': [1, -1, 1, 1, -1, -1, 1, -1, -1, 1], 'B2u': [1, -1, 1, -1, 1, -1, 1, -1, 1, -1],
	'Eu':  [2, 0, -2, 0, 0, -2, 0, 2, 0, 0]}
_C4V_CLASSES = ['E', 'C4', 'C2', 'sv', 'sd']
_C4V_CHAR = {'A1': [1, 1, 1, 1, 1], 'A2': [1, 1, 1, -1, -1], 'B1': [1, -1, 1, 1, -1],
	'B2': [1, -1, 1, -1, 1], 'E': [2, 0, -2, 0, 0]}
# class of i^inv sigma_xz^refl C4^k, as _CLASS[inv][refl][k]
_CLASS = [[['E', 'C4', 'C2', 'C4'], ['sv', 'sd', 'sv', 'sd']],
	[['i', 'S4', 'sh', 'S4'], ["C2'", "C2''", "C2'", "C2''"]]]

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
	if (PG == 'C4v' or PG == 'D4h'):
		# linear molecule (mol.symmetry = True -> Coov / Dooh; D4h needs Dooh).  PySCF linear ids: id % 10 is the C2v / D2h id,
		# id // 10 and the C2v / D2h id give Lambda, and the x / y partners differ in the last bit (id ^ 1).
		orbid = np.array(orbsym)
		d2 = orbid % 10
		if (mol.groupname == 'Coov'):
			sxz = (-1.0)**d2 # sigma xz (C2v id)
			par = np.ones(NAO)
			lam = 2 * (orbid // 10) + (d2 >= 2)
		if (mol.groupname == 'Dooh'): # also for C4v on a Dooh molecule (i not used)
			sxz = (-1.0)**(d2 + d2 // 4) # sigma xz = i C2y (D2h id)
			par = (-1.0)**(d2 // 4) # i
			lam = 2 * (orbid // 10) + np.isin(d2, [2, 3, 6, 7])
		C4 = np.zeros([NAO,NAO]) # C4z: Sigma -> 1, (x, y) -> [[cos L pi/2, -sin L pi/2], [sin L pi/2, cos L pi/2]]
		for p in range(NAO):
			if (lam[p] == 0):
				C4[p,p] = 1
		for d in sorted(set(orbid[(lam > 0) & (sxz > 0)])):
			xs = np.where(orbid == d)[0]
			ys = np.where(orbid == (d ^ 1))[0]
			c, s = np.cos(lam[xs[0]]*np.pi/2), np.sin(lam[xs[0]]*np.pi/2)
			# C4 maps the k-th x symmetry-adapted AO to c X_k + s Y_k (PySCF builds the x and y blocks of mol.symm_orb
			# in parallel).  Expand the x and y MOs in their blocks (cx, cy); then C4 x_i = c x_i + s sum_j T[j,i] y_j,
			# T = cy^-1 cx.  T = 1 when the x and y MOs are equivalent; otherwise (e.g. ROHF with only pi_y singly
			# occupied) the k-th x and k-th y MOs are not partners and T mixes them.
			Bx = mol.symm_orb[list(mol.irrep_id).index(d)]
			By = mol.symm_orb[list(mol.irrep_id).index(d ^ 1)]
			cx = np.linalg.lstsq(Bx, RHF.mo_coeff[:,xs], rcond=None)[0]
			cy = np.linalg.lstsq(By, RHF.mo_coeff[:,ys], rcond=None)[0]
			T = np.linalg.solve(cy, cx)
			C4[np.ix_(xs,xs)] = c * np.eye(len(xs))
			C4[np.ix_(ys,ys)] = c * np.eye(len(ys))
			C4[np.ix_(ys,xs)] = s * T
			C4[np.ix_(xs,ys)] = -s * np.linalg.inv(T)
		# elements i^inv sigma_xz^refl C4^k; characters by class
		if (PG == 'C4v'):
			classes, table, invs = _C4V_CLASSES, _C4V_CHAR, [0]
		if (PG == 'D4h'):
			classes, table, invs = _D4H_CLASSES, _D4H_CHAR, [0, 1]
		chi = dict(zip(classes, table[irrep]))
		nop = 8 * len(invs)
		ID = np.ones([nop,nci,nci])
		Ops = np.zeros([nop,NAO,NAO])
		n = 0
		for inv in invs:
			for refl in [0, 1]:
				for k in range(4):
					Ops[n,:,:] = np.linalg.matrix_power(C4, k)
					if (refl == 1):
						Ops[n,:,:] = np.diag(sxz) @ Ops[n,:,:]
					if (inv == 1):
						Ops[n,:,:] = np.diag(par) @ Ops[n,:,:]
					ID[n,:,:] = chi[_CLASS[inv][refl][k]]
					n += 1
	return nci, ID, Ops 


