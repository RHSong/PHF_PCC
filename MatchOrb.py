import numpy as np

def MatchOrb(MOs,Ref,nocc):
	occ0 = Ref[:,:nocc]
	occ1 = MOs[:,:nocc]
	vir0 = Ref[:,nocc:]
	vir1 = MOs[:,nocc:]
	m = occ0.T.conj() @ occ1
	u,s,vh = np.linalg.svd(m)
	v = vh.conj().T
	uh = u.conj().T
	print("occ ovlp=",s)
	occ = occ1 @ v @ uh
	m = vir0.T.conj() @ vir1
	u,s,vh = np.linalg.svd(m)
	v = vh.conj().T
	uh = u.conj().T
	print("vir ovlp=",s)
	vir = vir1 @ v @ uh
	newMO = 0 * Ref
	newMO[:,:nocc] = occ
	newMO[:,nocc:] = vir
	return newMO
