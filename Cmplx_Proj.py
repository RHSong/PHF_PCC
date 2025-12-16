import numpy as np
def CmplxProj(MO,NSO,ncik):
	Rk = np.zeros([ncik,NSO,NSO], dtype=complex)
	Rk[0,:,:] = np.eye(NSO)
	if (ncik == 2):
		Rk[1,:,:] = CmplxTrans(MO)
	return Rk

def CmplxTrans(MO):
	return MO.T.conj() @ MO.conj()

