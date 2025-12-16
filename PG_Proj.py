import numpy as np
from scipy.linalg import block_diag
from ao2mo import *
from PGSymm import *
def PG_Proj(mol, PG, irrep, ovlp, OrthAO, NAO, nx=0, ny=0, kx=1, ky=1):
	if (PG == 'Hub'):
		ID, nci, R = HubTable(nx, ny, irrep, kx, ky)
	if (PG == 'C2v'):
		ID, nci, R = C2vTable(mol,irrep)
	if (PG == 'C4v'):
		ID, nci, R = C4vTable(mol,irrep)
	if (PG == 'D2h'):
		ID, nci, R = D2hTable(mol,irrep)
	if (PG == 'D4h'):
		ID, nci, R = D4hTable(mol,irrep)
	if (PG == 'None'):
		nci = 1
		ID = np.ones([1,nci,nci])
		R = np.zeros([1,NAO,NAO])
		R[0,:,:] = np.eye(NAO)
	dim = np.shape(R)
	newR = np.zeros([dim[0], 2*dim[1], 2*dim[2]], dtype=complex)
	for i in range(len(ID)):
		tmp = ao2mo(R[i,:,:] @ ovlp, OrthAO, 2)
		newR[i,:,:] = block_diag(tmp,tmp)
	return nci, ID, newR
