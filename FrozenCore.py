import numpy as np
from parameter import *
#use rRHF MO energy as reference
# E0 = 2 * <i|h|i> + 2 * [ii|jj] - [ij|ji], h2 is in Mulliken order
# <p|h|q>  => <p|h|q> + 2 [ii|pq] - [pi|iq]
def FrozenCore(H1,H2,MO,NFC,NFV):
	HOne = ao2mo(H1,MO,2)
	HTwo = ao2mo(H2,MO,4)
	NMO = NAO - NFC - NFV
	E0 = 0
	H1new = np.zeros([NMO,NMO])
	H2new = np.zeros([NMO,NMO,NMO,NMO])
	for i in range(NFC):
		E0 += 2 * HOne[i,i]
		for j in range(NFC):
			E0 += 2 * HTwo[i,i,j,j] - HTwo[i,j,j,i]
	for p in range(NFC,NAO-NFV):
		for q in range(NFC,NAO-NFV):
			H1new[p-NFC,q-NFC] += HOne[p,q]
			for i in range(NFC):
				H1new[p-NFC,q-NFC] += 2 * HTwo[i,i,p,q] - HTwo[p,i,i,q]
	for p in range(NFC,NAO-NFV):
		for q in range(NFC,NAO-NFV):
			for r in range(NFC,NAO-NFV):
				for s in range(NFC,NAO-NFV):
					H2new[p-NFC,q-NFC,r-NFC,s-NFC] = HTwo[p,q,r,s]
	return E0, H1new, H2new
