import numpy as np
import copy
def DIIS(Errlist):
# build B matrix
	length = len(Errlist)
	dim = length+1
	B=np.ones([dim,dim])
	for i in range(length):
		for j in range(length):
			B[i,j] = np.dot(Errlist[i],Errlist[j])
	B[length,length]=0
# build V
	V=np.zeros(dim)
	V[length]=1
# Solve B*C=V
#	print(B)
	C = np.linalg.solve(B, V)
	return C

def FockDIIS(Flist,Dlist):
	Errlist=[]
	for i in range(len(Flist)):
		err = np.matmul(Flist[i],Dlist[i])-np.matmul(Dlist[i],Flist[i])
		Errlist.append(copy.deepcopy(np.ravel(err)))
	coeff = DIIS(Errlist)
	F = 0
	for i in range(len(coeff)-1):
		F = F + coeff[i]*Flist[i]
	return F

def UHFDIIS(FAlist,FBlist,DMAlist,DMBlist):
	Errlist=[]
	for i in range(len(FAlist)):
		err1 = np.matmul(FAlist[i],DMAlist[i])-np.matmul(DMAlist[i],FAlist[i])
		err2 = np.matmul(FBlist[i],DMBlist[i])-np.matmul(DMBlist[i],FBlist[i])
		Errlist.append(copy.deepcopy(np.ravel([err1,err2])))
	coeff = DIIS(Errlist)
	FA = 0
	FB = 0
	for i in range(len(coeff)-1):
		FA = FA + coeff[i]*FAlist[i]
		FB = FB + coeff[i]*FBlist[i]
	return FA,FB

