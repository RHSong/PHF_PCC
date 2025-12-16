import numpy as np
from WignerMat import *
from scipy import linalg

def RootWeightSy(ngrid):
	rootb = np.zeros(ngrid)
	weightb = np.zeros(ngrid)
	A = np.zeros(ngrid)
	B = np.zeros(ngrid)
	C = np.zeros(ngrid)
	W = 2
	for K in range(ngrid):
		A[K] = 2 - 1/(K+1)
		C[K] = 1 - 1/(K+1)
	for K in range(ngrid-1):
		weightb[K] = np.sqrt(C[K+1]/(A[K]*A[K+1]))
	mat = np.diag(rootb)
	for i in range(ngrid-1):
		mat[i,i+1]=mat[i+1,i]=weightb[i]
	vals,vecs = linalg.eigh(mat,lower=False)
	for i in range(ngrid):
		weightb[i] = W*vecs[0,i]*vecs[0,i]
	rootb = (vals + 1)*np.pi/2
	weightb = weightb*np.pi/2
	for i in range(ngrid):
		weightb[i] = weightb[i] * np.sin(rootb[i])
	if (ngrid == 1):
		rootb = np.zeros(1)
		weightb = np.ones(1)
	return rootb, weightb

def RootWeight(ngrid):
	nl, ny = ngrid
# first fill gamma
	rooty = np.zeros(ny)
	weighty = np.zeros(ny)
	if (ny ==1):
		rooty[0] = 0
		weighty[0] = 1
	else:
		dx = 2*np.pi / (ny-1)
		rooty = [i*dx for i in range(ny)]
		weighty = [dx for i in range(ny)]
		weighty[0] = dx / 2
		weighty[ny-1] = dx / 2
# Read Lebedev
	if (nl == 1):
		data = np.loadtxt('lebedev_000.txt')
	elif (nl == 6):
		data = np.loadtxt('lebedev_003.txt')
	elif (nl == 14):
		data = np.loadtxt('lebedev_005.txt')
	elif (nl == 26):
		data = np.loadtxt('lebedev_007.txt')
	elif (nl == 38):
		data = np.loadtxt('lebedev_009.txt')
	elif (nl == 50):
		data = np.loadtxt('lebedev_011.txt')
	else:
		print('Not support, check https://people.sc.fsu.edu/~jburkardt/datasets/sphere_lebedev_rule/sphere_lebedev_rule.html')
	roota = np.transpose(data)[0] / 180 * np.pi
	rootb = np.transpose(data)[1] / 180 * np.pi
	weightlbd = np.transpose(data)[2]
	weight = np.zeros([nl,ny])
	for l in range(nl):
		for y in range(ny):
			weight[l,y] = weightlbd[l] * weighty[y]
	return roota, rootb, rooty, weight

def BuildRotMat(a,b,y,NAO,NSO):
	Ca = np.cos(a/2)
	Sa = np.sin(a/2)
	Cb = np.cos(b/2)
	Sb = np.sin(b/2)
	Cy = np.cos(y/2)
	Sy = np.sin(y/2)
	SaRotMat = np.zeros([NSO,NSO],dtype=complex)
	SbRotMat = np.zeros([NSO,NSO],dtype=complex)
	SyRotMat = np.zeros([NSO,NSO],dtype=complex)
	for i in range(NAO):
		SaRotMat[i][i] = complex(Ca,Sa)
		SaRotMat[i+NAO][i+NAO] = complex(Ca,-Sa)
		SbRotMat[i][i] = SbRotMat[i+NAO][i+NAO] = Cb
		SbRotMat[i][i+NAO] = Sb
		SbRotMat[i+NAO][i] = -Sb
		SyRotMat[i][i] = complex(Cy,Sy)
		SyRotMat[i+NAO][i+NAO] = complex(Cy,-Sy)
	RotMat = SaRotMat @ SbRotMat @ SyRotMat
	return RotMat

