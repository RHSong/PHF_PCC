import numpy as np
from scipy import linalg
def HubTable(nx, ny, irrep, kx, ky):
	NAO = nx * ny
	if (irrep == 'PBCX'):
		ID, nci, R = HubSymmx(nx, ny, kx, ky, Transx=True, sigmax=False)
	if (irrep == 'PBCY'):
		ID, nci, R = HubSymmy(nx, ny, kx, ky, Transy=True, sigmay=False)
	if (irrep == 'SIGX'):
		ID, nci, R = HubSymmx(nx, ny, kx, ky, Transx=False, sigmax=True)
	if (irrep == 'SIGY'):
		ID, nci, R = HubSymmy(nx, ny, kx, ky, Transy=False, sigmay=True)
	if (irrep == 'PBCXY'):
		IDx, ncix, Rx = HubSymmx(nx, ny, kx, ky, Transx=True, sigmax=False)
		IDy, nciy, Ry = HubSymmy(nx, ny, kx, ky, Transy=True, sigmay=False)
		nop = len(IDx) * len(IDy)
		nci = 1
		ID = np.ones([nop,nci,nci])
		R = np.zeros([nop,NAO,NAO], dtype=complex)
		n = 0
		for i in range(len(IDx)):
			for j in range(len(IDy)):
				n += 1
				R[n,:,:] = Rx[i,:,:] @ Ry[j,:,:]
	return ID, nci, R

def HubSymmx(nx, ny, kx, ky, Transx=True, sigmax=False):
	nopx = 1
	nci = 1
	NAO = nx * ny
	F = np.exp(2j * np.pi / kx)
# symm in x-direction
	if (Transx and nx > 1):
		nopx = nx
		Tx = np.zeros([NAO,NAO], dtype=complex)
		for i in range(nx-1):
			for j in range(ny):
				ni = ny * i + j
				nj = ny * (i+1) + j
				Tx[ni,nj] = F
				Tx[ny * (nx-1) + j,j] = F
	if (sigmax and nx > 1):
		nopx = 2
		Tx = np.zeros(NAO,NAO)
		for i in range(nx // 2):
			for j in range(ny):
				ni = ny * i + j
				nj = ny * (nx-i-1) + j
				Tx[ni,nj] = 1
	Rx = np.zeros([nopx,NAO,NAO], dtype=complex)
	Rx[0,:,:] = np.eye(NAO)
	for i in range(1,nopx):
		Rx[i,:NAO,:NAO] = Tx @ Rx[i-1,:NAO,:NAO]
	ID = np.ones([nopx,nci,nci])
	return ID, nci, Rx

def HubSymmy(nx, ny, kx, ky, Transy=True, sigmay=False):
	nopy = 1
	nci = 1
	NAO = nx * ny
	F = np.exp(2j * np.pi / ky)
	if (Transy and ny > 1):
		nopy = ny
		Ty = np.zeros([NAO,NAO], dtype=complex)
		for i in range(nx):
			for j in range(ny-1):
				ni = ny * i + j
				nj = ny * i + j + 1
				Ty[ni,nj] = F
				Ty[ny * i + ny - 1, ny * i] = F
	if (sigmay and ny > 1):
		nopy = 2
		Ty = np.zeros(NAO,NAO)
		for i in range(nx):
			for j in range(ny // 2):
				ni = ny * i + j
				nj = ny * (i+1) - 1 - j
				Ty[ni,nj] = 1
	Ry = np.zeros([nopy,NAO,NAO], dtype=complex)
	Ry[0,:,:] = np.eye(NAO)
	for i in range(1,nopy):
		Ry[i,:NAO,:NAO] = Ty @ Ry[i-1,:NAO,:NAO]
	ID = np.ones([nopy,nci,nci])
	return ID, nci, Ry


def C2vTable(mol,irrep):
	nao = mol.nao
	nop = 4
	nci = 1
	ID = np.ones([nop,nci,nci])
	Ops = np.zeros([nop,nao,nao])
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
	ops = SymmMat(0,2,mol)
	for i in range(nop):
		Ops[i,:mol.nao,:mol.nao] = ops[i,:,:]
	return ID, nci, Ops

def C4vTable(mol,irrep):
	nao = mol.nao
	nop = 8
	nci = 1
	ID = np.ones([nop,nci,nci])
	Ops = np.zeros([nop,nao,nao])
	if (irrep == 'A1'):
		ID = ID
	if (irrep == 'A2'):
		ID[4:,:,:] = -1
	if (irrep == 'B1'):
		ID[1,:,:] = -1
		ID[3,:,:] = -1
		ID[6:,:,:] = -1
	if (irrep == 'B2'):
		ID[1,:,:] = -1
		ID[3,:,:] = -1
		ID[4:6,:,:] = -1
	if (irrep == 'E'):
		nci = 2
		ID = np.zeros([nop,nci,nci])
		ID[0,:,:] = np.eye(2)
		ID[1,:,:] = np.array([[0,1],[-1,0]])
		ID[2,:,:] = -np.eye(2)
		ID[3,:,:] = np.array([[0,-1],[1,0]])
		ID[4,:,:] = np.array([[1,0],[0,-1]]) # sigmaxz
		ID[5,:,:] = np.array([[-1,0],[0,1]]) # sigmayz
		ID[6,:,:] = np.array([[0,-1],[-1,0]])
		ID[7,:,:] = np.array([[0,1],[1,0]])
	ops = SymmMat(0,4,mol)
	for i in range(nop):
		Ops[i,:mol.nao,:mol.nao] = ops[i,:,:]
	return ID, nci, Ops

def D2hTable(mol,irrep):
	nao = mol.nao
	nop = 8
	nci = 1
	ID = np.ones([nop,nci,nci])
	Ops = np.zeros([nop,nao,nao])
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
	ops = SymmMat(1,2,mol)
	for i in range(nop):
		Ops[i,:mol.nao,:mol.nao] = ops[i,:,:]
	return ID, nci, Ops

def D4hTable(mol,irrep):
	nao = mol.nao
	nop = 16
	nci = 1
	ID = np.ones([nop,nci,nci])
	Ops = np.zeros([nop,nao,nao])
	if (irrep == 'A1g'):
		ID = ID
	if (irrep == 'A2g'):
		ID[4:8,:,:] = -1
		ID[12:,:,:] = -1
	if (irrep == 'Eg'):
		nci = 2
		ID = np.zeros([nop,nci,nci])
		ID[0,:,:] = np.eye(2)
		ID[1,:,:] = np.array([[0,1],[-1,0]])
		ID[2,:,:] = -np.eye(2)
		ID[3,:,:] = np.array([[0,-1],[1,0]])
		ID[4,:,:] = np.array([[-1,0],[0,1]]) # c2y
		ID[5,:,:] = np.array([[1,0],[0,-1]]) # c2x
		ID[6,:,:] = np.array([[0,1],[1,0]])
		ID[7,:,:] = np.array([[0,-1],[-1,0]])
		ID[8,:,:] = np.eye(2)
		ID[9,:,:] = np.array([[0,-1],[1,0]])
		ID[10,:,:] = np.array([[0,1],[-1,0]])
		ID[11,:,:] = -np.eye(2)
		ID[12,:,:] = np.array([[1,0],[0,-1]]) # sgxz
		ID[13,:,:] = np.array([[-1,0],[0,1]]) # sgyz
		ID[14,:,:] = np.array([[0,1],[1,0]])
		ID[15,:,:] = np.array([[0,-1],[-1,0]])
	ops = SymmMat(1,4,mol)
	for i in range(nop):
		Ops[i,:mol.nao,:mol.nao] = ops[i,:,:]
	return ID, nci, Ops

'''
Only for diatomic molecules
PG = 0, Cnv
for linear mol, z axis is the main rotation axis.
for c2v, the generator is c2(z), sigma(xz)
for c4v, the generator is c4(z), c2(z), sigma(xz)
for d2h, the generator is c2(z), c2(y), i
for d4h, the generator is c4(z), c2(z), c2(y), i
Build mat for sites and orbitals
'''
def SymmMat(PG,n,mol):
	s1, s2, p1, p2, d1, d2, f1, f2 = GetAO(mol)
	NAO = mol.nao
	OffID = np.zeros([2,2])
	OffID[0,1] = OffID[1,0] = 1
	if (PG == 0 and n == 2):
		nop = 4
		Ops = np.zeros([nop,NAO,NAO])
		es = np.eye(2)
		c2zs = np.eye(2)
		sgxzs = np.eye(2)
		e1 = OrbSymm(s1, p1, d1, f1, 'E')
		c2z1 = OrbSymm(s1, p1, d1, f1, 'c2z')
		sgxz1 = OrbSymm(s1, p1, d1, f1, 'sgxz')
		e2 = OrbSymm(s2, p2, d2, f2, 'E')
		c2z2 = OrbSymm(s2, p2, d2, f2, 'c2z')
		sgxz2 = OrbSymm(s2, p2, d2, f2, 'sgxz')
		Ops[0,:,:] = linalg.block_diag(e1,e2)
		Ops[1,:,:] = linalg.block_diag(c2z1,c2z2)
		Ops[2,:,:] = linalg.block_diag(sgxz1,sgxz2)
		Ops[3,:,:] = Ops[1,:,:] @ Ops[2,:,:] # Sgyz
	if (PG == 0 and n == 4):
		nop = 8
		Ops = np.zeros([nop,NAO,NAO])
		es = np.eye(2)
		c4zs = np.eye(2)
		c2zs = np.eye(2)
		sgxzs = np.eye(2)
		e1 = OrbSymm(s1, p1, d1, f1, 'E')
		c2z1 = OrbSymm(s1, p1, d1, f1, 'c2z')
		c4z1 = OrbSymm(s1, p1, d1, f1, 'c4z')
		sgxz1 = OrbSymm(s1, p1, d1, f1, 'sgxz')
		e2 = OrbSymm(s2, p2, d2, f2, 'E')
		c2z2 = OrbSymm(s2, p2, d2, f2, 'c2z')
		c4z2 = OrbSymm(s2, p2, d2, f2, 'c4z')
		sgxz2 = OrbSymm(s2, p2, d2, f2, 'sgxz')
		Ops[0,:,:] = linalg.block_diag(e1,e2)
		Ops[1,:,:] = linalg.block_diag(c4z1,c4z2)
		Ops[2,:,:] = linalg.block_diag(c2z1,c2z2)
		Ops[3,:,:] = Ops[1,:,:] @ Ops[2,:,:] # c43
		Ops[4,:,:] = linalg.block_diag(sgxz1,sgxz2) # sigmaxz v
		Ops[5,:,:] = Ops[2,:,:] @ Ops[4,:,:] # sigmayz v
		Ops[6,:,:] = Ops[1,:,:] @ Ops[4,:,:] # sigma d
		Ops[7,:,:] = Ops[3,:,:] @ Ops[4,:,:] # sigma d
	if (PG == 1 and n == 2):
		nop = 8
		Ops = np.zeros([nop,NAO,NAO])
		es = np.eye(2)
		c2zs = np.eye(2)
		c2ys = OffID
		Is = OffID
		e1 = OrbSymm(s1, p1, d1, f1, 'E')
		c2z1 = OrbSymm(s1, p1, d1, f1, 'c2z')
		c2y1 = OrbSymm(s1, p1, d1, f1, 'c2y')
		i1 = OrbSymm(s1, p1, d1, f1, 'I')
		Ops[0,:,:] = linalg.block_diag(e1,e1)
		Ops[1,:,:] = linalg.block_diag(c2z1,c2z1)
		Ops[2,:,:] = np.kron(c2ys,c2y1)
		Ops[3,:,:] = Ops[1,:,:] @ Ops[2,:,:] # c2x
		Ops[4,:,:] = np.kron(Is,i1)
		Ops[5,:,:] = Ops[4,:,:] @ Ops[1,:,:] # sigma xy
		Ops[6,:,:] = Ops[4,:,:] @ Ops[2,:,:] # sigma xz
		Ops[7,:,:] = Ops[4,:,:] @ Ops[3,:,:] # sigma yz
	if (PG == 1 and n == 4):
		nop = 16
		Ops = np.zeros([nop,NAO,NAO])
		es = np.eye(2)
		c4zs = np.eye(2)
		c2zs = np.eye(2)
		c2ys = OffID
		Is = OffID
		e1 = OrbSymm(s1, p1, d1, f1, 'E')
		c2z1 = OrbSymm(s1, p1, d1, f1, 'c2z')
		c4z1 = OrbSymm(s1, p1, d1, f1, 'c4z')
		c2y1 = OrbSymm(s1, p1, d1, f1, 'c2y')
		i1 = OrbSymm(s1, p1, d1, f1, 'I')
		Ops[0,:,:] = linalg.block_diag(e1,e1)
		Ops[1,:,:] = linalg.block_diag(c4z1,c4z1)
		Ops[2,:,:] = linalg.block_diag(c2z1,c2z1)
		Ops[3,:,:] = Ops[1,:,:] @ Ops[2,:,:] # c43
		Ops[4,:,:] = np.kron(c2ys,c2y1) # c'2
		Ops[5,:,:] = Ops[2,:,:] @ Ops[4,:,:] # c'2
		Ops[8,:,:] = np.kron(Is,i1)
		Ops[11,:,:] = Ops[8,:,:] @ Ops[2,:,:] # sigma h
		Ops[12,:,:] = Ops[8,:,:] @ Ops[4,:,:] # sigma xz v
		Ops[13,:,:] = Ops[8,:,:] @ Ops[5,:,:] # sigma yz v
		Ops[14,:,:] = Ops[1,:,:] @ Ops[12,:,:] # sigma d
		Ops[15,:,:] = Ops[3,:,:] @ Ops[12,:,:] # sigma d
		Ops[6,:,:] = Ops[1,:,:] @ Ops[4,:,:] # c''2
		Ops[7,:,:] = Ops[3,:,:] @ Ops[4,:,:] # c''2
		Ops[9,:,:] = Ops[11,:,:] @ Ops[1,:,:] # s4
		Ops[10,:,:] = Ops[11,:,:] @ Ops[3,:,:] # s4
	return Ops

'''
Get AO for diatomic molecules
the orbitals are sorted as
s
px, py, pz
dxy, dyz, dz2, dxz, dx2-y2
fy(y2-3x2), fxyz, fyz2, fz3, fxz2, fz(x2-y2), fx(x2-3y2)
'''
def GetAO(mol):
	s0 = mol.search_ao_label('%4s %4s .s' %(0,mol.atom_symbol(0)))
	p0 = mol.search_ao_label('%4s %4s .p' %(0,mol.atom_symbol(0)))
	d0 = mol.search_ao_label('%4s %4s .d' %(0,mol.atom_symbol(0)))
	f0 = mol.search_ao_label('%4s %4s .f' %(0,mol.atom_symbol(0)))
	s1 = mol.search_ao_label('%4s %4s .s' %(1,mol.atom_symbol(1)))
	p1 = mol.search_ao_label('%4s %4s .p' %(1,mol.atom_symbol(1)))
	d1 = mol.search_ao_label('%4s %4s .d' %(1,mol.atom_symbol(1)))
	f1 = mol.search_ao_label('%4s %4s .f' %(1,mol.atom_symbol(1)))
	return s0, s1, p0, p1, d0, d1, f0, f1

'''
spdf belongs to the same atom, let the first s orbital has the index 0
'''
def OrbSymm(s, p, d, f, Symm):
	N = len(s) + len(p) + len(d) + len(f)
	ind = s[0]
	s = s - ind
	px, py, pz = p.reshape([3,-1],order='F') - ind
	dxy, dyz, dz2, dxz, dx2y2 = d.reshape([5,-1],order='F') - ind
	fy3, fxyz, fyz2, fz3, fxz2, fzx2, fx3 = f.reshape([7,-1],order='F') - ind
	R = np.eye(N)
	if (Symm=='c2z'):
		for i in np.concatenate((px,py,dyz,dxz,fy3,fyz2,fxz2,fx3)):
			R[i,i] = -1
	if (Symm=='c4z'):
		for i in np.concatenate((dxy,dx2y2,fxyz,fzx2)):
			R[i,i] = -1
		for i in range(len(px)):
			n = px[i]
			m = py[i]
			R[n,n] = 0
			R[m,m] = 0
			R[n,m] = 1
			R[m,n] = -1
		for i in range(len(dyz)):
			n = dxz[i]
			m = dyz[i]
			R[n,n] = 0
			R[m,m] = 0
			R[n,m] = 1
			R[m,n] = -1
		for i in range(len(fy3)):
			n = fx3[i]
			m = fy3[i]
			R[n,n] = 0
			R[m,m] = 0
			R[n,m] = 1
			R[m,n] = -1
			n = fxz2[i]
			m = fyz2[i]
			R[n,n] = 0
			R[m,m] = 0
			R[n,m] = 1
			R[m,n] = -1
	if (Symm=='c2y'):
		for i in np.concatenate((px,pz,dxy,dyz,fz3,fxz2,fzx2,fx3)):
			R[i,i] = -1
	if (Symm=='sgxz'):
		for i in np.concatenate((py,dxy,dyz,fy3,fxyz,fyz2)):
			R[i,i] = -1
	if (Symm=='I'):
		for i in np.concatenate((px, py, pz, fy3, fxyz, fyz2, fz3, fxz2, fzx2, fx3)):
			R[i,i] = -1
	if (Symm=='E'):
		R = R
	return R


