import numpy as np
import cmath


def solve_ci(Hmat, Nmat):
    """
      solve_ci :

        Solve the generalized eigenvalue problem associated with the
        diagonalization of the Hamiltonian matrix among a set of
        non-orthogonal states. The Hamiltonian (Hmat) and overlap
        matrices (Nmat) should be provided on input.


      input arguments

        Hmat  - Hamiltonian matrix among all states
        Nmat  - overlap matrix among all states

      output arguments

        E     - Hamiltonian eigenvalue corresponding to the g.s.
        civec - Hamiltonian eigenvector corresponding to the g.s.



      some parameters

        thresh - threshold to eliminate linearly independent states
    """

    thresh = 1.0e-12

    # determine number of states stored in array

    idim = len(Hmat)

    # diagonalize overlap matrix

    Nmat = 1/2 * (Nmat + Nmat.conj().T)

    valN, vecN = np.linalg.eigh(-Nmat)

    maxk = idim

    for k in range(1, idim+1):
        if (abs(valN[k-1])/abs(valN[0]) < thresh):
            maxk = k-1
            break

    X = np.zeros((idim, maxk), dtype=complex)

    for k in range(1, maxk+1):
        fac = 1/cmath.sqrt(-valN[k-1])
        X[:, k-1] = fac * vecN[:, k-1]

    # project Hamiltonian matrix into linearly independent subspace

    Hlin = X.conj().T @ Hmat @ X
    Hlin = 1/2 * (Hlin + Hlin.conj().T)

    # diagonalize linearly-independent Hamiltonian matrix

    valH, vecH = np.linalg.eigh(Hlin)

    E = valH[0]

    # transform eigenvectors to original basis

    vecH1 = X @ vecH
    civec = vecH1[:, 0]

    return E, civec

def Lanczos(A,m,nb=1):
	dim = len(A)
	n = nb
	nmax = dim // 2
	Kspace = np.zeros([dim,nmax])
	v0 = np.eye(dim,n)
	Kspace[:,:n] = v0[:,:]
	tmp = v0
	for i in range(n,nmax-n,n):
		tmp = A @ tmp
		Kspace[:,i:i+n] = tmp[:,:]
	Q = np.linalg.qr(Kspace)[0]
	mat = Q.T @ A @ Q
	evals, evecs = np.linalg.eigh(mat)
	return evals[:m], Q @ evecs[:,:m]

def Davidson(A,n):
	dim = len(A)
	k = 8
	kmax = dim // 2
	t = np.eye(dim,k)
	V = np.zeros([dim,kmax])
	I = np.eye(dim)
	for m in range(k,kmax-k,k):
		if m <= k:
			V[:,:m] = t[:,:m]
			eold = 1
		elif m > k:
			eold = evals[:n]
		V[:,:m], R = np.linalg.qr(V[:,:m])
		T = V[:,:m].T @ A @ V[:,:m]
		evals, evecs = np.linalg.eigh(T)
		for j in range(0,k):
			w = (A-evals[j]*I) @ V[:,:m] @ evecs[:,j]
			q = w / (evals[j] - A[j,j])
			V[:,m+j] = q
		norm = np.linalg.norm(evals[:n] - eold)
		if norm < 1e-6:
			break
	return evals[:n]

