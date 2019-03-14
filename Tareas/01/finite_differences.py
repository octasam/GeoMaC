
import numpy as np

def initRHS(Nx, Ny):
	return np.zeros((Ny,Nx))

def initSol(Nx, Ny):
	return np.zeros((Ny+2, Nx+2))

def u_face(u1, u2):
    return 0.5 * (u1 + u2)

def laplacian2D_stat(Nx, Ny, r, k, ns = 0):
    """ Esta funcion calcula los coeficientes del 
    sistema lineal producido por el operador de 
    Laplace en 2D. Estos coeficientes son almacenados 
    en la matriz pentadiagonal correspondiente."""
    N = Nx * Ny
    A = np.zeros((N,N))

# Primero llena los bloques tridiagonales
    for j in range(1,Ny+1):
        ofs = Nx * (j-1)     
        # Primer renglón del bloque, considera BC en la pared izq.
        k1 = u_face(k[0,j  ], k[1,j]) # k_(i-1/2, j)
        k2 = u_face(k[2,j  ], k[1,j]) # k_(i+1/2, j)
        k3 = u_face(k[1,j+1], k[1,j]) # k_(i, j-1/2)
        k4 = u_face(k[1,j-1], k[1,j]) # k_(i, j+1/2)
        A[ofs    , ofs] = 1*ns + r * (k1 + k2 + k3 + k4) 
        A[ofs + 1, ofs] = -r * k2

        # Renglones intermedios del bloque 
        for i in range(2,Nx):
            k1 = u_face(k[i-1,j], k[i,j]) # k_(i-1/2, j)
            k2 = u_face(k[i+1,j], k[i,j]) # k_(i+1/2, j)
            k3 = u_face(k[i,j-1], k[i,j]) # k_(i, j-1/2)
            k4 = u_face(k[i,j+1], k[i,j]) # k_(i, j+1/2)
            I = ofs + i - 1
            A[I  , I] = 1*ns + r * (k1 + k2 + k3 + k4)
            A[I-1, I] = -r * k1
            A[I+1, I] = -r * k2

        # Último renglón del bloque, considera BC en la pared der.
        k1 = u_face(k[Nx-1,j  ], k[Nx,j]) # k_(i-1/2, j)
        k2 = u_face(k[Nx+1,j  ], k[Nx,j]) # k_(i+1/2, j)
        k3 = u_face(k[Nx  ,j-1], k[Nx,j]) # k_(i, j-1/2)
        k4 = u_face(k[Nx  ,j+1], k[Nx,j]) # k_(i, j+1/2)
        I = ofs + Nx - 1
        A[I-1,I] = -r * k1 
        A[I  ,I] = 1*ns + r * (k1 + k2 + k3 + k4) 

       
# Despues llena las dos diagonales externas
    I = 0
    for j in range(1, Ny):
        for i in range(1,Nx+1):
            k3 = u_face(k[i,j-1+1], k[i,j+1]) # k_(i, j-1/2)
            k4 = u_face(k[i,j+1], k[i,j]) # k_(i, j+1/2)
            A[I   , I+Nx] = -r * k3 # South, 3, down
            A[I+Nx, I   ] = -r * k4 # North, 4, up
            I += 1
    
    return A


def dirichlet(f, u, boundA, boundB, boundC, boundD):
	Nx = f.shape[1]
	Ny = f.shape[0]

	f[Ny-1,:   ] += boundA # Top wall
	f[0   ,:   ] += boundB # Bot wall
	f[:   ,0   ] += boundC # Left wall
	f[:   ,Nx-1] += boundD # Right wall

	u[Ny+1,:   ] = boundB # Top wall
	u[0   ,:   ] = boundA # Bot wall
	u[:   ,0   ] = boundC # Left wall
	u[:   ,Nx+1] = boundD # Right wall

