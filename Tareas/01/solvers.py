
import numpy as np
import time

def solve(A,u,f):
	"""
	Se utiliza un algoritmo del paquete linalg para obtener la solucion del sistema de N x N
	"""
	Nx = f.shape[1]
	Ny = f.shape[0]

	ut = np.copy(u[1:Ny+1,1:Nx+1])

	ut.shape = ut.size   # Cambiamos los arreglos a formato unidimensional
	f.shape = f.size     # Cambiamos los arreglos a formato unidimensional

	t1_start = time.perf_counter()
	ut = np.linalg.solve(A,f)
	t1_stop = time.perf_counter()
	print(time.ctime(), '\n CPU time: {:0.6f} '.format(t1_stop-t1_start))

	ut.shape = (Ny, Nx) # Regresamos el arreglo a formato bidimensional
	f.shape = (Ny, Nx)
	u[1:Ny+1,1:Nx+1] = ut

