import read_data as rd
import meshes as ms
import phys_data as pd 
import finite_differences as fd 
import solvers as sol
import visflow as vis
import numpy as np
#
# Lee la información de la conductividad térmica
#
k = pd.readImage('TUX50.png',0.01)
#
# Lee los datos de la simulación numérica
#
(a,b,c,d,Nx,Ny,bA,bB,bC,bD,ofilename) = rd.read()

Nx = k.shape[0]-2
Ny = k.shape[1]-2
#k = np.ones((Ny+2,Nx+2))
#
# Crea la malla del dominio
#
Lx, Ly, hx, hy, x, y = ms.createMesh(a,b,c,d,Nx,Ny)
#
# Prepara el sistema de ecuaciones a resolver
#
ht = 0.001
Tmax = 1.0
Nmax = int(Tmax/ht)
r = ht /(hx * hy)
f = fd.initRHS(Nx,Ny)
u = fd.initSol(Nx,Ny)
fd.dirichlet(f,u,bA,bB,bC,bD)
A = fd.laplacian2D_stat(Nx,Ny,r,k,1)
#
# Condición inicial
#
u[:] = np.copy(k[:])
#
# Fuente
#
si = int(Nx/2)
sj = int(Ny/2)
source = 1000.0 * ht
f[si,sj] = source
#source = 0.0

for i in range(Nmax):
#
# Resuelve el sistema lineal
#
	f = np.copy(u[1:Ny+1,1:Nx+1])
	f[si,sj] = source
	sol.solve(A,u,f)
	fd.dirichlet(f,u,bA,bB,bC,bD)
#
# Visualiza la solución
#
#	vis.colormap(u,'inferno') #hot, cool, rainbow, ...
	if (i % 20) == 0:
#		vis.surface(x,y,u,'inferno')
		vis.contour(x,y,u,'inferno')
#		vis.colormap(u,'inferno') #hot, cool, rainbow, ...

#
# Escribe datos a un archivo
#
#vis.write(x,y,u,ofilename)
