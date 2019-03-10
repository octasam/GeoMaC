import read_data as rd
import meshes as ms
import phys_data as pd 
import finite_differences as fd 
import solvers as sol
import visflow as vis
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
#
# Crea la malla del dominio
#
Lx, Ly, hx, hy, x, y = ms.createMesh(a,b,c,d,Nx,Ny)
#
# Prepara el sistema de ecuaciones a resolver
#
f = fd.initRHS(Nx,Ny)
u = fd.initSol(Nx,Ny)
fd.dirichlet(f,u,bA,bB,bC,bD)
A = fd.laplacian2D_stat(Nx,Ny,1,k)
#
# Resuelve el sistema lineal
#
sol.solve(A,u,f)
#
# Visualiza la solución
#
vis.colormap(u,'inferno') #hot, cool, rainbow, ...
vis.surface(x,y,u,'inferno','superficie.png')
vis.contour(x,y,u,'inferno','solucion.pdf')
#
# Escribe datos a un archivo
#
vis.write(x,y,u,ofilename)
