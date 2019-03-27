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
k = pd.readImage('TUX50.png',0.1)
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
Tmax = .25
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
val_source = 0.0
source = val_source * ht
f[si,sj] = source

import matplotlib.pyplot as plt
from matplotlib.animation import FuncAnimation
plt.style.use('ggplot')
#
# Figura y ejes
#
fig = plt.figure(figsize=(5,4))           # Figuras
ax = plt.axes(xlim=(0, 1), ylim=(0, 1)) # Ejes
#
# Se dibuja el primer conjunto de datos usando contornos
#
contour_opts = {'levels': np.linspace(0, 1, 20), 'cmap':'inferno', 'lw': 2}
cax = ax.contourf(x, y, u, **contour_opts)
fig.colorbar(cax)

title_graf = '$\partial T / \partial t = \partial (\Gamma \partial T/\partial x)/\partial x$ Time = {:>8.5f}'

plt.show()

input('Presiona <enter> para comenzar la ejecución')

def animate(i):
#
# Resuelve el sistema lineal
#
#	global source
	time_step = i * ht
#	if time_step > Tmax / 10:
#		source = 1000 * ht
	f = np.copy(u[1:Ny+1,1:Nx+1])
	f[si,sj] = source
	sol.solve(A,u,f)
	fd.dirichlet(f,u,bA,bB,bC,bD)
#
# Visualiza la solución
#
	ax.collections = []
	ax.contourf(x, y, u, **contour_opts)
#	label.set_text('Time = {:>8.5f}'.format(time_step))
	ax.set_title(title_graf.format(time_step))
#
# Función que controla la animación
#
anim = FuncAnimation(fig,           # La figura
                     animate,       # la función que cambia los datos
                     interval=100,  # Intervalo entre cuadros en milisegundos
                     frames=Nmax+1,     # Cuadros por segundo
                     repeat=False)   # Permite poner la animación en un ciclo

plt.show()

