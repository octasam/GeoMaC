#!/usr/bin/env pyhton

import numpy as np
import poisson2D
import time

import matplotlib.pyplot as plt
#from mpl_toolkits.mplot3d import Axes3D
from matplotlib.animation import FuncAnimation

filename = input('Nombre del archivo:')
k = poisson2D.leeImagen(filename,0.075)

ax = 0.0
bx = 1.0
ay = 0.0
by = 1.0
Nx = k.shape[0]-2
Ny = k.shape[1]-2
boundA = -1
boundB = -1
boundC = 1
boundD = 1

# TESTING ...
#Nx = 20
#Ny = 20
#k = np.ones((Nx+2,Ny+2))

#print(k)
## Calcula Delta x y Delta y
hx = (bx-ax)/(Nx+1)
hy = (by-ay)/(Ny+1)

ht = 0.0001
Tmax = .01
Nt = int(Tmax / ht)

r = ht / (hx*hy)

poisson2D.ImprimeDatos(ax,bx,ay,by,Nx,Ny,hx,hy,
                       boundA,"Dirichlet",
                       boundB,"Dirichlet",
                       boundC,"Dirichlet",
                       boundD,"Dirichlet")

## Definicion del sistema lineal de N+1 x N+1
f = np.zeros((Ny,Nx)) # RHS
A = poisson2D.Laplaciano2D_T(Nx, Ny,r,k) # Matriz del sistema


#for j in range(Ny):
#	for i in range(Nx):
#		f[i,j] = ht * k[i+1,j+1]

## Aplicacion de las condiciones de frontera Dirichlet

f[Ny-1,:   ] += boundB # Top wall
f[0   ,:   ] += boundA # Bot wall
f[:   ,0   ] += boundC # Left wall
f[:   ,Nx-1] += boundD # Right wall

#print(f)

## La solucion sera guardada en el arreglo u, que es de tamanio Ny+2 x Nx+2, pues incluye las fronteras
u = np.zeros((Ny+2, Nx+2))
ut = np.copy(u[1:Ny+1,1:Nx+1])

## Los valores en los lados del dominio son conocidos debido a las cond. Dirichlet
u[Ny+1,:   ] = boundB # Top wall
u[0   ,:   ] = boundA # Bot wall
u[:   ,0   ] = boundC # Left wall
u[:   ,Nx+1] = boundD # Right wall

x = np.linspace(ax,bx,Nx+2)
y = np.linspace(ay,by,Ny+2)
xg, yg = np.meshgrid(x,y)

fig = plt.figure(figsize=(4,4))           # Figuras
ax = plt.axes(xlim=(0, 1), ylim=(0, 1)) # Ejes

print('Solution starting ...')
for i in range(Nt):
	print(i,sep=' ',end=' ')

## Se utiliza un algoritmo del paquete linalg para obtener la solucion del sistema de N x N
	ut.shape = ut.size   # Cambiamos los arreglos a formato unidimensional
	f.shape = f.size     # Cambiamos los arreglos a formato unidimensional

	f += ut

#	t1_start = time.perf_counter()
	ut = np.linalg.solve(A,f)
#	t1_stop = time.perf_counter()
#	print(time.ctime(), '\n CPU time: {:0.6f} '.format(t1_stop-t1_start))

	ut.shape = (Ny, Nx) # Regresamos el arreglo a formato bidimensional
	u[1:Ny+1,1:Nx+1] = ut

	f.shape = (Ny, Nx) # Regresamos el arreglo a formato bidimensional
	f[:] = 0
	f[Ny-1,:   ] += r * boundB # Top wall
	f[0   ,:   ] += r * boundA # Bot wall
	f[:   ,0   ] += r * boundC # Left wall
	f[:   ,Nx-1] += r * boundD # Right wall

	if( (i % 5)==0):
		ax.contourf(x, y, u, 30, cmap='inferno')
#		ax.contour(x, y, u, 10, colors='black', alpha=0.01)
		plt.savefig('hot%s.png' % i)

#plt.show()

