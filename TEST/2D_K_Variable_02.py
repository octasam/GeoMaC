#!/usr/bin/env pyhton

import numpy as np
import poisson2D
import time

import matplotlib.pyplot as plt

filename = input('Nombre del archivo:')
k = poisson2D.leeImagen(filename,0.075)

ax = 0.0
bx = 1.0
ay = 0.0
by = 1.0
Nx = k.shape[0]-2
Ny = k.shape[1]-2
boundA = -0
boundB = -0
boundC = 0
boundD = 0

## Calcula Delta x y Delta y
hx = (bx-ax)/(Nx+1)
hy = (by-ay)/(Ny+1)

ht = 1
r = ht / (hx*hy)

poisson2D.ImprimeDatos(ax,bx,ay,by,Nx,Ny,hx,hy,
                       boundA,"Dirichlet",
                       boundB,"Dirichlet",
                       boundC,"Dirichlet",
                       boundD,"Dirichlet")

## Definicion del sistema lineal de N+1 x N+1
f = np.zeros((Ny,Nx)) # RHS
A = poisson2D.Laplaciano2D(Nx, Ny,1,k) # Matriz del sistema


for j in range(Ny):
	for i in range(Nx):
		f[i,j] = k[i+1,j+1]

## Aplicacion de las condiciones de frontera Dirichlet

f[Ny-1,:   ] += boundB # Top wall
f[0   ,:   ] += boundA # Bot wall
f[:   ,0   ] += boundC # Left wall
f[:   ,Nx-1] += boundD # Right wall


## La solucion sera guardada en el arreglo u, que es de tamanio Ny+2 x Nx+2, pues incluye las fronteras
u = np.zeros((Ny+2, Nx+2))

## Se utiliza un algoritmo del paquete linalg para obtener la solucion del sistema de N x N
ut = np.copy(u[1:Ny+1,1:Nx+1])
ut.shape = ut.size   # Cambiamos los arreglos a formato unidimensional
f.shape = f.size     # Cambiamos los arreglos a formato unidimensional

t1_start = time.perf_counter()
ut = np.linalg.solve(A,f)
t1_stop = time.perf_counter()
print(time.ctime(), '\n CPU time: {:0.6f} '.format(t1_stop-t1_start))

## Los valores en los lados del dominio son conocidos debido a las cond. Dirichlet
u[Ny+1,:   ] = boundB # Top wall
u[0   ,:   ] = boundA # Bot wall
u[:   ,0   ] = boundC # Left wall
u[:   ,Nx+1] = boundD # Right wall

poisson2D.ImprimeSistema(A,ut,f)

ut.shape = (Ny, Nx) # Regresamos el arreglo a formato bidimensional
u[1:Ny+1,1:Nx+1] = ut

x = np.linspace(ax,bx,Nx+2)
y = np.linspace(ay,by,Ny+2)
xg, yg = np.meshgrid(x,y)

poisson2D.GuardaSolucion('imagenX', x, y, k)

plt.imshow(u,cmap='inferno')
plt.show()

#poisson2D.GuardaSolucion('SALIDA', x, y, u)

# Post-procesamiento ...
NNX = u.shape[0]
unew = np.copy(u)
for j in range(u.shape[0]):
    for i in range(u.shape[1]):
        unew[i,j] = u[NNX-i-1,j]

poisson2D.GraficaSuperficieC(xg,yg,unew,'inferno') #hot, cool, rainbow, ...

