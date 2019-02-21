#!/usr/bin/env pyhton

import numpy as np
import poisson2D 
import sys
import time

# Los datos del problema se leen de un archivo. Se debe pasar
# el nombre del archivo de entrada; se pasa tambien el nombre de
# un archivo de salida donde se almacenara la solucion.
try:
    in_file_name = sys.argv[1]; out_file_name = sys.argv[2]
except:
    mensaje = """ Error: La ejecucion de este programa requiere de 2 argumentos. 
        Ejecucion correcta: python 2D_Poisson_01.py entrada salida
        donde "entrada" es el nombre de un archivo que contiene los
        datos del problema :  a? y b? son las coordenadas 
        inicial y final del dominio respectivamente las direcciones 
        x y y; Nx y Ny son el numero de incognitas en las 
        direcciones correspondientes; A, B, C y D son las 
        condiciones de frontera. El nombre "salida" se usa para 
        almacenar la solucion del problema.

        Por ejemplo: python 2D_Poisson_01.py INPUT_03 SALIDA"""

    print(mensaje)
    sys.exit(1)
        
ax, bx, ay, by, Nx, Ny, boundA, boundB, boundC, boundD = poisson2D.LeeDatos(in_file_name)
        
# Calcula Delta x y Delta y
hx = (bx-ax)/(Nx+1)
hy = (by-ay)/(Ny+1)

poisson2D.ImprimeDatos(ax,bx,ay,by,Nx,Ny,hx,hy,
                       boundA,"Dirichlet",
                       boundB,"Dirichlet",
                       boundC,"Dirichlet",
                       boundC,"Dirichlet")

# Definicion del sistema lineal de N+1 x N+1
f = np.zeros((Ny,Nx)) # RHS
A = poisson2D.Laplaciano2D(Nx, Ny,-4) # Matriz del sistema

# Aplicacion de las condiciones de frontera Dirichlet
f[0   ,:] -= boundA # Bottom wall    
f[Ny-1,:] -= boundB # Upper wall
f[:,0   ] -= boundC # Left wall 
f[:,Nx-1] -= boundD # Right wall

# La solucion sera guardada en el arreglo u, que es de tamanio Ny+2 x Nx+2, pues incluye las fronteras
u = np.zeros((Ny+2, Nx+2))

# Se utiliza un algoritmo del paquete linalg para obtener la solucion del sistema de N x N
ut = np.copy(u[1:Ny+1,1:Nx+1])
ut.shape = ut.size   # Cambiamos los arreglos a formato unidimensional
f.shape = f.size     # Cambiamos los arreglos a formato unidimensional

t1_start = time.perf_counter()
ut = np.linalg.solve(A,f)
t1_stop = time.perf_counter()
print(time.ctime(), '\n CPUT time: {:0.6f} \n '.format(t1_stop-t1_start))

# Los valores en los lados del dominio son conocidos debido a las cond. Dirichlet
u[Ny+1,:   ] = boundB 
u[:   ,0   ] = boundC 
u[:   ,Nx+1] = boundD
u[0   ,:   ] = boundA 

poisson2D.ImprimeSistema(A,ut,f)


# Calcula el error con respecto a la solucion analitica

input('\n Presione <enter> para generar la grafica ')

ut.shape = (Ny, Nx) # Regresamos el arreglo a formato bidimensional
u[1:Ny+1,1:Nx+1] = ut

x = np.linspace(ax,bx,Nx+2)
y = np.linspace(ay,by,Ny+2)
xg, yg = np.meshgrid(x,y)

poisson2D.GuardaSolucion(out_file_name, x, y, u)

poisson2D.GraficaSuperficieC(xg,yg,u,'hot') #hot, cool, rainbow, ...




