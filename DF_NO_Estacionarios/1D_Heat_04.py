#!/usr/bin/env pyhton

#
# Autor : Luis M. de la Cruz Salas
# Fecha : 2015
#

import numpy as np
import matplotlib.pyplot as plt
import time

def Laplaciano(N,r,k):
    A = np.zeros((N,N))
    A[0,0] = 1 + r * (k[1]+k[2]) 
    A[0,1] = -r * k[2]
    for i in range(1,N-1):
        A[i,i] = 1 + r * (k[i+1]+k[i+2])
        A[i,i+1] = -r * k[i+2]
        A[i,i-1] = -r * k[i+1]
    A[N-1,N-2] = -r * k[N]; 
    A[N-1,N-1] = 1 + r * (k[N]+k[N+1]) 
    return A

# Tolerancia : criterio de termino anticipado
tolerancia = 1e-6

# Datos del dominio y la malla espacial
a = 0.0 
b = 1.0
N = 10
h = (b-a)/(N+1)

# Datos del tiempo de simulacion y paso de tiempo
dt = 0.1
Tmax = 2.0
Nt = int(Tmax / dt)

# Definicion de los arreglos donde se guardara la solucion
x = np.linspace(a,b,N+2)
u = np.zeros(N+2)

# Dato fisico (P. ej. difusividad termica)
#Ejemplo 1.
k = np.ones(N+2)

#Ejemplo 2.
k[0:int(N/2)] = 0.02

#Ejemplo 3.
#k[0:int(N/4)] = 1.0
#k[int(N/4):int(N/2)] = 0.5
#k[int(N/2):int(3*N/2)] = 1.5
#k[int(3*N/4):] = 2.0

#Ejemplo 4.
#k = np.random.rand(N+2)

# Condiciones de frontera
boundA = -1 # Dirichlet en a
boundB = 1  # Direchlet en b
u[0] = boundA
u[N+1] = boundB

# Este es lado derecho de la ecuación, que contiene la condicion inicial
# es decir la solucion en el paso 0. Por eso hacemos una copia de u.
f = np.copy(u[1:N+1])
uold = np.copy(u)

# Graficación de la condición inicial
plt.plot(x,u,'o--k',label='Inicial')

# Construccion de la matriz
r = dt / (h*h) 
A = Laplaciano(N,r,k)

sumat = 0.0
error = []

# Ciclo en el tiempo, desde 0 hasta Nt-1
for n in range(Nt+1):
    
    # Solucion en el espacio con Euler hacia atras: IMPLICITO
    t1_start = time.perf_counter()    
    f[0] += r * k[1] * boundA
    f[N-1] += r * k[N+1] * boundB
    u[1:N+1] = np.linalg.solve(A,f) # Sol. del sistema lineal
    t1_stop = time.perf_counter()
    sumat += (t1_stop - t1_start)
    
    err = np.sqrt(h) * np.linalg.norm(uold-u)
    error.append(err)
    
    print("n = ", n, ' Error = %12.10g' % err)
    if ((n % 20)==0):
        plt.plot(x,u,'-')
        
    # Actualizacion de la solucion para dar el siguiente paso
    t1_start = time.perf_counter()    
    f = np.copy(u[1:N+1])
    uold = np.copy(u)
    t1_stop = time.perf_counter()
    sumat += (t1_stop - t1_start)
    
    if (err < tolerancia):
        break

mensaje = 'Sol Num. Error = %5.3g, Pasos = %d, CPU = %g segs' % (err,n,sumat)
plt.plot(x,u,'o--k',lw=3,label='Final')
plt.xlabel('$x$')
plt.ylabel('$u(x)$')
plt.grid()
plt.title(mensaje)
plt.legend()
plt.savefig('implicit.pdf')
plt.show()

plt.plot(x,k,'s-')
plt.ylabel('$\Gamma$')
plt.xlabel('$x$')
plt.grid()
plt.show()

plt.plot(error)
plt.yscale('log')
plt.xlabel('$n$')
plt.ylabel('$RMS$')
plt.grid()
plt.savefig('impError.pdf')
plt.show()
