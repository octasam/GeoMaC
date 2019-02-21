#!/usr/bin/env pyhton

#
# Autor : Luis M. de la Cruz Salas
# Fecha : lun feb 18 12:39:48 CST 2019
#

import numpy as np
import matplotlib.pyplot as plt
import time

# Tolerancia : criterio de termino anticipado
tolerancia = 1e-6

# Datos del dominio y la malla espacial
a = 0.0 
b = 1.0
N = 10
h = (b-a)/(N+1)

# Datos del tiempo de simulacion y paso de tiempo
dt = 0.001
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
boundB = 1  # Dirichlet en b
u[0] = boundA
u[N+1] = boundB

# Graficación de la condición inicial
plt.plot(x,u,'o--k',label='Inicial')

sumat = 0.0
r = dt / (h*h)
error = []

# Ciclo en el tiempo, desde 0 hasta Nt-1
for n in range(Nt+1):
    err = 0.0  
    
    # Ciclo para resolver en el espacio
    # Euler hacia adelante: EXPLICITO
    t1_start = time.perf_counter()    
    for i in range(1,N+1):
        unew = u[i] + r * (k[i] * u[i-1] - (k[i+1]+k[i]) * u[i] + k[i+1] * u[i+1])
        err += (unew - u[i])**2 
        u[i] = unew
    t1_stop = time.perf_counter()
    sumat += (t1_stop - t1_start)

    err = np.sqrt(h*err)
    error.append(err)
    print("n = ", n, ' Error = %12.10g' % err)
    if ((n % 20)==0):
        plt.plot(x,u,'-')
    if (err < tolerancia):
        break

mensaje = 'Sol Num. Error = %5.3g, Pasos = %d, CPU = %g segs' % (err,n,sumat)
plt.plot(x,u,'o--k',lw=3,label='Final')
plt.xlabel('$x$')
plt.ylabel('$u(x)$')
plt.grid()
plt.title(mensaje)
plt.legend()
plt.savefig('explicit.pdf')
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
plt.savefig('expError.pdf')
plt.show()


