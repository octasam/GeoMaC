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
N = 50
h = (b-a)/(N+1)

# Datos del tiempo de simulacion y paso de tiempo
dt = 0.0001
Tmax = 1.0
Nt = int(Tmax / dt)

# Dato fisico (P. ej. difusividad termica)
k = 1

print(" h = ", h)
print(" dx^2/2k = ", h*h/(2.0*k))
print(" dt = ", dt, " Nt = ", Nt)

# Definicion de los arreglos donde se guardara la solucion
x = np.linspace(a,b,N+2)
u = np.zeros(N+2)

# Condiciones de frontera
boundA = -1 # Dirichlet en a
boundB = 1  # Dirichlet en b
u[0] = boundA
u[N+1] = boundB

# Graficación de la condición inicial
plt.plot(x,u,'--k',label='Inicial')

sumat = 0.0
r = dt * k / (h*h)
error = []

# Ciclo en el tiempo, desde 0 hasta Nt-1
for n in range(Nt):
    err = 0.0  
    
    # Ciclo para resolver en el espacio
    # Euler hacia adelante: EXPLICITO
    t1_start = time.perf_counter()    
    for i in range(1,N+1):
        unew = u[i] + r * (u[i-1] - 2 * u[i] + u[i+1])
        err += (unew - u[i])**2 
        u[i] = unew
    t1_stop = time.perf_counter()
    sumat += (t1_stop - t1_start)

    err = np.sqrt(h*err)
    error.append(err)
    print("n = ", n, ' Error = %12.10g' % err)
    if ((n % 200)==0):
        plt.plot(x,u,'-')
    if (err < tolerancia):
        break

mensaje = 'Sol Num. Error = %5.3g, Pasos = %d, CPU = %g segs' % (err,n,sumat)
plt.plot(x,u,'--k',lw=3,label='Final')
plt.xlabel('$x$')
plt.ylabel('$u(x)$')
plt.grid()
plt.title(mensaje)
plt.legend()
plt.show()

plt.plot(error)
plt.yscale('log')
plt.xlabel('$n$')
plt.ylabel('$RMS$')
plt.grid()
plt.show()


