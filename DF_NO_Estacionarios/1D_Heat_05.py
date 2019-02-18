#!/usr/bin/env pyhton

#
# Autor : Luis M. de la Cruz Salas
# Fecha : Sat Feb 16 20:36:42 CST 2019
#

import numpy as np
import matplotlib.pyplot as plt
from matplotlib.animation import FuncAnimation

import time

# Tolerancia : criterio de termino anticipado
tolerancia = 1e-6

# Datos del dominio y la malla espacial
a = 0.0 
b = 1.0
N = 10
h = (b-a)/(N+1)

# Datos del tiempo de simulacion y paso de tiempo
ht = 0.001
Tmax = 2.0
Nt = int(Tmax / ht)

print('Nt =', Nt)

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

sumat = 0.0
r = ht / (h*h) 

fig = plt.figure(figsize=(8,4))           # Figuras
ax = plt.axes(xlim=(0, 1), ylim=(-1.2, 1.2)) # Ejes
line, = ax.plot(x, u, 'o-b', label='FDM')
label = ax.text(0.75, -.75, 'Time = {:>8.5f}'.format(0),
                ha='center', va='center',
                fontsize=12)

ax.plot(x,u,'--k',label='Inicial',lw=2)
plt.xlabel('$x$')
plt.ylabel('$T$')
plt.grid()
ax.legend()

title_graf = '$\partial T / \partial t = \partial (\Gamma \partial T/\partial x)/\partial x$'

def explicitEuler(i):
    time_step = i * ht
    err = 0.0  
    for i in range(1,N+1): # Ciclo espacial
        unew = u[i] + r * (k[i] * u[i-1] - (k[i+1]+k[i]) * u[i] + k[i+1] * u[i+1])
        err += (unew - u[i])**2 
        u[i] = unew
    err = np.sqrt(h*err)
#
# Actualizamos la gráfica
#
    line.set_ydata(u) # cambia los datos en la dirección y

    if(err < tolerancia):
        label.set_text('Time = {:>8.5e}\n Time = {:>8.5f}'.format(sumat))
    label.set_text('Error = {:>8.5e}\n Time = {:>8.5f}'.format(err, time_step))
    ax.set_title(title_graf)

#
# Función que controla todo el cálculo y la animación
#
anim = FuncAnimation(fig,            # La figura
                     explicitEuler, # la función que cambia los datos
                     interval=1,     # Intervalo entre cuadros en milisegundos
                     frames=Nt+1,   # Cuadros
                     repeat=False)   # Permite poner la animación en un ciclo


plt.show()



