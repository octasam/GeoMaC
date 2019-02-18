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
N = 50
h = (b-a)/(N+1)

# Datos del tiempo de simulacion y paso de tiempo
ht = 0.0001
Tmax = 1.0
Nt = int(Tmax / ht)

# Dato fisico (P. ej. difusividad termica)
k = 1

print(" h = ", h)
print(" dx^2/2k = ", h*h/(2.0*k))
print(" ht = ", ht, " Nt = ", Nt)

# Condiciones de frontera
boundA = -1 # Dirichlet en a
boundB = 1  # Direchlet en b

# Definicion de los arreglos donde se guardar la solucion
x = np.linspace(a,b,N+2)
u = np.zeros(N+2)
u[0] = boundA
u[N+1] = boundB

#plt.plot(x,u,'-b')

sumat = 0.0

r = ht * k / (h*h) 

fig = plt.figure(figsize=(8,4))           # Figuras
ax = plt.axes(xlim=(0, 1), ylim=(-1.2, 1.2)) # Ejes
line, = ax.plot(x, u, '--', label='FDM')
label = ax.text(0.75, -.75, 'Time = {:>8.5f}'.format(0),
                ha='center', va='center',
                fontsize=12)

ax.plot(x,u,'b-',label='Sol. Exac',lw=2)
plt.xlabel('$x$ [m]')
plt.ylabel('$\phi$ [...]')
ax.legend()

title_graf = '$\partial \phi / \partial t + \partial(p u \phi)/\partial x= \partial (\Gamma \partial\phi/\partial x)/\partial x$ con FVM'

def explicitEuler(i):
    time_step = i * ht
    err = 0.0  
    for i in range(1,N+1): # Ciclo espacial
        u_next = u[i] + r * (u[i+1] - 2 * u[i] + u[i-1])
        err += (u_next - u[i])**2 
        u[i] = u_next
    err = np.sqrt(h*err)
#    t1 = time.clock()
#    sumat += (t1 - t0)
#    print("n = ", n, ' Error = %12.10g' % err)
#    if ((i % 100)==0):
#        print('time step = {:.5f} \t error = {.5f}'.format(time_step,err))

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



