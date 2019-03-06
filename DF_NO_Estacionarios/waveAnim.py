#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sat Sep 16 17:03:18 2017

@author: luiggi
"""

import numpy as np
import matplotlib.pyplot as plt

def f(x):
    """
    Calcula la condición inicial del problema en todo el dominio.
    x : arreglo que contiene las coordenadas de la malla.
    """
    return np.sin(np.pi * x)

def g(x):
    """
    Calcula la velocidad inicial del problema.
    x : arreglo que contiene las coordenadas de la malla.
    """
    return 0

def solExacta(x, t):
    """
    Calcula la solución exacta del problema.
    x : arreglo que contiene las coordenadas de la malla.
    t : instante de tiempo en que se calcula la solución exacta.
    """
    return np.sin(np.pi * x) * np.cos(2 * np.pi * t)

def calcError(sol_n,sol_e):
    """
    Calcula el error de la solución aproximada con respecto a la exacta.
    sol_n : arreglo que contiene la solución numérica.
    sol_e : arreglo que contiene la solución exacta.
    """
    return np.abs(sol_n-sol_e)

def condicionesIniciales(l, ht, u, x):
    """
    Calcula la condición inicial en el instante n = 1
    l : representa a lambda = alpha * ht / h.
    ht : paso de tiempo.
    u : condición inicial en el instante correspondiente a n = 0.
    x : arreglo que contiene las coordenadas de la malla.
    """
    N = len(u)
    w = np.zeros(N)
    for i in range(1,N-1):
        w[i] = (1 - l**2) * u[i] + 0.5 * l**2 * (u[i+1] + u[i-1]) + ht * g(x[i])
#        w[i] = u[i] + ht * g(x[i])
    return w

def solver(j):
    """
    Calcula la solución del problema para n = 1, 0, ..., Nt
    u : condición inicial en el instante correspondiente a n = 0.
    w : condición inicial en el instante correspondiente a n = 1
    N : Número total de incógnitas internas.
    x : arreglo que contiene las coordenadas de la malla.
    Nt : Número total de pasos en el tiempo.
    l : representa a lambda = alpha * ht / h.
    """
    global u, w
    global ht, lamb
    
    time_step = j * ht
    print('time step = {}'.format(time_step))
    
    s = np.zeros(N+2)
    for i in range(1,N+1):
        s[i] = 2 * (1 - lamb**2) * w[i] + lamb**2 * (w[i+1] + w[i-1]) - u[i]
    u = w.copy()
    w = s.copy()
        
    line.set_ydata(s) # cambia los datos en la dirección y
    label.set_text('Step = {:>5d} \n Time = {:>3.2f}'.format(j, time_step))
    ax.set_title(title_graf)

#
# Inicializamos algunos módulos para la animación
#
from matplotlib.animation import FuncAnimation
fig = plt.figure(figsize=(8,4))           # Figuras

title_graf = '$\partial^2 u / \partial t^2 - a^2 (\partial^2 u /\partial x^2) = 0$ con FDM'

L = 1        # Longitud del dominio
N = 9        # Número de incógnitas internas 9
Tmax = 1.0   # Tiempo máximo de simulación
ht = 0.05    # Paso de tiempo .005
alpha = 2    # Dato físico

h = L / (N+1)        # Tamaño de la malla
Nt = int(Tmax / ht)  # Numero total de pasos
lamb = alpha * ht / h     # Parámetro lambda
Tmax = Nt * ht       # Tiempo total de simulación

print(" N = %d, h = %g, lambda = %g, Nt = %g, ht = %g, Tmax = %g" % (N,h,lamb,Nt,ht,Tmax))

x = np.linspace(0,L,N+2)  # Definición de la malla

u = f(x) # Condición inicial 
w = condicionesIniciales(lamb, ht, u, x)

ax = plt.axes(xlim=(-0.2, L + 0.2), ylim=(-1.2, 1.2)) # Ejes
line, = ax.plot(x, u, '--ko', lw=2, label='FDM')
label = ax.text(-.05, 0.97, 'Time = {:>3.2f}'.format(0),
                ha='center', va='center',
                fontsize=12)
ax.plot(x, solExacta(x,Tmax),'^b-', label = "Exacta")  
plt.xlabel('$x$ [m]')
plt.ylabel('$u$ [...]')
plt.grid()
ax.legend()

#
# Función que controla todo el cálculo y la animación
#
anim = FuncAnimation(fig,            # La figura
                     solver, # la función que cambia los datos
                     interval=500,     # Intervalo entre cuadros en milisegundos
                     frames=Nt+1,   # Cuadros
                     repeat=True)   # Permite poner la animación en un ciclo

plt.show()