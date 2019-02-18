#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sun Sep 17 15:45:44 2017

@author: luiggi
"""

import numpy as np


def Euler(a,b,f,alpha,Nt):
    """
    Esta funcion resuelve el problema 
    u'=f(t,u), u(0) = alpha en a <= t <=  b 
    usando el método de Euler hacia adelante
    """
    w = np.zeros(Nt+1)
    t = np.zeros(Nt+1)
    ht = (b - a)/Nt
    t[0] = a
    w[0] = alpha
    print(' w[t = ',0.0,'] = w[',0,'] = ', w[0])
    for i in range(Nt):
        w[i+1] = w[i] + ht * f(t[i], w[i])
        t[i+1] = a + (i+1) * ht
        print(' w[t = ',t[i+1],'] = w[',i+1,'] = ', w[i+1])
        
    return t, w, ht

def f(t,y):
    return y - t**2 +1

def Exacta(t):
    return (t+1)**2 - 0.5 * np.exp(t)

Nt = 4 # Valores a usar en el ejemplo: 4 y 10
a = 0
b = 2
alpha = 0.5

t,w,ht = Euler(a,b,f,alpha,Nt)

y = Exacta(t)
import matplotlib.pyplot as plt
plt.plot(t,w,'o-', label = 'Numérica')
plt.plot(t,y,'s-', label = 'Exacta')
plt.grid()
plt.xlabel('$t$')
plt.ylabel('$y(t)$')
plt.legend()
plt.show()

input('Presiona enter para continuar')

Nt = 10 # Valores a usar en el ejemplo: 4 y 10

t,w,ht = Euler(a,b,f,alpha,Nt)
y = Exacta(t)
import matplotlib.pyplot as plt
plt.plot(t,w,'o-', label = 'Numérica')
plt.plot(t,y,'s-', label = 'Exacta')
plt.grid()
plt.xlabel('$t$')
plt.ylabel('$y(t)$')
plt.legend()
plt.show()

input('Presiona enter para continuar')

e = np.abs(y - w)

from pandas import DataFrame

#A = np.array((t,w,y,e))

#tabla = DataFrame(A.T,columns=list('twye'))
tabla = DataFrame(np.array((t,w,y,e)).T,columns=['t_i','w_i','y_i','|y_i-w_i|'])
tabla.to_latex('tabla.tex')
tabla.to_csv('tabla.csv')
tabla.to_excel('tabla.xlsx')
print(tabla)

input('Presiona enter para continuar')

plt.plot(t,e,'vr-')
plt.grid()
plt.xlabel('$t$')
plt.ylabel('$error$')
plt.show()

