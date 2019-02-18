#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Sep 18 13:15:14 2017

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
    for i in range(Nt):
        w[i+1] = w[i] + ht * f(t[i], w[i])
        t[i+1] = a + (i+1) * ht        
    return t, w, ht


def TaylorO2(a,b,f,alpha,Nt):
    """
    Esta funcion resuelve el problema 
    u'=f(t,u), u(0) = alpha en a <= t <=  b 
    usando el método de Taylor de orden 2
    """
    w = np.zeros(Nt+1)
    t = np.zeros(Nt+1)
    ht = (b - a)/Nt
    t[0] = a
    w[0] = alpha
    for i in range(Nt):
        w[i+1] = w[i] + ht * ((1+0.5*ht)*(w[i]-t[i]**2+1) - ht*t[i])
        t[i+1] = a + (i+1) * ht        
    return t, w, ht

def f(t,y):
    return y - t**2 +1

def Exacta(t):
    return (t+1)**2 - 0.5 * np.exp(t)

Nt = 10
a = 0
b = 2
alpha = 0.5

t,w,ht = Euler(a,b,f,alpha,Nt)
t,wT,ht = TaylorO2(a,b,f,alpha,Nt)
y = Exacta(t)

import matplotlib.pyplot as plt
plt.plot(t,w,'o-', label = 'Euler')
plt.plot(t,wT,'v-', label = 'Taylor O2')
plt.plot(t,y,'-', label = 'Exacta')
plt.grid()
plt.xlabel('$t$')
plt.ylabel('$y(t)$')
plt.legend()
plt.show()

input('Presiona enter para continuar')

eE = np.abs(y - w)
eT = np.abs(y - wT)

from pandas import DataFrame

tabla = DataFrame(np.array((t,w,wT,y,eE,eT)).T,columns=['t_i','wE_i','wT_i','y_i','e_Euler','e_Taylor'])
tabla.to_csv('tabla.csv')
print(tabla)

input('Presiona enter para continuar')

plt.plot(t,eE,'vr-', label='Error Euler')
plt.plot(t,eT,'sb-', label='Error Taylor')
plt.grid()
plt.xlabel('$t$')
plt.ylabel('$error$')
plt.legend()
plt.show()

