#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sat Sep 16 17:03:18 2017

@author: luiggi
"""

import numpy as np
import math
import matplotlib.pyplot as plt

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

def RK2(a,b,f,alpha,Nt):
    """
    Esta funcion resuelve el problema 
    u'=f(t,u), u(0) = alpha en a <= t <=  b 
    usando el método Runge-Kutta de orden 2
    """
    w = np.zeros(Nt+1)
    t = np.zeros(Nt+1)
    ht = (b - a)/Nt
    t[0] = a
    w[0] = alpha
    for i in range(Nt):
        k1 = ht * f(t[i], w[i])
        w[i+1] = w[i] + ht * f(t[i] + ht * 0.5, w[i] + k1 * 0.5)
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

t, wE, ht = Euler(a,b,f, alpha,Nt)
t, wT, ht = TaylorO2(a,b,f,alpha,Nt)
t, wR, ht = RK2(a,b,f, alpha,Nt)
y = Exacta(t)

plt.plot(t,wE,'b-', label="Euler")
plt.plot(t,wT,'g-', label="T2")
plt.plot(t,wR,'k-', label="RK2")
plt.plot(t, y, 'r-', label="Exacta")
plt.grid()
plt.legend()
plt.show()

input('Presiona enter para continuar')

eE = np.abs(y - wE)
eT = np.abs(y - wT)
eR = np.abs(y - wR)

from pandas import DataFrame

tabla = DataFrame(np.array((t,wE,wT,wR,y,eE,eT,eR)).T,columns=['t_i','Euler','Taylor2','RK2', 'y','e_Euler','e_Taylor','e_Rk2'])
tabla.to_csv('tabla.csv')
print(tabla)

input('Presiona enter para continuar')

plt.plot(t,eE,'vr-', label='Error Euler')
plt.plot(t,eT,'sb-', label='Error Taylor')
plt.plot(t,eR,'og-', label='Error RK2')
plt.grid()
plt.xlabel('$t$')
plt.ylabel('$error$')
plt.legend()
plt.show()
