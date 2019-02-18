#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sat Sep 16 16:57:20 2017

@author: luiggi
"""
import numpy as np
#import pylab as pl      
import matplotlib.pyplot as pl
    
def sol(t, y0, l):
    return y0 * np.exp(l * t)

def FE(a, b, Nt, y0, l):
    y = np.linspace(a,b,Nt+1)
    ht = (b-a) / Nt    
    y[0] = y0
    for i in range(Nt):
        y[i+1] = (1 + ht * l) * y[i]
    return y

def BE(a, b, Nt, y0, l):
    y = np.linspace(a,b,Nt+1)
    ht = (b-a) / Nt    
    y[0] = y0
    for i in range(Nt):
        y[i+1] = (1 /(1 - ht * l)) * y[i] 
    return y

# Usar diferentes valores de l (lambda) y de Nt para ver la convergencia

Nt = 5
a = 0
b = 2
t = np.linspace(a,b,Nt+1)

y0 = 2
l = -6
y1 = FE(a, b, Nt, y0, l)
y2 = BE(a, b, Nt, y0, l)
ye = sol(t, y0 ,l )

Titulo = ' $y(t) = y_0 e^{\lambda t}$'
pl.plot(t,ye, 'r-', label='Sol. Exa.')
pl.plot(t,y1, '^b-', label='Sol. Num. FE')
pl.plot(t,y2, 'vg-', label='Sol. Num. BE')
pl.title(Titulo)
pl.xlabel('$t$')
pl.ylabel('$y(t)$')
pl.legend()
pl.grid()
pl.show()



