#!/usr/bin/env python
# coding: utf-8

# # Cálculo de derivadas numéricas
# ### Autor: Luis M. de la Cruz Salas
# ### Rev: 16/02/2019


import numpy as np

# * Aproximación de la primera derivada usando diferencias finitas hacía adelante (Forward):
# 
# $\displaystyle
# \frac{\partial u(x)}{\partial x} \approx \lim\limits_{h\to 0} \frac{u(x+h) - u(x)}{h}
# $

def forwardFD(u,x,h):
    """ Esquema de diferencias finitas hacia adelante
    u : es la función a evaluar
    x : es la posición
    h : es el tamaño de malla.
    """
    return (u(x+h)-u(x))/h

# * Aproximación de la primera derivada usando diferencias finitas hacía atrás (Backward):
# 
# $\displaystyle
# \frac{\partial u(x)}{\partial x} \approx \lim\limits_{h\to 0} \frac{u(x) - u(x-h)}{h}
# $

def backwardFD(u,x,h):
    """ Esquema de diferencias finitas hacia atrás
    u : es la función a evaluar
    x : es la posición
    h : es el tamaño de malla.
    """
    return (u(x)-u(x-h))/h


# * Aproximación de la primera derivada usando diferencias finitas hacía centradas (Centered):
# 
# $\displaystyle
# \frac{\partial u(x)}{\partial x} \approx \lim\limits_{h\to 0} \frac{u(x+h) - u(x-h)}{2h}
# $

def centeredFD(u,x,h):
    """ Esquema de diferencias finitas centradas
    u : es la función a evaluar
    x : es la posición
    h : es el tamaño de malla.
    """
    return (u(x+h)-u(x-h))/(2*h)


def error(ue, un):
    """ Cálculo del error
    ue : valor exacto
    un : valor numérico aprximado
    """
    return np.fabs(un - ue)

