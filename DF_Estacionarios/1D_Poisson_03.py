#!/usr/bin/env pyhton
"""
    Construye y resuelve el sistema Ax = b, producto de aplicar \n
    el método de diferencias finitas a la ecuación de Poisson   \n
    en una dimensión, con una condicion tipo Dirichlet(lado izquierdo)
    y una condición tipo Neumman(lado derecho).

Ejecución:

    python 1D_Poisson_03.py a b N A B
    donde a y b son las fronteras izq. y der. del dominio 
    respetivamente N es el numero de incognitas, A es una cond. 
    Dirichlet en a, mientras que B es una cond. Neumman en b.
    Por ejemplo: python 1D_Poisson_03.py 0 1 10 1 0.1
"""


import numpy as np

def Laplaciano1D(N, diagonal):
    """
        Devuelve la matriz A, del sistema Ax = b, 
        considerando las condiciones de frontera. 
        Params:
            N: int
                Número de incognitas. \n
            diagonal: escalar
                Valor de la diagonal principal. \n
        Return:
            A: ndarray
                Matriz NxN, con tres bandas distintas a cero.
    """
    A = np.zeros((N,N))
    A[0,0] = diagonal; A[0,1] = 1
    for i in range(1,N-1):
        A[i,i] = diagonal
        A[i,i+1] = 1
        A[i,i-1] = 1
    A[N-1,N-2] = -1; A[N-1,N-1] = 1 
    return A

def ImprimeDatos(a,b,N,h,A,B):
    """ 
        Muestra en pantalla los datos del problema.
        Params:
            a: escalar
                Punto inicial del dominio
            b: escalar
                Punto final del dominio
            N: int
                Número de nodos(incognitas)
            A: escalar
                Condicion de frontera tipo Dirichlet
            B: escalar
                Condicion de frontera tipo Neumman
    """
    print()
    print("+----------------------------------------------------+")
    print("|      Solucion de la ecuacion de Laplace en 1D      |")
    print("+----------------------------------------------------+")
    print("| Autor: Luis M. de la Cruz S.                       |")
    print("+----------------------------------------------------+")
    print("|                 Datos de entrada                   |")
    print("+----------------------------------------------------+")
    print("| Punto inicial del dominio      : a = %g" % a)
    print("| Punto final del dominio        : b = %g" % b)
    print("| Numero total de incognitas     : N = %d" % N)
    print("| El tamanio de la malla es      : h = %g " % h)
    print("| Cond. de front. Dirichlet en a : A = %g" % A)
    print("| Cond. de front. Neumman en   b : B = %g" % B)
    print("+----------------------------------------------------+")

def ImprimeSistema(A,u,f):
    """
        Muestra en pantalla el sistema lineal de ecuaciones, su 
        solución si y sólo si su longitud es menor igual a 10.
        Params:
            A: ndarray
                Matriz del laplaciano.
            u: ndarray
                Solución del sistema de ecuaciones.
            f: ndarray
                Lado derecho del sistema.

    """
    print("\n Lado derecho del sistema : \n", f)
    print("\n Matriz del sistema : \n", A)
    print("\n Solucion del sistema : \n", u)

# Los datos del problema se leen de la linea de comando como cadenas de caracteres, luego estas se transforman al formato requerido. Se usa el modulo sys.
import sys
try:
    a = float(sys.argv[1])
    b = float(sys.argv[2])
    N = int(sys.argv[3])
    boundA = float(sys.argv[4])
    boundB = float(sys.argv[5])
except:
    mensaje = """ Error: La ejecucion de este programa requiere de 5 argumentos. 
        Ejecucion correcta: python 1D_Poisson_03.py a b N A B
        donde a y b son las fronteras izq. y der. del dominio 
        respetivamente N es el numero de incognitas, A es una cond. 
        Dirichlet en a, mientras que B es una cond. Neumman en b.

        Por ejemplo: python Poisson1D_02.py 0 1 10 1 0.1 """
    print(mensaje)
    sys.exit(1)
    
# Calcula la Delta x
h = (b-a)/(N+1)
ImprimeDatos(a,b,N,h,boundA,boundB)

# Definicion del sistema lineal de N x N
f = np.zeros(N+1)         # Lado derecho
A = Laplaciano1D(N+1, -2) # Matriz del sistema

# Aplicacion de las condiciones de frontera
f[0] -= boundA     # Dirichlet
f[N] = h*boundB # Neumman

# La solucion sera guardada en el arreglo u, que es de tamanio N+2, pues incluye las fronteras
u = np.zeros(N+2)

# Se utiliza un algoritmo del paquete linalg para obtener la solucion del sistema de N x N
u[1:N+2] = np.linalg.solve(A,f)

# El valor en el extremo izquierdo es conocido debido a las cond. Dirichlet
u[0] = boundA

ImprimeSistema(A,u,f)

input('\n Presione <enter> para generar la grafica ')

import matplotlib.pyplot as pp

x = np.linspace(a,b,N+2)
pp.plot(x,u,'-bo')
pp.grid()
pp.show()


