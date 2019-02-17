#!/usr/bin/env pyhton
"""
    Construye y resuelve el sistema Ax = b, producto de aplicar \n
    el método de diferencias finitas a la ecuación de Poisson   \n
    en una dimensión, con condiciones de frontera tipo Dirchelet.

Ejecución:

    Solicita del usuario las dimensiones del dominio, \n
    el numero de nodos y el valor de las condiciones.
"""
import numpy as np

def Laplaciano1D(N, diagonal):
    """
        Devuelve la matriz A, del sistema Ax = b. 
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
    A[N-1,N-2] = 1; A[N-1,N-1] = diagonal
    return A

# Comienza la ejecución del programa
print()
print("+----------------------------------------------------+")
print("|      Solucion de la ecuacion de Laplace en 1D      |")
print("+----------------------------------------------------+")
print("| Autor: Luis M. de la Cruz S.                       |")
print("+----------------------------------------------------+")
print("|                 Datos de entrada                   |")
print("+----------------------------------------------------+")

# Tamanio del dominio (se lee del teclado, se toma como cadena, luego se transforma a flotante)
a = float(input('| Punto inicial del dominio      : a = '))
b = float(input('| Punto final del dominio        : b = '))

# Numero de incognitas y delta x
N = int(input('| Numero total de incognitas     : N = '))
h = (b-a)/(N+1)
print("| El tamanio de la malla es      : h = %g " % h)

# Condiciones Dirichlet
boundA = float(input('| Cond. de front. Dirichlet en a : A = '))
boundB = float(input('| Cond. de front. Dirichlet en b : B = '))
print("+----------------------------------------------------+")

# Definicion del sistema lineal de N x N
f = np.zeros(N)         # Lado derecho
A = Laplaciano1D(N, -2) # Matriz del sistema

# Aplicacion de las cond. Dirichlet para el calculo correcto de la sol.
f[0] -= boundA
f[N-1] -= boundB 

# La solucion sera guardada en el arreglo u, que es de tamanio N+2, pues incluye las fronteras
u = np.zeros(N+2)

# Se utiliza un algoritmo del paquete linalg para obtener la solucion del sistema de N x N
u[1:N+1] = np.linalg.solve(A,f)

# Los valores en los extremos son conocidos debido a las cond. Dirichlet
u[0] = boundA
u[N+1] = boundB 

print("\n Lado derecho del sistema : \n", f)
print("\n Matriz del sistema : \n", A)
print("\n Solucion del sistema : \n", u)

input('\n Presione <enter> para generar la grafica ')

import matplotlib.pyplot as plt

x = np.linspace(a,b,N+2)
plt.plot(x,u,'-bo')
plt.show()


