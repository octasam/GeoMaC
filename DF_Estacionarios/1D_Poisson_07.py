#!/usr/bin/env pyhton

import numpy as np
import poisson1D 
import sys

def Evalua(x):
    return np.exp(x) - x - np.exp(1) + 4

# Los datos del problema se leen de un archivo. Se debe pasar
# el nombre del archivo de entrada; se pasa tambien el nombre de
# un archivo de salida donde se almacenara la solucion.
try:
    in_file_name = sys.argv[1]
    out_file_name = sys.argv[2]
except:
    mensaje = """ Error: La ejecucion de este programa requiere de 2 argumentos. 
        Ejecucion correcta: python 1D_Poisson_07.py entrada salida
        donde "entrada" es el nombre de un archivo que contiene los
        datos del problema : a y b, que son las cotas del dominio;
        N el numero de incognitas; A y B las cond. de frontera 
        Dirichlet. El nombre "salida" se usa para almacenar la 
        solucion del problema.

        Por ejemplo: python 1D_Poisson_07.py INPUT_02 SALIDA"""
    print(mensaje)
    sys.exit(1)
        
a, b, N, boundA, boundB = poisson1D.LeeDatos(in_file_name)
        
# Calcula la Delta x
h = (b-a)/(N+1)

poisson1D.ImprimeDatos(a,b,N,h,boundA,"Neumman",boundB,"Dirichlet")

x = np.linspace(a,b,N+2)

# Definicion del sistema lineal de N+1 x N+1
f = np.zeros(N+1)             
f = h*h*np.exp(x[0:N+1])         # Lado derecho
A = poisson1D.Laplaciano1D(N+1,-2) # Matriz del sistema

# Aplicacion de las condiciones de frontera
f[0] = h*boundA   # Neumman
f[N] -= boundB    # Dirichlet
A[0,0] = -1; A[0,1] = 1; # Ajuste de la matrix debido a Neumman

# La solucion sera guardada en el arreglo u, que es de tamanio N+2, pues incluye las fronteras
u = np.zeros(N+2)

# Se utiliza un algoritmo del paquete linalg para obtener la solucion del sistema de N x N
u[0:N+1] = np.linalg.solve(A,f)

# El valor en el extremo derecho es conocido debido a la cond. Dirichlet
u[N+1] = boundB 

poisson1D.ImprimeSistema(A,u,f)

# Calcula el error con respecto a la solucion analitica

Error = np.sqrt(h) * np.linalg.norm(Evalua(x) - u)
print(" Error = %12.10g " % Error)

input('\n Presione <enter> para generar la grafica ')

xa = np.linspace(a,b,100)
ua = Evalua(xa)

GTit = 'Eq. de Poisson : $\partial^2 u(x)/\partial x^2 = e^x, \,\,\, u^\prime(0) = 0; \,\,\, u(1) = 3$'
NSolTit = 'Sol. Num.: $E(h=%2.5g) = %2.5g$' % (h,Error)
ESolTit = 'Sol. Exa.: $u(x) = e^x-x+4-e$'
poisson1D.GraficaSolucion(x,u,xa,ua,GTit,NSolTit,ESolTit,"1D_Poisson_07.pdf")

poisson1D.GuardaSolucion(out_file_name, x, u)



