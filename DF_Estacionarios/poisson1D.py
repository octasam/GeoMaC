#!/usr/bin/env pyhton
"""
    Construye y resuelve el sistema Ax = b, producto de aplicar \n
    el método de diferencias finitas a la ecuación de Poisson   \n
    en una dimensión, con condiciones de frontera tipo Dirchelet.

Example:

    Solicita del usuario las dimensiones del dominio, \n
    el numero de nodos y el valor de las condiciones.
"""



import numpy as np
import matplotlib.pyplot as plt
import time

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

def LeeDatos(filename):
    """
        Lee información de un archivo de texto.
        Param
            filename: string
                Nombre del archivo de texto a ser leido
        Return:
            a: escalar
                Punto inicial del dominio
            b: escalar
                Punto final del dominio
            N: int
                Número de nodos(incognitas)
            A: escalar
                Condicion de frontera  
            B: escalar
                Condicion de frontera  

    """
    ifile = open(filename, 'r')     # abre el archivo de entrada
    file_lines = ifile.readlines()  # lee las lineas del archivo
    ifile.close();                  # cierra el archivo de entrada
    a, b, N, A, B = file_lines[0].split() # separa las columnas de la primera linea
    a = float(a); b = float(b); N = int(N); A = float(A); B = float(B); # convierte los datos
    return a, b, N, A, B

def ImprimeDatos(a,b,N,h,A,cd1,B,cd2):
    """ 
        Muestra en pantalla los datos del problema.
        Params:
            a: escalar
                Punto inicial del dominio
            b: escalar
                Punto final del dominio
            N: int
                Número de nodos(incognitas)
            h: escalar
                Tamanio de la malla
            A: escalar
                Valor de la condicion de frontera  
            cd1: string
                Tipo de condicion nodo a.
            B: escalar
                Condicion de frontera 
            cd2: string
                Tipo de condicion nodo b.
    """
    print()
    print("+----------------------------------------------------+")
    print("|      Solucion de la ecuacion de Laplace en 1D      |")
    print("+----------------------------------------------------+")
    print("| Autor: Luis M. de la Cruz S.                       |")
    print("+----------------------------------------------------+")
    print("|                 Datos de entrada                   |")
    print("+----------------------------------------------------+")
    print("| Punto inicial del dominio        : a = %g" % a)
    print("| Punto final del dominio          : b = %g" % b)
    print("| Numero total de incognitas       : N = %d" % N)
    print("| El tamanio de la malla es        : h = %g " % h)
    print("| Cond. de front. ", cd1, " en a : A = %g" % A)
    print("| Cond. de front. ", cd2, " en b : B = %g" % B)
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

def GraficaSolucion(x,u,xa,ua,Gtitle,Ntitle,Etitle,filename):
    """
        Grafica la solución.
        Param
            x: ndarray
                Coordenada x, de cada nodo de la malla.
            u: ndarray
                Solución aproximada.
            xa: ndarray
                Coordenada x, para la solución exacta.
            ua: ndarray
                Solución exacta.
            Gtitle: string 
                Nombre del grafico.
            Ntitle: string
                Etiqueta para la solución aproximacimada.
            Etitle: string
                Etiqueta para la solución exacta. 
            filename: string
                Nombre de la imagen para ser almacenada.
    """
    umax = np.max(u); umin = np.min(u); 
    umaxa = np.max(ua); umina = np.min(ua); 
    umax = max(umax,umaxa); umin = min(umin,umina); 
    percent = 0.1*np.abs(umax - umin)
    umax += percent; umin -= percent;
    
    plt.plot(x,u,'--o',label=Ntitle)
    plt.plot(xa,ua,'-r',label=Etitle)
    plt.title(Gtitle)
    plt.xlabel("$x$")
    plt.ylabel("$u(x)$")
    plt.ylim(umin,umax)
    plt.legend(loc='upper right')
    plt.savefig(filename)
    plt.show()
    
def GuardaSolucion(filename, x, u):
    """
        Almacena la solución en un archivo de texto.
        Param
            filename: string
                Nombre del archivo de texto donde
                se escribira el resultado.
            x: ndarray
                Coordenada x, de cada nodo de la malla.
            u: ndarray
                Solución obtenida.
    """
    ofile = open(filename, 'w') 
    for i in range(0,len(u)):
        ofile.write('%g\t%12.5e\n' % (x[i],u[i]))
    ofile.close()

if __name__ == "__main__":

    import sys
    # Los datos del problema se leen de un archivo. Se debe pasar
    # el nombre del archivo de entrada; se pasa tambien el nombre de
    # un archivo de salida donde se almacenara la solucion.
    try:
        in_file_name = sys.argv[1]; out_file_name = sys.argv[2]
    except:
        mensaje = """ Error: La ejecucion de este programa requiere de 2 argumentos. 
        Ejecucion correcta: python poisson1D.py entrada salida
        donde "entrada" es el nombre de un archivo que contiene los
        datos del problema : a y b, que son las cotas del dominio;
        N el numero de incognitas; A y B las cond. de frontera 
        Dirichlet. El nombre "salida" se usa para almacenar la 
        solucion del problema.
        Ejemplo de archivo entrada: 0 1 20 1 1 """; print(mensaje); sys.exit(1)
        
    a, b, N, boundA, boundB = LeeDatos(in_file_name)        
        
    # Calcula la Delta x
    h = (b-a)/(N+1)

    ImprimeDatos(a,b,N,h,boundA,"Dirichlet",boundB,"Dirichlet")

    # Definicion del sistema lineal de N x N
    f = np.zeros(N)         # Lado derecho
    A = Laplaciano1D(N, -2) # Matriz del sistema

    # Aplicacion de las condiciones de frontera Dirichlet
    f[0] -= boundA   
    f[N-1] -= boundB 

    # La solucion sera guardada en el arreglo u, que es de tamanio N+2, pues incluye las fronteras
    u = np.zeros(N+2)

    # Se utiliza un algoritmo del paquete linalg para obtener la solucion del sistema de N x N
    u[1:N+1] = np.linalg.solve(A,f)

# Los valores en los extremos son conocidos debido a las cond. Dirichlet
    u[0] = boundA
    u[N+1] = boundB 

    ImprimeSistema(A,u,f)

# Calcula el error con respecto a la solucion analitica
    x = np.linspace(a,b,N+2)
    Error = np.sqrt(h) * np.linalg.norm(x - u)
    print(" Error = %12.10g " % Error)

    input('\n Presione <enter> para generar la grafica ')

    xa = np.linspace(a,b,100)
    ua = xa
    GraficaSolucion(x,u,xa,ua,"Prueba con Laplace","Numerica","Exacta","laplace1D.pdf")

    GuardaSolucion(out_file_name, x, u)



