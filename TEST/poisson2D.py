#!/usr/bin/env pyhton

import numpy as np
from mpl_toolkits.mplot3d import Axes3D
import pylab as pl
import matplotlib.image as mpimg

def leeImagen(filename, factor):
    img=mpimg.imread(filename)
    k = img[:,:,0]
    for j in range(k.shape[1]):
        for i in range(k.shape[0]):
            if k[i,j] == 0.0:
                k[i,j] = 0.1
            if k[i,j] < 1.0:
                k[i,j] *= factor                

    return k.transpose()


def u_face(u1, u2):
    return 0.5 * (u1 + u2)

def Laplaciano2D(Nx, Ny, r, k):
    """ Esta funcion calcula los coeficientes del 
    sistema lineal producido por el operador de 
    Laplace en 2D. Estos coeficientes son almacenados 
    en la matriz pentadiagonal correspondiente."""
    N = Nx * Ny
    A = np.zeros((N,N))

# Primero llena los bloques tridiagonales
    for j in range(1,Ny+1):
        ofs = Nx * (j-1)     
        # Primer renglón del bloque, considera BC en la pared izq.
        k1 = u_face(k[0,j  ], k[1,j]) # k_(i-1/2, j)
        k2 = u_face(k[2,j  ], k[1,j]) # k_(i+1/2, j)
        k3 = u_face(k[1,j+1], k[1,j]) # k_(i, j-1/2)
        k4 = u_face(k[1,j-1], k[1,j]) # k_(i, j+1/2)
        A[ofs    , ofs] = r * (k1 + k2 + k3 + k4) 
        A[ofs + 1, ofs] = -r * k2

        # Renglones intermedios del bloque 
        for i in range(2,Nx):
            k1 = u_face(k[i-1,j], k[i,j]) # k_(i-1/2, j)
            k2 = u_face(k[i+1,j], k[i,j]) # k_(i+1/2, j)
            k3 = u_face(k[i,j-1], k[i,j]) # k_(i, j-1/2)
            k4 = u_face(k[i,j+1], k[i,j]) # k_(i, j+1/2)
            I = ofs + i - 1
            A[I  , I] = r * (k1 + k2 + k3 + k4)
            A[I-1, I] = -r * k1
            A[I+1, I] = -r * k2

        # Último renglón del bloque, considera BC en la pared der.
        k1 = u_face(k[Nx-1,j  ], k[Nx,j]) # k_(i-1/2, j)
        k2 = u_face(k[Nx+1,j  ], k[Nx,j]) # k_(i+1/2, j)
        k3 = u_face(k[Nx  ,j-1], k[Nx,j]) # k_(i, j-1/2)
        k4 = u_face(k[Nx  ,j+1], k[Nx,j]) # k_(i, j+1/2)
        I = ofs + Nx - 1
        A[I-1,I] = -r * k1 
        A[I  ,I] = r * (k1 + k2 + k3 + k4) 

       
# Despues llena las dos diagonales externas
    I = 0
    for j in range(1, Ny):
        for i in range(1,Nx+1):
            k3 = u_face(k[i,j-1+1], k[i,j+1]) # k_(i, j-1/2)
            k4 = u_face(k[i,j+1], k[i,j]) # k_(i, j+1/2)
            A[I   , I+Nx] = -r * k3 # South, 3, down
            A[I+Nx, I   ] = -r * k4 # North, 4, up
            I += 1
    
    return A


def Laplaciano2D_T(Nx, Ny, r, k):
    """ Esta funcion calcula los coeficientes del 
    sistema lineal producido por el operador de 
    Laplace en 2D. Estos coeficientes son almacenados 
    en la matriz pentadiagonal correspondiente."""
    N = Nx * Ny
    A = np.zeros((N,N))

# Primero llena los bloques tridiagonales
    for j in range(1,Ny+1):
        ofs = Nx * (j-1)     
        # Primer renglón del bloque, considera BC en la pared izq.
        k1 = u_face(k[0,j  ], k[1,j]) # k_(i-1/2, j)
        k2 = u_face(k[2,j  ], k[1,j]) # k_(i+1/2, j)
        k3 = u_face(k[1,j+1], k[1,j]) # k_(i, j-1/2)
        k4 = u_face(k[1,j-1], k[1,j]) # k_(i, j+1/2)
        A[ofs    , ofs] = 1 + r * (k1 + k2 + k3 + k4) 
        A[ofs + 1, ofs] = -r * k2

        # Renglones intermedios del bloque 
        for i in range(2,Nx):
            k1 = u_face(k[i-1,j], k[i,j]) # k_(i-1/2, j)
            k2 = u_face(k[i+1,j], k[i,j]) # k_(i+1/2, j)
            k3 = u_face(k[i,j-1], k[i,j]) # k_(i, j-1/2)
            k4 = u_face(k[i,j+1], k[i,j]) # k_(i, j+1/2)
            I = ofs + i - 1
            A[I  , I] = 1 + r * (k1 + k2 + k3 + k4)
            A[I-1, I] = -r * k1
            A[I+1, I] = -r * k2

        # Último renglón del bloque, considera BC en la pared der.
        k1 = u_face(k[Nx-1,j  ], k[Nx,j]) # k_(i-1/2, j)
        k2 = u_face(k[Nx+1,j  ], k[Nx,j]) # k_(i+1/2, j)
        k3 = u_face(k[Nx  ,j-1], k[Nx,j]) # k_(i, j-1/2)
        k4 = u_face(k[Nx  ,j+1], k[Nx,j]) # k_(i, j+1/2)
        I = ofs + Nx - 1
        A[I-1,I] = -r * k1 
        A[I  ,I] = 1 + r * (k1 + k2 + k3 + k4) 

       
# Despues llena las dos diagonales externas
    I = 0
    for j in range(1, Ny):
        for i in range(1,Nx+1):
            k3 = u_face(k[i,j-1+1], k[i,j+1]) # k_(i, j-1/2)
            k4 = u_face(k[i,j+1], k[i,j]) # k_(i, j+1/2)
            A[I   , I+Nx] = -r * k3 # South, 3, down
            A[I+Nx, I   ] = -r * k4 # North, 4, up
            I += 1
    
    return A


def LeeDatos(filename):
    """ Esta funcion lee los datos de un archivo. a? y b? son 
    las coordenadas inicial y final del dominio respectivamente
    las direcciones x y y;
    Nx y Ny son el numero de incognitas en las direcciones 
    correspondientes; A, B, C y D son las condiciones de
    frontera. Se regresa la tupla (ax,bx,ay,by,Nx,Ny,A,B,C,D)."""
    ifile = open(filename, 'r')     # abre el archivo de entrada
    file_lines = ifile.readlines()  # lee las lineas del archivo
    ifile.close();                  # cierra el archivo de entrada
    ax, bx, ay, by, Nx, Ny, A, B, C, D = file_lines[0].split() # separa las columnas de la primera linea
    ax = float(ax); bx = float(bx); ay = float(ay); by = float(by); 
    Nx = int(Nx); Ny = int(Ny); 
    A = float(A); B = float(B); C = float(C); D = float(D); 
    return ax, bx, ay, by, Nx, Ny, A, B, C, D

def ImprimeDatos(ax,bx,ay,by,Nx,Ny,hx,hy,A,cd1,B,cd2,C,cd3,D,cd4):
    """ Esta funcion imprime los datos del problema a resolver."""
    print()
    print("+----------------------------------------------------+")
    print("|      Solucion de la ecuacion de Laplace en 2D      |")
    print("+----------------------------------------------------+")
    print("| Autor: Luis M. de la Cruz S.                       |")
    print("+----------------------------------------------------+")
    print("|                 Datos de entrada                   |")
    print("+----------------------------------------------------+")
    print("| Punto inicial del dominio en x   : ax = %g" % ax)
    print("| Punto final del dominio en x     : bx = %g" % bx)
    print("| Punto inicial del dominio en y   : ax = %g" % ay)
    print("| Punto final del dominio en y     : bx = %g" % by)
    print("| Numero total de incognitas en x  : Nx = %d" % Nx)
    print("| Numero total de incognitas en y  : Ny = %d" % Ny)
    print("| Numero total de incognitas       : N = %d" % (Nx*Ny))
    print("| El tamanio de la malla en x es   : hx = %g " % hx)
    print("| El tamanio de la malla en y es   : hy = %g " % hy)
    print("| Cond. de front. ", cd1, "en ax : A = %g" % A)
    print("| Cond. de front. ", cd2, "en bx : B = %g" % B)
    print("| Cond. de front. ", cd3, "en ay : C = %g" % C)
    print("| Cond. de front. ", cd4, "en by : D = %g" % D)
    print("+----------------------------------------------------+")

def ImprimeSistema(A,u,f):
    """ Esta funcion imprime el sistema lineal asi como la solucion
    del mismo, siempre y cuando su longitud sea menor o igual a 10."""
    print("\n Lado derecho del sistema : size = %d \n" % f.size, f)
    print("\n Matriz del sistema : \n", A)
    print("\n Solucion del sistema : size = %d \n" % u.size, u)

def GraficaSuperficieC(xg,yg,u,colormap):
    pl.contourf(xg, yg, u, 100, alpha=.95, cmap=colormap)
    C = pl.contour(xg, yg, u, 100, colors='black', alpha=0.01, linewidth=.5)
    pl.clabel(C, inline=1, fontsize=10)
    
    fig = pl.figure()
    ax = Axes3D(fig)
    ax.plot_surface(xg, yg, u, rstride=2, cstride=2, alpha=.95, cmap=colormap)

    pl.show()

def GuardaSolucion(filename, x, y, u):
    """ Esta funcion guarda la solucion en un archivo para su
    posterior analisis, en un archivo de nombre filename."""    
    ofile = open(filename, 'w')
    for i in range(0,x.size):
        for j in range(0,y.size):
            ofile.write('%12.10g \t %12.10g \t %12.10g\n' % (x[i], y[j],u[j,i]))
    ofile.close()



#if __name__ == "__main__":





