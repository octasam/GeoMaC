import numpy as np

def createMesh(a,b,c,d,Nx,Ny):
    """
    Esta funcion hace la malla ...
    """
    
    x = np.linspace(a,b,Nx+2)
    y = np.linspace(c,d,Ny+2)

    hx = (b-a)/(Nx+1)
    hy = (d-c)/(Ny+1)

    Lx = b - a
    Ly = d - c

    return Lx, Ly, hx, hy, x, y

#if __name__ == "__main__":