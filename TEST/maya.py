import numpy as np

def mesh(a,b,c,d,Nx,Ny):
	x = np.linspace(a,b,Nx)
	y = np.linspace(c,d,Ny)
	return x, y

def deltas(a,b,c,d,Nx,Ny):
	hx =  (b-a)/(Nx+1)
	hy = (d-c)/(Ny+1)
	return hx, hy

if __name__ == "__main__":

	a = 0
	b = 1
	c = 0
	d = 1
	Nx = 9
	Ny = 9

	x, y = mesh(a,b,c,d,Nx,Ny)
	hx, hy = deltas(a,b,c,d,Nx,Ny)

	print(x, y)
	print(hx, hy)