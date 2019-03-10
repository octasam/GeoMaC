import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
import pylab as pl

def colormap(u,colormap,name = 'False'):
	plt.imshow(u,colormap)
	if name != 'False':
		plt.savefig(name)
	plt.show()

def surface(x,y,u,colormap,name = 'False'):
	xg, yg = np.meshgrid(x,y)
	fig = pl.figure()
	ax = Axes3D(fig)
	ax.plot_surface(xg, yg, u, rstride=2, cstride=2, alpha=.95, cmap=colormap)
	if name != 'False':
		plt.savefig(name)
	plt.show()

def contour(x,y,u,colormap,name = 'False'):
	xg, yg = np.meshgrid(x,y)
	pl.contourf(xg, yg, u, 100, alpha=.95, cmap=colormap)
	C = pl.contour(xg, yg, u, 100, colors='black', alpha=0.01, linewidth=.5)
	pl.clabel(C, inline=1, fontsize=10)

	if name != 'False':
		pl.savefig(name)

	pl.show()

def GuardaSolucion(filename, x, y, u):
    """ Esta funcion guarda la solucion en un archivo para su
    posterior analisis, en un archivo de nombre filename."""    
    ofile = open(filename, 'w')
    for i in range(0,x.size):
        for j in range(0,y.size):
            ofile.write('%12.10g \t %12.10g \t %12.10g\n' % (x[i], y[j],u[j,i]))
    ofile.close()


def write(x,y,u,out_file_name):
	ofile = open(out_file_name,'w') # abre archivo para salida
	for i in range(0,x.size):
		for j in range(0,y.size):
			ofile.write('%12.10g \t %12.10g \t %12.10g\n' % (x[i], y[j],u[j,i]))
	ofile.close()
