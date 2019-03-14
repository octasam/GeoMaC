import sys

def read():
	try:
		in_file_name = sys.argv[1]
		out_file_name = sys.argv[2]
	except:
		mensaje = """ Error: La ejecucion de este programa es como sigue: 

		python programa.py entrada.txt salida.txt

		programa.py: es el programa en Python a ejecutar.
		entrada.txt: es el nombre de un archivo que contiene los datos del problema. 
		salida.txt: se usa para almacenar la solucion del problema."""
		print(mensaje)
		sys.exit(1)

	ifile = open(in_file_name, 'r') # abre el archivo de entrada
	ifile_lines = ifile.readlines() # lee las lineas del archivo
	ifile.close()                   # cierra el archivo de entrada
	a, b, c, d, Nx, Ny, A, B, C, D = ifile_lines[0].split() # separa las columnas de la primera linea
	a = float(a)
	b = float(b)
	c = float(c)
	d = float(d)
	Nx = int(Nx)
	Ny = int(Ny)
	A = float(A)
	B = float(B)
	C = float(C)
	D = float(D)

	return (a,b,c,d,Nx,Ny,A,B,C,D,out_file_name)

