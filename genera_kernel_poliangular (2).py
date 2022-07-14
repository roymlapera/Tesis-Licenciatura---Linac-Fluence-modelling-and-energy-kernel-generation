import numpy as np
import matplotlib.pyplot as plt
import matplotlib.cm as cm
import matplotlib.colors as colors
import scipy.ndimage
import math
from glob           import glob 
from os.path        import splitext, split, dirname

def plotPolarContour(kernel2draw, rAxis2draw, thetaAxis2draw, Rmax, drawType):

	if drawType == 'full':
		# flippedkernel2draw = np.fliplr(kernel2draw)
		# kernel2draw = np.concatenate((kernel2draw, flippedkernel2draw) ,axis=1)
		thetaAxis2draw = np.concatenate((thetaAxis2draw ,180.0 + thetaAxis2draw),axis=0)

	deltaR = np.zeros(len(rAxis2draw))
	deltaR[0] = rAxis2draw[0]
	k=0
	while k < len(rAxis2draw)-1:
		deltaR[k+1] = rAxis2draw[k+1]-rAxis2draw[k]
		k+=1

	rAxis2draw=rAxis2draw-deltaR/2
	deltaTheta = np.zeros(len(thetaAxis2draw))
	deltaTheta[0] = thetaAxis2draw[0] 
	k=0

	while k < len(thetaAxis2draw)-1:
		deltaTheta[k+1] = thetaAxis2draw[k+1]-thetaAxis2draw[k]
		k+=1

	thetaAxis2draw=thetaAxis2draw-deltaTheta/2
	thetaAxis2draw = np.radians(thetaAxis2draw)

	thetaAxisMesh, rAxisMesh = np.meshgrid(thetaAxis2draw, rAxis2draw)

	# kernel2draw = scipy.ndimage.zoom(kernel2draw, 3)
	# rAxisMesh = scipy.ndimage.zoom(rAxisMesh, 3)
	# thetaAxisMesh = scipy.ndimage.zoom(thetaAxisMesh, 3)

	# fig = plt.figure()
	# ax = fig.add_subplot(111, polar = True)

	# polarImg = ax.pcolor(thetaAxisMesh, rAxisMesh, kernel2draw, cmap=cm.jet, norm=colors.LogNorm()) #
	# polarImg = ax.pcolormesh(thetaAxisMesh, rAxisMesh, kernel2draw, cmap=cm.jet, norm=colors.LogNorm()) #
	# polarImg = ax.pcolormesh(thetaAxis2draw, rAxis2draw, kernel2draw, cmap=cm.jet, norm=colors.LogNorm())
	# kernel2draw.max()*0.01
	# L = [kernel2draw.max()*0.001, kernel2draw.max()*0.01]

	L_min = -12
	L_max = 7
	L = np.logspace(L_min, L_max, 1000)

	l_exp = np.arange(L_min, L_max+1)

	l = 10.0**l_exp

	# print l
	# polarImg = ax.contour(thetaAxisMesh, rAxisMesh, kernel2draw,  cmap=cm.jet, levels = L)
	#plt.title('kernel')

	fig, ax = plt.subplots(subplot_kw=dict(projection='polar'))
	polarImg = ax.contourf(thetaAxisMesh, rAxisMesh, kernel2draw,L, cmap=cm.jet, norm=colors.LogNorm())

	ax.set_theta_zero_location("S")
	#ax.set_theta_direction(-1)
	ax.set_rmax(Rmax)
	ax.grid(True)
	polarImg.set_clim(np.amin(L), np.amax(L))
	cb = plt.colorbar(polarImg) 
	cb.set_ticks(l)
	path = 'C:\Users\Install\Desktop\Tesis Roy\Imagenes finales/'
	fig.savefig(path + 'kernel_rotado_45grados.png', dpi = 200, bbox_inches='tight')
	plt.show()

def plotPolarContour_full(kernel2draw, rAxis2draw, thetaAxis2draw, Rmax):

	
	flippedkernel2draw = np.fliplr(kernel2draw)
	kernel2draw = np.concatenate((kernel2draw, flippedkernel2draw) ,axis=1)
	thetaAxis2draw = np.concatenate((thetaAxis2draw ,180.0 + thetaAxis2draw),axis=0)
	# print theta_axis
	thetaAxis2draw_ext = np.append(-thetaAxis2draw[0], thetaAxis2draw)  

	# print theta_axis_ext

	kernel2draw_ext = np.zeros([kernel2draw.shape[0],kernel2draw.shape[1]+1])
	kernel2draw_ext[:,0] = kernel2draw[:,-1]
	kernel2draw_ext[:,1:] = kernel2draw

	deltaR = np.zeros(len(rAxis2draw))
	deltaR[0] = rAxis2draw[0]
	k=0
	while k < len(rAxis2draw)-1:
		deltaR[k+1] = rAxis2draw[k+1]-rAxis2draw[k]
		k+=1

	rAxis2draw=rAxis2draw-deltaR/2
	deltaTheta = np.zeros(len(thetaAxis2draw_ext))
	deltaTheta[0] = thetaAxis2draw_ext[0] 
	k=0

	while k < len(thetaAxis2draw)-1:
		deltaTheta[k+1] = thetaAxis2draw_ext[k+1]-thetaAxis2draw_ext[k]
		k+=1

	thetaAxis2draw_ext=thetaAxis2draw_ext-deltaTheta/2
	thetaAxis2draw_ext = np.radians(thetaAxis2draw_ext)

	thetaAxisMesh, rAxisMesh = np.meshgrid(thetaAxis2draw_ext, rAxis2draw)

	# kernel2draw = scipy.ndimage.zoom(kernel2draw, 3)
	# rAxisMesh = scipy.ndimage.zoom(rAxisMesh, 3)
	# thetaAxisMesh = scipy.ndimage.zoom(thetaAxisMesh, 3)

	# fig = plt.figure()
	# ax = fig.add_subplot(111, polar = True)

	# polarImg = ax.pcolor(thetaAxisMesh, rAxisMesh, kernel2draw, cmap=cm.jet, norm=colors.LogNorm()) #
	# polarImg = ax.pcolormesh(thetaAxisMesh, rAxisMesh, kernel2draw, cmap=cm.jet, norm=colors.LogNorm()) #
	# polarImg = ax.pcolormesh(thetaAxis2draw, rAxis2draw, kernel2draw, cmap=cm.jet, norm=colors.LogNorm())
	# kernel2draw.max()*0.01
	# L = [kernel2draw.max()*0.001, kernel2draw.max()*0.01]

	L_min = -7
	L_max = 7
	L = np.logspace(L_min, L_max, 1000)

	l_exp = np.arange(L_min, L_max+1)

	l = 10.0**l_exp

	# print l
	# polarImg = ax.contour(thetaAxisMesh, rAxisMesh, kernel2draw,  cmap=cm.jet, levels = L)
	#plt.title('kernel')

	fig, ax = plt.subplots(subplot_kw=dict(projection='polar'))
	polarImg = ax.contourf(thetaAxisMesh, rAxisMesh, kernel2draw_ext,L, cmap=cm.jet, norm=colors.LogNorm())

	ax.set_theta_zero_location("S")
	#ax.set_theta_direction(-1)
	ax.set_rmax(Rmax)
	ax.grid(True)
	polarImg.set_clim(np.amin(L), np.amax(L))
	cb = plt.colorbar(polarImg) 
	cb.set_ticks(l)
	# path = 'C:\Users\Install\Desktop\Tesis Roy\Imagenes finales/'
	# fig.savefig(path + 'kernel_poliangular.png', dpi = 200, bbox_inches='tight')
	plt.show()

def plotPolarContour_full_360(kernel2draw, rAxis2draw, thetaAxis2draw, Rmax):

	
	# flippedkernel2draw = np.fliplr(kernel2draw)
	# kernel2draw = np.concatenate((kernel2draw, flippedkernel2draw) ,axis=1)
	thetaAxis2draw = np.concatenate((thetaAxis2draw ,180.0 + thetaAxis2draw),axis=0)
	# print theta_axis
	thetaAxis2draw_ext = np.append(-thetaAxis2draw[0], thetaAxis2draw)  

	# print theta_axis_ext

	kernel2draw_ext = np.zeros([kernel2draw.shape[0],kernel2draw.shape[1]+1])
	kernel2draw_ext[:,0] = kernel2draw[:,-1]
	kernel2draw_ext[:,1:] = kernel2draw

	deltaR = np.zeros(len(rAxis2draw))
	deltaR[0] = rAxis2draw[0]
	k=0
	while k < len(rAxis2draw)-1:
		deltaR[k+1] = rAxis2draw[k+1]-rAxis2draw[k]
		k+=1

	rAxis2draw=rAxis2draw-deltaR/2
	deltaTheta = np.zeros(len(thetaAxis2draw_ext))
	deltaTheta[0] = thetaAxis2draw_ext[0] 
	k=0

	while k < len(thetaAxis2draw)-1:
		deltaTheta[k+1] = thetaAxis2draw_ext[k+1]-thetaAxis2draw_ext[k]
		k+=1

	thetaAxis2draw_ext=thetaAxis2draw_ext-deltaTheta/2
	thetaAxis2draw_ext = np.radians(thetaAxis2draw_ext)

	thetaAxisMesh, rAxisMesh = np.meshgrid(thetaAxis2draw_ext, rAxis2draw)

	# kernel2draw = scipy.ndimage.zoom(kernel2draw, 3)
	# rAxisMesh = scipy.ndimage.zoom(rAxisMesh, 3)
	# thetaAxisMesh = scipy.ndimage.zoom(thetaAxisMesh, 3)

	# fig = plt.figure()
	# ax = fig.add_subplot(111, polar = True)

	# polarImg = ax.pcolor(thetaAxisMesh, rAxisMesh, kernel2draw, cmap=cm.jet, norm=colors.LogNorm()) #
	# polarImg = ax.pcolormesh(thetaAxisMesh, rAxisMesh, kernel2draw, cmap=cm.jet, norm=colors.LogNorm()) #
	# polarImg = ax.pcolormesh(thetaAxis2draw, rAxis2draw, kernel2draw, cmap=cm.jet, norm=colors.LogNorm())
	# kernel2draw.max()*0.01
	# L = [kernel2draw.max()*0.001, kernel2draw.max()*0.01]

	L_min = -7
	L_max = 7
	L = np.logspace(L_min, L_max, 1000)

	l_exp = np.arange(L_min, L_max+1)

	l = 10.0**l_exp

	# print l
	# polarImg = ax.contour(thetaAxisMesh, rAxisMesh, kernel2draw,  cmap=cm.jet, levels = L)
	#plt.title('kernel')

	fig, ax = plt.subplots(subplot_kw=dict(projection='polar'))
	polarImg = ax.contourf(thetaAxisMesh, rAxisMesh, kernel2draw_ext,L, cmap=cm.jet, norm=colors.LogNorm())

	ax.set_theta_zero_location("S")
	#ax.set_theta_direction(-1)
	ax.set_rmax(Rmax)
	ax.grid(True)
	polarImg.set_clim(np.amin(L), np.amax(L))
	cb = plt.colorbar(polarImg) 
	cb.set_ticks(l)
	# path = 'C:\Users\Install\Desktop\Tesis Roy\Imagenes finales/'
	# fig.savefig(path + 'kernel_rotado_90grados.png', dpi = 200, bbox_inches='tight')
	plt.show()


def GeneraKernel360(kernel180):                                                                   # A partir de un kernel de 180 grados genera uno de 360 (espejando)
	flippedkernel180 = np.fliplr(kernel180)
	kernel360 = np.concatenate((kernel180, flippedkernel180) ,axis=1)
	return kernel360

def RotaKernel(kernel, nc_rotacion, NR, NC):                                                      # Rota el kernel "nc_rotacion" grados 

	kernel_aux = np.copy(kernel)

	for i in range(nc_rotacion):
		for fila in range(NR):
			kernel_aux[fila,:] = np.concatenate((kernel_aux[fila,-1:],kernel_aux[fila,:-1]))

	return kernel_aux  

#------------------------------------------------------------------------------------------------------------------------------------------------------------

def findFiles(phsp_path):                                                                           # FINDS FILES INSIDE FOLDER WITH IAEAphsp EXTENSION
    find_names  = phsp_path + '/*.IAEAphsp'                                                       # Names are assumed to be .......... .IAEAheader
    all_names   = glob(find_names)                                                                  # Found files with the searched names
    
    phsp_files  = []                                                                                # List to put correct names
    for name in all_names: phsp_files.append(splitext(name)[0])                                     # Appends name without extension to files
    
    print 'Espacios de fases encontrados:'

    print '\n'.join(phsp_files)

    phsp_files = np.delete(phsp_files,[0,1,3,4,5,6])

    print 'Abriendo:'

    print '\n'.join(phsp_files)

    return phsp_files

def parseFile(phsp_file):                                                                           # PARSES IAEA FILE FOR PARTICLES AND PACK FORMAT
    phsp_header = phsp_file + '.IAEAheader'                                                         # Header file with full information

    if (phsp_file == 'C:\Users\Install\Documents\EGSnrc\Phase Space\VarianClinaciX_6MV_20x20_w1'):
        with open(phsp_header, 'r') as fop: lines = fop.readlines()                                     # Opens file to operate. IT will closed automatically
        for index, line in enumerate(lines):                                                            # Search in every line for keywords

            if ('$BYTE_ORDER' in line) and ('1234' in lines[index + 1]):                                # Determins the type of unpack format
                packFormat  = np.dtype([('Type', '<i1'), ('Data', '<f4', (7)), ('Extra', '<i4', (2))])   # Assuming IAEA format with ZLAST
            if ('$BYTE_ORDER' in line) and ('4321' in lines[index + 1]):
                packFormat  = np.dtype([('Type', '>i1'), ('Data', '>f4', (7)), ('Extra', '>i4', (2))])   # Assuming IAEA format with ZLAST
            if '$PARTICLES' in line:  particles = int(lines[index + 1])                                 # Total number of particles in file
    else:
        with open(phsp_header, 'r') as fop: lines = fop.readlines()                                     # Opens file to operate. IT will closed automatically
        for index, line in enumerate(lines):                                                            # Search in every line for keywords

            if ('$BYTE_ORDER' in line) and ('1234' in lines[index + 1]):                                # Determins the type of unpack format
                packFormat  = np.dtype([('Type', '<i1'), ('Data', '<f4', (7)), ('Extra', '<i4', (1))])   # Assuming IAEA format with ZLAST
            if ('$BYTE_ORDER' in line) and ('4321' in lines[index + 1]):
                packFormat  = np.dtype([('Type', '>i1'), ('Data', '>f4', (7)), ('Extra', '>i4', (1))])   # Assuming IAEA format with ZLAST
            if '$PARTICLES' in line:  particles = int(lines[index + 1])                                 # Total number of particles in file
    
    return particles, packFormat

#------------------------------------------------------------------------------------------------------------------------------------------------------------

def DevuelveAngulos(data, type, maxEnergy):                                                          # Devuelve los angulos directores respecto a z y los pesos estadisticos de los e-  
    energy      = np.abs(data[:, 0])
    weight      = data[:, 5]
    zCos          = (1 - (data[:, 3])**2 - (data[:, 4])**2)**.5

    isElectron    = 1*(type == 2) * (energy < maxEnergy) 

    ang = np.arccos(zCos)
    weight = weight*isElectron*energy
    ang = ang[weight>0]
    weight = weight[weight>0]

    return ang, weight

def DevuelveAngulosDeArchivo(data_path, maxEnergy):                                                  # HAce lo mismo que DevuelveAngulos pero para cada pack y junta todos los datos
    phsp_files   = findFiles(data_path)                                        
    for phsp in phsp_files:
        particles, packFormat = parseFile(phsp)                                                 
        
        with open(phsp + '.IAEAphsp', 'rb') as IAEAphsp:                                        
            pack       = np.fromfile(IAEAphsp, dtype=packFormat, count=-1)                  
            data, weight   = DevuelveAngulos(pack['Data'], pack['Type'], maxEnergy)           
    return data, weight  	

def GeneraHistograma(data, weight, bins):
	histo, bin_edges = np.histogram(data, bins=bins, range=[0,np.pi/2.0], weights=weight, normed=True)
	return histo

def Genera_PesoAngular(path, filename, bins, maxEnergy):                       # Toma el path y el nombre de archivo, lo abre, grafica y devuelve los pesos estadisticos ya normalizados
	peso_angular = AbreHistograma(path + filename)

	fig = plt.figure() 

	plt.plot(peso_angular, "r.")

	plt.show()

	print peso_angular.sum() * np.pi/180.0

	return peso_angular

#-----------------------------------------------------------------------------------------------------------------------------------------------------------

def GuardaHistograma(archivo,data): np.savetxt(archivo, data)

def AbreHistograma(file):
    with open(file) as f:
        histo = np.loadtxt(f, delimiter="\n", dtype='double', comments="#", skiprows=0, usecols=None)
    return histo

#-------------------------------------------------------------------------------- MAIN ------------------------------------------------------------------------------

#Rotacion del kernel

NC = 180
NR = 210

nc_rotacion = 1 #Guardo los kernels rotados cada 2 grados

# archivo = np.load("kernel_polienergetico_pablo.npz")

# kernel_poli = archivo['kernel']
# rAxis     = archivo['rAxis']
# thetaAxis = archivo['thetaAxis']

# plotPolarContour_full(kernel_poli, rAxis, thetaAxis, 4)

# kernel_poli = GeneraKernel360(kernel_poli)

# kernels = []

# kernels.append(kernel_poli)

# for i in range(90/nc_rotacion-1):					# Genero los kernels rotados a partir de rotar el kernel_poli y los meto todos en un vector de kernels llamado "kernels"
# 	kernel_rotado = RotaKernel(kernel_poli,nc_rotacion*i,NR,NC)
# 	kernels.append(kernel_rotado)

#Grafico un kernel para mostrar como se ve no mas

# kernel_rotado = RotaKernel(kernel_poli,nc_rotacion*45,NR,NC)

# plotPolarContour_full_360(kernel_rotado, rAxis, thetaAxis, 4)

# np.savez("kernel_monoangular_45.npz", kernel = kernel_rotado, rAxis = rAxis, thetaAxis = thetaAxis)



archivo = np.load("kernel_poliangular_pablo_con_pesos.npz")

kernel_poli = archivo['kernel']
rAxis     = archivo['rAxis']
thetaAxis = archivo['thetaAxis']

plotPolarContour_full(kernel_poli, rAxis, thetaAxis, 4)

#-------------------------------------------------------------------------------------------------------------------------------------------------------------------

#Genero los pesos angulares

path = 'C:\Users\Install\Documents\EGSnrc\Phase Space/'

filename =  'Hist_angulos.txt'

filename_pablo  =  'espectro_incidencia.txt'

bins = 90

maxEnergy = 6.0

# Angulos, weight = DevuelveAngulosDeArchivo(path, maxEnergy)

# Histograma = GeneraHistograma(Angulos, weight, bins)

# fig, ax = plt.subplots()

# plt.hist(Angulos, bins=bins, histtype='step', normed=True, color='b', weights=weight)

# plt.axis([-1.1, 1.1, 0, 1.7])

# plt.xlabel('Angulo director respecto al eje Z', fontsize=14)

# plt.grid(True)

# fig.savefig('C:\Users\Install\Desktop\Tesis Roy\Imagenes finales/pesos_ang_dir_en_z.png', dpi = 200, bbox_inches='tight')

# plt.show()

# GuardaHistograma(path + filename, Histograma)				

#Generacion del peso espectral

# peso_angular = Genera_PesoAngular(path, filename_pablo, bins, maxEnergy)

# kernel_poliangular = np.zeros((NR, NC*2))

# for i in range(len(peso_angular)):							# hago la suma pesada de kernels polienergeticos en diferentes incidencias
# 	kernel_poliangular += peso_angular[i] * kernels[i] * (np.pi/2.0)/90.0


# kernel_poliangular = kernel_poliangular[:,:-180]

# # print kernel_poliangular.shape

# plotPolarContour_full(kernel_poliangular, rAxis, thetaAxis, 4)

# np.savez("poliangular_pablo.npz", kernel = kernel_poliangular, rAxis = rAxis, thetaAxis = thetaAxis)