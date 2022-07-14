import numpy as np
import matplotlib.pyplot as plt
import matplotlib.cm as cm
import matplotlib.colors as colors
import scipy.ndimage
import math
from glob           import glob 
from os.path        import splitext, split, dirname

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



def calculatePackFluenciaEnergetica(data, type, maxEnergy, bins):                                                     # CALCULATES FULL SPECTRUM FROM PACK DATA
    spectrum    = np.zeros(bins)                                                                       # Energy spectrum (E)
        
    energy      = np.abs(data[:, 0])
    weight      = data[:, 5]                                                                        # Weight                     of every photon in Phase plane of current unpack
    
    isElectron    = 1*(type == 2) * (energy < maxEnergy)					 							# Only photons inside calculation limit are proccessed
    eneInd      = isElectron * (energy * bins / maxEnergy).astype(int)                                # Index for spectrum sumation
    sumation    = isElectron * weight * energy

    for ind, value in enumerate(eneInd): 
    	spectrum[value] += sumation[ind]                           # Spectrum photon countingt 
                                                                  
    return spectrum

def CalculaEspectroFluenciaEnergetica(data_path, maxEnergy, bins):                                    # OPENS IAEA FILE AND CALCULATES PARTICLE SPECTRUM  IN small DISC
    spectrum     = np.zeros(bins)                                                                      # Energy spectrum (E)
    phsp_files   = findFiles(data_path)                                                  # Finds IAEA files to work with

    for phsp in phsp_files:
        particles, packFormat = parseFile(phsp)                                                     # Specifies format and amount of particles to read
        
        with open(phsp + '.IAEAphsp', 'rb') as IAEAphsp:                                            # Opens binary file to read data
            for x in range(particles/100000):
                pack        = np.fromfile(IAEAphsp, dtype=packFormat, count=100000)                   # Gets full pack data from one file
                spectrum   += calculatePackFluenciaEnergetica(pack['Data'], pack['Type'], maxEnergy, bins)            # Calculates spectrum part from one file

    spectrum /= (spectrum.sum() * (1.0*maxEnergy/bins))

    return spectrum

def GuardaEspectroFluenciaEnergetica(archivo,data): np.savetxt(archivo, data)

def AbreEspectroTXT(file):
    with open(file) as f:
        spectrum = np.loadtxt(f, delimiter="\n", dtype='double', comments="#", skiprows=0, usecols=None)
    return spectrum

def Genera_PesoEspectral(path, filename, bins, maxEnergy):
	peso_espectral = AbreEspectroTXT(path + filename)

	fig = plt.figure() 

	plt.plot(peso_espectral, "r.")

	plt.show()

	print peso_espectral.sum()*6.0/8

	return peso_espectral

def GraficaHistogramaEnergiaCinetica(data, maxEnergy, bins, path):									
    
    masa_en_reposo_electron = 0.511  #MeV

    #normalizacion
    # data /= ( data.sum() * (maxEnergy) )

    abcisas = np.arange(bins)

    abcisas = abcisas /float(bins-1) * maxEnergy

    fig = plt.figure() 

    abcisas = np.ravel(zip(abcisas,abcisas+maxEnergy/(bins-1)))
    data= np.ravel(zip(data,data))

    plt.plot(abcisas, data, "r")

    plt.axis([0, 6.1, 0.0, 0.38])

    plt.xlabel('Energia [MeV]', fontsize=14)
    plt.ylabel('Fluencia energetica [MeV]', fontsize=14)

    plt.grid(True)

    plt.show()

    fig.savefig(path + 'espectro_fluencia_energetica.png', dpi = 200, bbox_inches='tight')


def plotPolarContour(kernel2draw, rAxis2draw, thetaAxis2draw, Rmax, drawType):

	if drawType == 'full':
		flippedkernel2draw = np.fliplr(kernel2draw)
		kernel2draw = np.concatenate((kernel2draw, flippedkernel2draw) ,axis=1)
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

	L_min = -7
	L_max = 7
	L = np.logspace(L_min, L_max, 1000)

	l_exp = np.arange(L_min, L_max+1)

	l = 10.0**l_exp

	# print l
	# polarImg = ax.contour(thetaAxisMesh, rAxisMesh, kernel2draw,  cmap=cm.jet, levels = L)
	# plt.title('kernel')

	fig, ax = plt.subplots(subplot_kw=dict(projection='polar'))
	polarImg = ax.contourf(thetaAxisMesh, rAxisMesh, kernel2draw,L, cmap=cm.jet, norm=colors.LogNorm())

	ax.set_theta_zero_location("S")
	#ax.set_theta_direction(-1)
	ax.set_rmax(Rmax)
	ax.grid(True)
	polarImg.set_clim(np.amin(L), np.amax(L))
	cb = plt.colorbar(polarImg) 
	cb.set_ticks(l)
	path = 'C:\EGSnrc\EGS_HOME\edknrc/'
	# fig.savefig(path + 'fullWater_50M_poli.png', dpi = 200, bbox_inches='tight')
	plt.show()

	# plt.figure()	

	# plt.plot(rAxis, kernel2draw[:,0], "b.")
	# plt.plot(rAxis, kernel2draw[:,10], "g.")
	# plt.plot(rAxis, kernel2draw[:,20], "r.")
	# plt.plot(rAxis, kernel2draw[:,30], "k.")
	# plt.plot(rAxis, kernel2draw[:,40], "y.")
	# plt.plot(rAxis, kernel2draw[:,50], "m.")
	# plt.plot(rAxis, kernel2draw[:,60], "c.")
	# plt.plot(rAxis, kernel2draw[:,70], "b.")

	# plt.axis([0, 4, 1e-6, 1e7])

	# plt.yscale('log')

	# plt.show()

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
	path = 'C:\Users\Install\Desktop\Tesis Roy\Imagenes finales/'
	fig.savefig(path + 'kernel_polienergetico.png', dpi = 200, bbox_inches='tight')
	plt.show()

#----------------------------------------------------------------------- MAIN --------------------------------------------------------------------

#Espectro de fluencia energetica

path = 'C:\Users\Install\Documents\EGSnrc\Phase Space/'

# filename =  'Espectro_fluencia_energetica_w1.txt'

filename =  'Espectro_fluencia_energetica_w1_100bins.txt'

bins = 100

maxEnergy = 6.0

# Espectro = CalculaEspectroFluenciaEnergetica(path, maxEnergy, bins)

# GuardaEspectroFluenciaEnergetica(path + filename,Espectro)

Espectro = AbreEspectroTXT(path + filename)

GraficaHistogramaEnergiaCinetica(Espectro, maxEnergy, bins, path)							

# # #Generacion del peso espectral

# peso_espectral = Genera_PesoEspectral(path, filename, bins, maxEnergy)

# # #Kernels polienergeticos

# kernel_full_0375MeV = np.load("fullWater_50M_0375MeV.npz")
# kernel_full_1125MeV = np.load("fullWater_50M_1125MeV.npz")
# kernel_full_1875MeV = np.load("fullWater_50M_1875MeV.npz")
# kernel_full_2625MeV = np.load("fullWater_50M_2625MeV.npz")
# kernel_full_3375MeV = np.load("fullWater_50M_3375MeV.npz")
# kernel_full_4125MeV = np.load("fullWater_50M_4125MeV.npz")
# kernel_full_4875MeV = np.load("fullWater_50M_4875MeV.npz")
# kernel_full_5625MeV = np.load("fullWater_50M_5625MeV.npz")

# kernel_full_0375MeV_kernel = kernel_full_0375MeV['kernel']
# kernel_full_1125MeV_kernel = kernel_full_1125MeV['kernel']
# kernel_full_1875MeV_kernel = kernel_full_1875MeV['kernel']
# kernel_full_2625MeV_kernel = kernel_full_2625MeV['kernel']
# kernel_full_3375MeV_kernel = kernel_full_3375MeV['kernel']
# kernel_full_4125MeV_kernel = kernel_full_4125MeV['kernel']
# kernel_full_4875MeV_kernel = kernel_full_4875MeV['kernel']
# kernel_full_5625MeV_kernel = kernel_full_5625MeV['kernel']


# kernel_poli = peso_espectral[0] * kernel_full_0375MeV_kernel + peso_espectral[1] * kernel_full_1125MeV_kernel + peso_espectral[1] * kernel_full_1125MeV_kernel + peso_espectral[2] * kernel_full_1875MeV_kernel +  peso_espectral[3] * kernel_full_2625MeV_kernel + peso_espectral[4] * kernel_full_3375MeV_kernel + peso_espectral[5] * kernel_full_4125MeV_kernel + peso_espectral[6] * kernel_full_4875MeV_kernel + peso_espectral[7] * kernel_full_5625MeV_kernel

# kernel_poli = kernel_poli * (maxEnergy/8)

# rAxis     = kernel_full_0375MeV['rAxis']
# thetaAxis = kernel_full_0375MeV['thetaAxis']

# np.savez("fullWater_50M_polienergetico.npz", kernel = kernel_poli, rAxis = rAxis, thetaAxis = thetaAxis)

# plotPolarContour_full(kernel_poli, rAxis, thetaAxis, 4)

# #-------------------------------------------------------------------------------- 5/12 --------------------------------------------------------------------------

# #Checkeo normalizacion del kernel

# simulatedTheta = 180 #cantidad de angulos equiespaciados en la generacion MC del kernel
# simulatedRadii = 210

# header_rows = 216
# footer_rows = 8

# monoEDK = np.genfromtxt('fullWater_50M_0375MeV.egslst', skip_header=header_rows, skip_footer=footer_rows, usecols=(0, 1, 2), unpack=True)

# rAxis = monoEDK[1,0:simulatedRadii]
# thetaAxis = np.zeros(simulatedTheta)

# deltaR = np.zeros(len(rAxis))
# deltaR[0] = rAxis[0]
# k=0
# while k < simulatedRadii-1:
#     deltaR[k+1] = rAxis[k+1]-rAxis[k]
#     k+=1
    
# # shift rAxis to the center of the voxel

# rAxis_mid = rAxis - deltaR/2

# for k in range(simulatedTheta):
# #    print k
#     start = k*simulatedRadii
#     # print start
#     stop = (k+1)*simulatedRadii
#     # print stop
#     # kernel_poli[:,k] = monoEDK[2,start:stop]
#     thetaAxis[k] = monoEDK[0,start]

# RAxis = np.concatenate((np.array([0.0]), rAxis[0:-1]) ,axis=0)
# ThetaAxis = np.concatenate((np.array([0.0]), thetaAxis[0:-1]) ,axis=0)

# volume = np.zeros(kernel_poli.shape)

# for k in range(simulatedTheta):
#    volume[:,k] = 2.0*np.pi*(np.cos(np.radians(ThetaAxis[k]))-np.cos(np.radians(thetaAxis[k])))*(RAxis*deltaR**2+deltaR*RAxis**2+(deltaR**3)/3)

# # flippedvolume = np.fliplr(volume)
# # volume = np.concatenate((volume, flippedvolume) ,axis=1)

# # flippedkernel_poli = np.fliplr(kernel_poli)
# # kernel_poli = np.concatenate((kernel_poli, flippedkernel_poli) ,axis=1)
   
# kernel_poli = kernel_poli*volume

# print kernel_poli.sum()