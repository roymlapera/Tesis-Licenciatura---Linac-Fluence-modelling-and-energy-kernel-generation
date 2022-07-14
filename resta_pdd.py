import numpy as np
import matplotlib.pyplot as plt
import os
from glob           import glob 
from os.path        import splitext, split, dirname
from mpl_toolkits.axes_grid1 import make_axes_locatable


def findFiles(phsp_path):                                                                           # FINDS FILES INSIDE FOLDER WITH IAEAphsp EXTENSION
    find_names  = phsp_path + '/*.IAEAphsp'                                                       # Names are assumed to be .......... .IAEAheader
    all_names   = glob(find_names)                                                                  # Found files with the searched names
    
    phsp_files  = []                                                                                # List to put correct names
    for name in all_names: phsp_files.append(splitext(name)[0])                                     # Appends name without extension to files
    
    print 'Archivos encontrados:'

    print phsp_files

    phsp_files = np.delete(phsp_files,[1,2,3])

    print 'Abriendo:'

    print phsp_files

    return phsp_files

def parseFile(phsp_file):                                                                           # PARSES IAEA FILE FOR PARTICLES AND PACK FORMAT
    phsp_header = phsp_file + '.IAEAheader'                                                         # Header file with full information

    if (phsp_file == 'C:\EGSnrc\EGS_HOME\dosxyznrc\VarianClinaciX_6MV_20x20_w1'):
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


#Calcula perfil de intensidad para X = 0 o para Y = 0

def CalculaPerfilIntensidad_Pack(data, type, maxEnergy, bins, axis):                                          # CALCULATES FULL SPECTRUM FROM PACK DATA
    perfil    = np.zeros(bins)                                                                       # Espectro de coseno director de ELECTRONES
    epsilon   = 2.0																					# Considero particulas entre - epsilon y + epsilon (esto esta en cm)

    energy      = np.abs(data[:, 0])
    PosX        = data[:, 1] 
    PosY        = data[:, 2]
    weight      = data[:, 5]                                                                       	 # Weight of every photon in Phase plane of current unpack
    
    if(axis=='x'):
		isElectron    = 1*(type == 2) * (energy < maxEnergy) * (np.abs(PosX) < epsilon) 				# Se seleccionan los e- que cumplen estar por debajo de la E nominal y con X = 0
		eneInd      = isElectron * ( (PosY + 30.0) * bins / 60.0).astype(int)                            # Vector de bins de posicion, se suma 30 cm  para que no haya numeros negativos y 
		sumation    = isElectron * weight 																# luego se normaliza al numero de bins deseado   
    elif(axis=='y'):
    	isElectron    = 1*(type == 2) * (energy < maxEnergy) * (np.abs(PosY) < epsilon) 				# Se seleccionan los e- que cumplen estar por debajo de la E nominal y con X = 0
    	eneInd      = isElectron * ( (PosX + 30.0) * bins / 60.0).astype(int)                            # Vector de bins de posicion, se suma 30 cm para que no haya numeros negativos y 
    	sumation    = isElectron * weight 
    else:
    	print 'ERROR'

    # print np.nanmin(energy)
    # print np.nanmax(energy)

    for ind, value in enumerate(eneInd): 															# Se enumeran los valores del array eneInd y se va sumando a perfil (en el valor 
    	perfil[value] += sumation[ind]                        									    # de eneInd) el aporte de sumation en el indice de la enumeracion de eneInd

    return perfil

def CalculaPerfilIntensidadTotal(data_path, maxEnergy, bins,axis):                                    				# OPENS IAEA FILE AND CALCULATES PARTICLE SPECTRUM  IN small DISC
    perfil     = np.zeros(bins)                                                                      			# Perfil de intensidad
    phsp_files   = findFiles(data_path)                                                  						# Finds IAEA files to work with
    for phsp in phsp_files:
        particles, packFormat = parseFile(phsp)                                                   				# Specifies format and amount of particles to read
        
        with open(phsp + '.IAEAphsp', 'rb') as IAEAphsp:                                            			# Opens binary file to read data
            for x in range(particles/100000):
                pack        = np.fromfile(IAEAphsp, dtype=packFormat, count=100000)                   			# COUNT = NUMBER OF ITEMS TO READ!!!
                perfil   += CalculaPerfilIntensidad_Pack(pack['Data'], pack['Type'], maxEnergy, bins,axis)           # Calcula la parte del perfil en el archivo

    return perfil  																		# Perfil tiene los valores del eje Y del histograma en funcion del numero de bin 

def GraficaPerfilIntensidad(data, bins):																# CAMBIAR ACA LA LEYENDA SEGUN X E Y
	
	abcisas = np.arange(bins)

	abcisas = abcisas * 600.0 /float(bins) - 300.0

	fig = plt.figure()

	plt.plot(abcisas,data,'b*')

	plt.xlabel('Posicion X [mm]', fontsize=14)

	plt.show()



def ComparaPDD_fotones_y_electrones (path, archivo_todas_las_part, archivo_fotones, archivo_electrones):
	data1 = np.loadtxt(archivo_todas_las_part, skiprows=1)

	#dosis con todo

	profundidad1 = data1[:,0].astype(np.float64) 
	dosis1 = data1[:,1].astype(np.float64) 

	#dosis solo fotones

	data2 = np.loadtxt(archivo_fotones, skiprows=1)

	profundidad2 = data2[:,0].astype(np.float64) 
	dosis2 = data2[:,1].astype(np.float64) 

	#dosis solo electrones

	data3 = np.loadtxt(archivo_electrones, skiprows=1)

	profundidad3 = data3[:,0].astype(np.float64) 
	dosis3 = data3[:,1].astype(np.float64) 


	fig = plt.figure()

	plt.plot(profundidad1, dosis1,'b', label = 'fotones + electrones')

	plt.plot(profundidad2, dosis2,'g', label = 'fotones')

	plt.plot(profundidad2, dosis1 - dosis2,'r', label = 'fotones+electrones - electrones')

	plt.plot(profundidad3, dosis3,'c', label = 'electrones')

	plt.xlabel('Profundidad [cm]', fontsize=14)
	plt.ylabel('Dosis [Gy]', fontsize=14)

	plt.legend(prop={'size':14}, loc=(0.3,0.38))

	plt.axis([0, 8, -0.1e-16, 1.2e-16])

	fig.savefig(path + 'PDDs.png', dpi = 200, bbox_inches='tight')

	plt.show()

def GraficaPDD_y_PDDpromedio (path):
	files = []

	PDDs = []

	fig = plt.figure()

	for i in os.listdir(path):
	    if os.path.isfile(os.path.join(path,i)) and 'perfil_z0.2_phsp1_electrones' in i:
	        files.append(i)

	        data = np.loadtxt(str(os.path.join(path,i)))
	        profundidad= data[:,0].astype(np.float64) 
	        dosis = data[:,1].astype(np.float64)
	        plt.plot(profundidad, dosis, "*")

	        PDDs.append(dosis)


	promedio = np.zeros(len(dosis))

	for x in range(len(dosis)):
		promedio[x] = np.average([PDD[x] for PDD in PDDs])

	print promedio
	print dosis

	# np.savetxt('C:\EGSnrc\EGS_HOME\dosxyznrc/perfil_promedio_z0.2_phsp1_electrones_1.txt', np.c_[profundidad,promedio])

	plt.plot(profundidad, promedio,label = 'promedio')

	plt.xlabel('Profundidad [cm]', fontsize=14)
	plt.ylabel('Dosis [Gy]', fontsize=14)

	plt.legend(prop={'size':14})

	plt.show()

def GraficaPerfiles (path):
	files = []

	PDDs = []

	fig = plt.figure()

	for file in os.listdir(path):
	    if os.path.isfile(os.path.join(path,file)) and ('perfil_z' in file) and file.endswith("phspjulieta_electrones.txt"):
	        files.append(file)
	        print file
	        data = np.loadtxt(str(os.path.join(path,file)))
	        x = data[:,0].astype(np.float64) 
	        intensidad = data[:,1].astype(np.float64)
	        plt.plot(x, intensidad, ".", label=str(i))

	plt.xlabel('Posicion x [cm]', fontsize=14)
	plt.ylabel('Intensidad', fontsize=14)

	plt.legend(prop={'size':14})

	plt.show()


def GraficaFluenciaYPerfilDosis ( archivo_fluencia, archivo_dosis):
	bins = 120
	fluencia      = np.loadtxt(archivo_fluencia)
	Posicionx_fluencia = np.arange(bins)
	Posicionx_fluencia = Posicionx_fluencia * 600.0 /float(bins) - 300.0

	data_dosis         = np.loadtxt(archivo_dosis)
	Posicionx_dosis    = data_dosis[:,0].astype(np.float64) * 10.0
	dosis              = data_dosis[:,1].astype(np.float64) 

	#Normalizacion
	fluencia /= fluencia.sum()
	dosis /= dosis.sum()

	fig = plt.figure()

	plt.plot(Posicionx_fluencia, fluencia, ".", label="Fluencia de particulas norm- Esp fase Disco ")
	plt.plot(Posicionx_dosis, dosis, ".", label="Dosis norm - Esp fase Disco ")

	plt.legend(prop={'size':14})
	plt.xlabel('Posicion x [mm]', fontsize=14)
	
	plt.show()



def GraficaPDDHomogeneoYNoHomogeneo(file_homog, file_no_homog):
	data1 = np.loadtxt(file_homog)
	profundidad1 = data1[:,0].astype(np.float64) 
	dosis1 = data1[:,1].astype(np.float64)

	data2 = np.loadtxt(file_no_homog)
	profundidad2 = data2[:,0].astype(np.float64)
	dosis2 = data2[:,1].astype(np.float64)

	plt.figure()

	plt.plot(profundidad1, dosis1, "b.", label = 'homog')
	plt.plot(profundidad2, dosis2, "g.", label = 'no homog')

	plt.xlabel('Profundidad [cm]', fontsize=14)
	plt.ylabel('Dosis [Gy]', fontsize=14)

	plt.legend(prop={'size':14})

	plt.show()

#----------------------------------------------------------------------------------- MAIN -----------------------------------------------------------------------------------

# archivo_todas_las_part   = 'dosis_phsp1_todo_1.txt'
# archivo_fotones          = 'dosis_phsp1_fotones_1.txt'
# archivo_electrones       = 'dosis_phsp1_electrones_1.txt'

path = 'C:\EGSnrc\EGS_HOME\dosxyznrc/'

# ComparaPDD_fotones_y_electrones (path, archivo_todas_las_part, archivo_fotones, archivo_electrones)

#-----------------------------------------------------------------------------------------------------------------------------------------------------------------------------

# path = 'C:\EGSnrc\EGS_HOME\dosxyznrc/'

# GraficaPDD_y_PDDpromedio (path)

#-----------------------------------------------------------------------------------------------------------------------------------------------------------------------------

# path = 'C:\EGSnrc\EGS_HOME\dosxyznrc/'

# GraficaPerfiles (path)

#-----------------------------------------------------------------------------------------------------------------------------------------------------------------------------

# path = 'C:\EGSnrc\EGS_HOME\dosxyznrc/'

# axis = 'y'
# Perfil = CalculaPerfilIntensidadTotal(path, 6.0, 120,axis)

# np.savetxt('C:\EGSnrc\EGS_HOME\dosxyznrc/fluencia_phsp1_1Melectrones.txt', Perfil)

# GraficaPerfilIntensidad(Perfil,120)

#-----------------------------------------------------------------------------------------------------------------------------------------------------------------------------

# archivo_fluencia_julieta =

# archivo_dosis_julieta = 'C:\EGSnrc\EGS_HOME\dosxyznrc/perfil_z0.1_phspjulieta_electrones.txt'

# archivo_fluencia_w1 = 'C:\EGSnrc\EGS_HOME\dosxyznrc/fluencia_phsp1_1Melectrones.txt'

# archivo_dosis_w1 = 'C:\EGSnrc\EGS_HOME\dosxyznrc/perfil_z0.1_phsp1_electrones_1centrado.txt'

# GraficaFluenciaYPerfilDosis (archivo_fluencia_w1, archivo_dosis_w1)

#-----------------------------------------------------------------------------------------------------------------------------------------------------------------------------

file_homog = 'PDD_homogeneo_.txt'

file_no_homog = 'PDD_no_homogeneo_.txt'

GraficaPDDHomogeneoYNoHomogeneo(file_homog, file_no_homog)










