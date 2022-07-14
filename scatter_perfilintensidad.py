import numpy 		as np
from glob           import glob 
from os.path        import splitext, split, dirname
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec


def findFiles(phsp_path):                                                                           # FINDS FILES INSIDE FOLDER WITH IAEAphsp EXTENSION
    find_names  = phsp_path + '/*.IAEAphsp'                                                       # Names are assumed to be .......... .IAEAheader
    all_names   = glob(find_names)                                                                  # Found files with the searched names
    
    phsp_files  = []                                                                                # List to put correct names
    for name in all_names: phsp_files.append(splitext(name)[0])                                     # Appends name without extension to files
    
    print 'Espacios de fases encontrados:'

    print '\n'.join(phsp_files)

    phsp_files = np.delete(phsp_files,[1,0,3,4,5,6])

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


#Crea un mapa de la posicion (x,y) de un tipo de particulas:

def PosXYPack(data, type, maxEnergy):                                                     #  POS X Y DE LAS PARTICULAS DE PACK DATA
    PosXparticula    = data[:, 1]                                                                  
    PosYparticula    = data[:, 2]    
    energy           = np.abs(data[:, 0])
    weight           = data[:, 5] 

    isElectron  = 1*(type == 2) * (energy < maxEnergy)
    
    PosX = PosXparticula[np.nonzero(isElectron)]											# El array PosX contiene la posicion X de las particulas que cumplen isElectron == 1
    PosY = PosYparticula[np.nonzero(isElectron)]

    weight = weight[np.nonzero(isElectron)]

    return PosX, PosY, weight

def PuntosXYparticulas(data_path, maxEnergy):                                    # ABRE EL ARCHIVO IAEA Y DEVUELVE DOS VECTORES CON POS X E Y DE LAS PARTICULAS  IN small DISC
    VectorX     = [0]
    VectorY     = [0]      
    W     = [0]                                                                 
    phsp_files   = findFiles(data_path)                                                  # Finds IAEA files to work with
    
    for phsp in phsp_files:
        particles, packFormat = parseFile(phsp)                                                     # Specifies format and amount of particles to read
    
        bigpack     = 100000
        smallpack   = particles % bigpack

        with open(phsp + '.IAEAphsp', 'rb') as IAEAphsp:                                            # Opens binary file to read data
            for x in range(particles/(bigpack)):
                pack        = np.fromfile(IAEAphsp, dtype=packFormat, count=bigpack)             # Gets full pack data from one file Cuanto mas chico el count menos electrones considera
                PosX,PosY, weight   = PosXYPack(pack['Data'], pack['Type'], maxEnergy) 
                VectorX = np.concatenate([VectorX,PosX])
                VectorY = np.concatenate([VectorY,PosY])
                W = np.concatenate([W,weight])

            pack        = np.fromfile(IAEAphsp, dtype=packFormat, count=smallpack)             # Gets full pack data from one file Cuanto mas chico el count menos electrones considera
            PosX,PosY, weight   = PosXYPack(pack['Data'], pack['Type'], maxEnergy) 
            VectorX = np.concatenate([VectorX,PosX])
            VectorY = np.concatenate([VectorY,PosY])
            W = np.concatenate([W,weight])
    return VectorX, VectorY, W

def GraficaPtosXY(VectorX,VectorY):
	plt.scatter(VectorX,VectorY)

	plt.show()

# Guardado de datos

def GuardaPerfilIntensidad(data, file):
	np.savetxt(file, data)

def GuardaPosicionXY(x, y, weight, file1, file2, file3):       #Guarda una estructura de dos vectores en dos txt diferentes
    np.savetxt(file1, x)
    np.savetxt(file2, y)
    np.savetxt(file3, weight)


def AbreVectorTXT(file):
	with open(file) as f:
  		vector = np.loadtxt(f, delimiter="\n", dtype='double', comments="#", skiprows=0, usecols=None)
  	return vector

def Grafica_PosXY_perfiles(x, y, weight, perfil_x, perfil_y, bins, path):
    abcisas = np.arange(bins)
    abcisas = abcisas * 600.0 /float(bins) - 300.0

    # perfil_x /= perfil_x.sum()
    # perfil_y /= perfil_y.sum()

    fig = plt.figure()

    gs = gridspec.GridSpec(2, 2,
                       width_ratios=[8,4],
                       height_ratios=[4,8])

    sub1 = fig.add_subplot(gs[0])
    sub1.set_xticks([])
    sub1.set_yticks(np.arange(0, max(perfil_x), 0.2))
    sub1.plot(abcisas, perfil_x, "g.")
    sub1.tick_params(labelsize=10)
    sub1.grid()

    sub2 = fig.add_subplot(gs[3])
    sub2.set_yticks([])
    sub2.set_xticks(np.arange(0, max(perfil_y), 0.2))
    sub2.plot(perfil_y, abcisas, "g.")
    sub2.tick_params(labelsize=10)
    sub2.grid()
 
    hist, xedges, yedges = np.histogram2d(x*10.0, y*10.0, bins=[121,121], range=[[-300, 300], [-300, 300]], normed=None, weights=weight)

    # np.savetxt('fluencia_energetica_imshow_ejes.txt', np.c_[xedges,yedges])
    # np.savetxt('fluencia_energetica_imshow_histo.txt', hist)

    binwidth = 2
    xymax = np.max( [np.max(np.fabs(perfil_x)), np.max(np.fabs(perfil_y))] )
    lim = ( int(xymax/binwidth) + 1) * binwidth

    sub3 = fig.add_subplot(gs[2], aspect='equal')
    im3 = sub3.imshow(hist, extent=[-300,300,-300,300], interpolation='None')
    sub3.axhline(linewidth=2, color='g')
    sub3.axvline(linewidth=2, color='g')
    sub3.axis([-300, 300, -300, 300])
    sub3.tick_params(labelsize=10)
    sub3.set_aspect('equal')
    # cm3 = plt.colorbar(im3)

    plt.tight_layout()
    plt.show()
    fig.savefig(path + 'fluencia energetica.png', dpi = 200, bbox_inches='tight')

#Fluencia energetica

def CalculaFluenciaEnergetica_Pack(data, type, maxEnergy, bins, axis):                                          # CALCULATES FULL SPECTRUM FROM PACK DATA
    perfil    = np.zeros(bins)                                                                       # Espectro de coseno director de ELECTRONES
    epsilon   = 2.0                                                                                 # Considero particulas entre - epsilon y + epsilon (esto esta en cm)

    energy      = np.abs(data[:, 0])
    PosX        = data[:, 1] 
    PosY        = data[:, 2]
    weight      = data[:, 5]                                                                         # Weight of every photon in Phase plane of current unpack
    
    if(axis=='x'):
        isElectron    = 1*(type == 2) * (energy < maxEnergy) * (np.abs(PosX) < epsilon)                 # Se seleccionan los e- que cumplen estar por debajo de la E nominal y con X = 0
        eneInd      = isElectron * ( (PosY + 30.0) * bins / 60.0).astype(int)                            # Vector de bins de posicion, se suma 30 cm  para que no haya numeros negativos y 
        sumation    = isElectron * weight * energy                                                            # luego se normaliza al numero de bins deseado   
    elif(axis=='y'):
        isElectron    = 1*(type == 2) * (energy < maxEnergy) * (np.abs(PosY) < epsilon)                 # Se seleccionan los e- que cumplen estar por debajo de la E nominal y con X = 0
        eneInd      = isElectron * ( (PosX + 30.0) * bins / 60.0).astype(int)                            # Vector de bins de posicion, se suma 30 cm para que no haya numeros negativos y 
        sumation    = isElectron * weight *energy
    else:
        print 'ERROR'

    for ind, value in enumerate(eneInd):                                                            # Se enumeran los valores del array eneInd y se va sumando a perfil (en el valor 
        perfil[value] += sumation[ind]                                                              # de eneInd) el aporte de sumation en el indice de la enumeracion de eneInd

    return perfil

def CalculaFluenciaEnergeticaTotal(data_path, maxEnergy, bins,axis):                                                  # OPENS IAEA FILE AND CALCULATES PARTICLE SPECTRUM  IN small DISC
    perfil     = np.zeros(bins)                                                                                 # Perfil de intensidad
    phsp_files   = findFiles(data_path)                                                                         # Finds IAEA files to work with
    for phsp in phsp_files:
        particles, packFormat = parseFile(phsp)                                                                 # Specifies format and amount of particles to read
        
        with open(phsp + '.IAEAphsp', 'rb') as IAEAphsp:                                                        # Opens binary file to read data
            for x in range(particles/100000):
                pack        = np.fromfile(IAEAphsp, dtype=packFormat, count=100000)                             # COUNT = NUMBER OF ITEMS TO READ!!!
                perfil   += CalculaPerfilIntensidad_Pack(pack['Data'], pack['Type'], maxEnergy, bins,axis)           # Calcula la parte del perfil en el archivo

    return perfil   

def PosXY_Fluenciaenergetica_Pack(data, type, maxEnergy):                                                     #  POS X Y DE LAS PARTICULAS DE PACK DATA
    PosXparticula    = data[:, 1]                                                                  
    PosYparticula    = data[:, 2]   
    energy           = np.abs(data[:, 0])
    weight           = data[:, 5] 

    isElectron  = 1*(type == 2) * (energy < maxEnergy)
    
    PosX = PosXparticula[np.nonzero(isElectron)]                                            # El array PosX contiene la posicion X de las particulas que cumplen isElectron == 1
    PosY = PosYparticula[np.nonzero(isElectron)]

    weight = weight[np.nonzero(isElectron)]* energy[np.nonzero(isElectron)]

    return PosX, PosY, weight

def PuntosXY_Fluenciaenergetica_particulas(data_path, maxEnergy):                                    # ABRE EL ARCHIVO IAEA Y DEVUELVE DOS VECTORES CON POS X E Y DE LAS PARTICULAS  IN small DISC
    VectorX     = [0]
    VectorY     = [0]      
    W           = [0]                                                                 
    phsp_files   = findFiles(data_path)                                                  # Finds IAEA files to work with
    
    for phsp in phsp_files:
        particles, packFormat = parseFile(phsp)                                                     # Specifies format and amount of particles to read
    
        bigpack     = 100000
        smallpack   = particles % bigpack

        with open(phsp + '.IAEAphsp', 'rb') as IAEAphsp:                                            # Opens binary file to read data
            for x in range(particles/(bigpack)):
                pack        = np.fromfile(IAEAphsp, dtype=packFormat, count=bigpack)             # Gets full pack data from one file Cuanto mas chico el count menos electrones considera
                PosX,PosY, weight   = PosXY_Fluenciaenergetica_Pack(pack['Data'], pack['Type'], maxEnergy) 
                VectorX = np.concatenate([VectorX,PosX])
                VectorY = np.concatenate([VectorY,PosY])
                W = np.concatenate([W,weight])

            pack        = np.fromfile(IAEAphsp, dtype=packFormat, count=smallpack)             # Gets full pack data from one file Cuanto mas chico el count menos electrones considera
            PosX,PosY, weight   = PosXY_Fluenciaenergetica_Pack(pack['Data'], pack['Type'], maxEnergy) 
            VectorX = np.concatenate([VectorX,PosX])
            VectorY = np.concatenate([VectorY,PosY])
            W = np.concatenate([W,weight])
            print W.max(),W.min()
    return VectorX, VectorY, W

#---------------------------------------------------------------------------- MAIN --------------------------------------------------------------------------------------

# Perfil_X = CalculaPerfilIntensidadTotal('C:\Users\Install\Documents\EGSnrc\Phase Space/', 6.0, 121,'x')
# Perfil_Y = CalculaPerfilIntensidadTotal('C:\Users\Install\Documents\EGSnrc\Phase Space/', 6.0, 121,'y')

# GuardaPerfilIntensidad( Perfil_X, 'C:\Users\Install\Documents\EGSnrc\Phase Space/Perfil_fluencia_X.txt')
# GuardaPerfilIntensidad( Perfil_Y, 'C:\Users\Install\Documents\EGSnrc\Phase Space/Perfil_fluencia_Y.txt')

# ParticulasX, ParticulasY, weight = PuntosXYparticulas('C:\Users\Install\Documents\EGSnrc\Phase Space/', 6.0)
# GuardaPosicionXY(ParticulasX, 
#                  ParticulasY,
#                     weight, 
#                     'C:\Users\Install\Documents\EGSnrc\Phase Space/PosicionesX_fluencia.txt', 
#                     'C:\Users\Install\Documents\EGSnrc\Phase Space/PosicionesY_fluencia.txt',
#                     'C:\Users\Install\Documents\EGSnrc\Phase Space/Pesos_fluencia.txt')


# Perfil_X = AbreVectorTXT('C:\Users\Install\Documents\EGSnrc\Phase Space/Perfil_fluencia_X.txt')
# Perfil_Y = AbreVectorTXT('C:\Users\Install\Documents\EGSnrc\Phase Space/Perfil_fluencia_Y.txt')

# x = AbreVectorTXT('C:\Users\Install\Documents\EGSnrc\Phase Space/PosicionesX_fluencia.txt')
# y = AbreVectorTXT('C:\Users\Install\Documents\EGSnrc\Phase Space/PosicionesY_fluencia.txt')
# weight = AbreVectorTXT('C:\Users\Install\Documents\EGSnrc\Phase Space/Pesos_fluencia.txt')

# bins = 121

# path = 'C:\Users\Install\Documents\EGSnrc\Phase Space/'

# Grafica_PosXY_perfiles(x, y, weight, Perfil_X, Perfil_Y, bins, path)

#----------------------------------------------------------------------------- 30/11 ---------------------------------------------------------------------------+

# Perfil_X = CalculaFluenciaEnergeticaTotal('C:\Users\Install\Documents\EGSnrc\Phase Space/', 6.0, 121,'x')
# Perfil_Y = CalculaFluenciaEnergeticaTotal('C:\Users\Install\Documents\EGSnrc\Phase Space/', 6.0, 121,'y')

# GuardaPerfilIntensidad( Perfil_X, 'C:\Users\Install\Documents\EGSnrc\Phase Space/Perfil_fluenciaenergetica_X.txt')
# GuardaPerfilIntensidad( Perfil_Y, 'C:\Users\Install\Documents\EGSnrc\Phase Space/Perfil_fluenciaenergetica_Y.txt')

# ParticulasX, ParticulasY, weight = PuntosXY_Fluenciaenergetica_particulas('C:\Users\Install\Documents\EGSnrc\Phase Space/', 6.0)
# GuardaPosicionXY(ParticulasX, 
#                  ParticulasY,
#                     weight, 
#                     'C:\Users\Install\Documents\EGSnrc\Phase Space/PosicionesX_fluenciaenergetica.txt', 
#                     'C:\Users\Install\Documents\EGSnrc\Phase Space/PosicionesY_fluenciaenergetica.txt',
#                     'C:\Users\Install\Documents\EGSnrc\Phase Space/Pesos_fluenciaenergetica.txt')


Perfil_X = AbreVectorTXT('C:\Users\Install\Documents\EGSnrc\Phase Space/Perfil_fluenciaenergetica_X.txt')
Perfil_Y = AbreVectorTXT('C:\Users\Install\Documents\EGSnrc\Phase Space/Perfil_fluenciaenergetica_Y.txt')

print Perfil_X.max(), Perfil_Y.max()


x = AbreVectorTXT('C:\Users\Install\Documents\EGSnrc\Phase Space/PosicionesX_fluenciaenergetica.txt')
y = AbreVectorTXT('C:\Users\Install\Documents\EGSnrc\Phase Space/PosicionesY_fluenciaenergetica.txt')
weight = AbreVectorTXT('C:\Users\Install\Documents\EGSnrc\Phase Space/Pesos_fluenciaenergetica.txt')

bins = 121

path = 'C:\Users\Install\Documents\EGSnrc\Phase Space/'

Grafica_PosXY_perfiles(x, y, weight, Perfil_X, Perfil_Y, bins, path)