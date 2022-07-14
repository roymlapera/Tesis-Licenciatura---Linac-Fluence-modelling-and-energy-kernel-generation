import numpy 		as np
from glob           import glob 
from os.path        import splitext, split, dirname
import matplotlib.pyplot as plt


def findFiles(phsp_path):                                                                			# FINDS FILES INSIDE FOLDER WITH IAEAphsp EXTENSION
    find_names  = phsp_path + '/*.IAEAheader'                    									# Names are assumed to be .......... .IAEAheader
    all_names   = glob(find_names)                                                                  # Found files with the searched names
    
    phsp_files  = []                                                                                # List to put correct names
    for name in all_names: phsp_files.append(splitext(name)[0])                                     # Appends name without extension to files
        
    return phsp_files

def parseFile(phsp_file):                                                                           # PARSES IAEA FILE FOR PARTICLES AND PACK FORMAT
    phsp_header = phsp_file + '.IAEAheader'                                                         # Header file with full information

    with open(phsp_header, 'r') as fop: lines = fop.readlines()                                     # Opens file to operate. IT will closed automatically
    for index, line in enumerate(lines):                                                            # Search in every line for keywords

        if ('$BYTE_ORDER' in line) and ('1234' in lines[index + 1]):                               	# Determins the type of unpack format
            packFormat  = np.dtype([('Type', '<i1'), ('Data', '<f4', (7)), ('Extra', '<i4', (2))])   # Assuming IAEA format with ZLAST
        if ('$BYTE_ORDER' in line) and ('4321' in lines[index + 1]):
            packFormat  = np.dtype([('Type', '>i1'), ('Data', '>f4', (7)), ('Extra', '>i4', (2))])   # Assuming IAEA format with ZLAST
        if '$PARTICLES' in line:  particles = int(lines[index + 1])                                 # Total number of particles in file
                
    return particles, packFormat


#Calcula histograma de cosenos directores u y v de un tipo de particulas para determinada energia de referencia

def CalculaHistogramaCosenos_Pack(data, type, maxEnergy, bins, energy_ref):                                          # CALCULATES FULL SPECTRUM FROM PACK DATA
    spectrum    = np.zeros(bins)                                                                       # Espectro de coseno director de ELECTRONES
    
    energy      = np.abs(data[:, 0])
    PosX        = data[:, 1] 
    PosY        = data[:, 2]
    weight      = data[:, 5]

    epsilon = 0.3 																			 # Weight of every photon in Phase plane of current unpack

    #Coseno director de X:

    # xCos = data[:,3]
    # isElectron    = 1*(type == 2) * (energy < maxEnergy) * (energy_ref - epsilon <energy) * (energy<energy_ref + epsilon)								
    # eneInd      = isElectron * ((xCos+1.0) * bins /2.0).astype(int)                               # Index for spectrum sumation (el 2 es porque el coseno director puede variar de -1 a 1)
    # sumation    = isElectron * weight

    #Coseno director de Y:

    yCos        = data[:, 4]                                                                    
    isElectron    = 1*(type == 2) * (energy < maxEnergy) * (energy_ref - epsilon <energy) * (energy<energy_ref + epsilon)
    eneInd      = isElectron * ((yCos + 1.0) * bins / 2.0).astype(int)							 # Index for spectrum sumation (el 2 es porque el coseno director puede variar de -1 a 1)                                 
    sumation    = isElectron * weight

    #Coseno director en Z:

    # xCos        = data[:, 3]
    # yCos        = data[:, 4]
    # zCos        = (1.0 - (data[:, 3])**2 - (data[:, 4])**2)**.5                                                                 
    # isElectron    = 1*(type == 2) * (energy < maxEnergy) 												# Only photons inside calculation limit are proccessed
    # eneInd      = isElectron * ( zCos * bins ).astype(int)                               			 # Index for spectrum sumation (el 2 es porque el coseno director puede variar de -1 a 1)
    # sumation    = isElectron * weight


    for ind, value in enumerate(eneInd): 
    	spectrum[value] += sumation[ind]                        									   # Spectrum photon countingt 

    return spectrum

def CalculaHistogramaTotalCosenos(data_path, maxEnergy, bins, energy_ref):                                    # OPENS IAEA FILE AND CALCULATES PARTICLE SPECTRUM  IN small DISC
    spectrum     = np.zeros(bins)                                                                      # Energy spectrum (E)
    phsp_files   = findFiles(data_path)                                                  # Finds IAEA files to work with
    for phsp in phsp_files:
        particles, packFormat = parseFile(phsp)                                                     # Specifies format and amount of particles to read
        
        with open(phsp + '.IAEAphsp', 'rb') as IAEAphsp:                                            # Opens binary file to read data
            for x in range(particles/100000):
                pack        = np.fromfile(IAEAphsp, dtype=packFormat, count=100000)                   # COUNT = NUMBER OF ITEMS TO READ!!!
                spectrum   += CalculaHistogramaCosenos_Pack(pack['Data'], pack['Type'], maxEnergy, bins, energy_ref)            # Calculates spectrum part from one file

    # spectrum /= spectrum.sum()        															#NO ESTA NORMALIZANDO

    return spectrum  																		# Spectrum tiene los valores del eje Y del histograma en funcion del numero de bin 


#Calcula el espectro de cosenos directores u y v de un tipo de particulas EN UNA DETERMINADA REGION

def CalculaHistogramaCosenos_Pack_por_region(data, type, maxEnergy, bins):                                          # CALCULATES FULL SPECTRUM FROM PACK DATA
    spectrum    = np.zeros(bins)                                                                       # Espectro de coseno director de ELECTRONES
    
    energy      = np.abs(data[:, 0])
    PosX        = data[:, 1] 
    PosY        = data[:, 2]
    weight      = data[:, 5] 																			 # Weight of every photon in Phase plane of current unpack

    #Comentar y descomentar esto segun si uno quiere coseno director en la direccion X, Y o Z

    #Coseno director de X:

    # xCos = data[:,3]
    # isElectron    = 1*(type == 2) * (energy < maxEnergy)  * (-25.0<PosX) * (PosX<-15.0) * (15.0<PosY) * (PosY<25.0)									
    # eneInd      = isElectron * ((xCos+1.0) * bins /2.0).astype(int)                               # Index for spectrum sumation (el 2 es porque el coseno director puede variar de -1 a 1)
    # sumation    = isElectron * weight

    #Coseno director de Y:

    yCos        = data[:, 4]                                                                    
    isElectron    = 1*(type == 2) * (energy < maxEnergy) * (-25.0<PosX) * (PosX<-15.0) * (15.0<PosY) * (PosY<25.0)	
    eneInd      = isElectron * ((yCos + 1.0) * bins / 2.0).astype(int)							 # Index for spectrum sumation (el 2 es porque el coseno director puede variar de -1 a 1)                                 
    sumation    = isElectron * weight

    #Coseno director en Z:

    # xCos        = data[:, 3]
    # yCos        = data[:, 4]
    # zCos        = (1.0 - (data[:, 3])**2 - (data[:, 4])**2)**.5                                                                 
    # isElectron    = 1*(type == 2) * (energy < maxEnergy) 												# Only photons inside calculation limit are proccessed
    # eneInd      = isElectron * ( zCos * bins ).astype(int)                               			 # Index for spectrum sumation (el 2 es porque el coseno director puede variar de -1 a 1)
    # sumation    = isElectron * weight


    for ind, value in enumerate(eneInd): 
    	spectrum[value] += sumation[ind]                        									   # Spectrum photon countingt 

    return spectrum

def CalculaHistogramaTotalCosenos_por_region(data_path, maxEnergy, bins):                                    # OPENS IAEA FILE AND CALCULATES PARTICLE SPECTRUM  IN small DISC
    spectrum     = np.zeros(bins)                                                                      # Energy spectrum (E)
    phsp_files   = findFiles(data_path)                                                  # Finds IAEA files to work with
    for phsp in phsp_files:
        particles, packFormat = parseFile(phsp)                                                     # Specifies format and amount of particles to read
        
        with open(phsp + '.IAEAphsp', 'rb') as IAEAphsp:                                            # Opens binary file to read data
            for x in range(particles/100000):
                pack        = np.fromfile(IAEAphsp, dtype=packFormat, count=100000)                   # COUNT = NUMBER OF ITEMS TO READ!!!
                spectrum   += CalculaHistogramaCosenos_Pack_por_region(pack['Data'], pack['Type'], maxEnergy, bins)            # Calculates spectrum part from one file

    spectrum /= spectrum.sum()

    return spectrum  																		# Spectrum tiene los valores del eje Y del histograma en funcion del numero de bin 


#Funciones que guardan y grafican histogramas

def GuardaHistogramaTotalCosenos(data, file):
	np.savetxt(file, data)

def GraficaHistogramaTotalCosenos_por_region(data, data_region_total, bins):																# CAMBIAR ACA LA LEYENDA SEGUN X E Y

	abcisas = np.arange(bins)

	abcisas = abcisas *2.0 /float(bins-1) - 1.0

	#Para zCos
	#abcisas = abcisas /float(bins-1)

	#fig = plt.figure()

	opacity = 0.5

	bar_width = 0.02

	#plt.bar(abcisas, data, width = 0.09,color='green', align='center', alpha=0.5, label='Energia = 2.5 MeV')

	fig, ax = plt.subplots()

	rects1 = plt.bar(abcisas, data, bar_width,
                 alpha=opacity,
                 color='b',
                 label='Region 1')
 
	rects2 = plt.bar(abcisas + bar_width, data_region_total, bar_width,
                 alpha=opacity,
                 color='g',
                 label='Region Total')

	plt.legend()

	plt.axis([-1.1, 1.1, 0, 0.2])

	plt.xlabel('Coseno director en X', fontsize=14)

	plt.grid(True)

	plt.show()

def GraficaHistogramaTotalCosenos_dif_energias(data1, data2, data3, bins):																# CAMBIAR ACA LA LEYENDA SEGUN X E Y

	abcisas = np.arange(bins)

	abcisas = abcisas *2.0 /float(bins-1) - 1.0

	#Para zCos
	#abcisas = abcisas /float(bins-1)

	#fig = plt.figure()

	opacity = 0.5

	bar_width = 0.02

	#plt.bar(abcisas, data, width = 0.09,color='green', align='center', alpha=0.5, label='Energia = 2.5 MeV')

	fig, ax = plt.subplots()

	rects1 = plt.bar(abcisas, data1, bar_width,
                 alpha=opacity,
                 color='b',
                 label='Energia = 0.46 MeV')
 
	rects2 = plt.bar(abcisas + bar_width, data2, bar_width,
                 alpha=opacity,
                 color='g',
                 label='Energia = 2.5 MeV')

	rects3 = plt.bar(abcisas + bar_width + bar_width, data3, bar_width,
                 alpha=opacity,
                 color='r',
                 label='Energia = 5 MeV')

	plt.legend()

	plt.axis([-1.1, 1.1, 0, 12])

	plt.xlabel('Coseno director en Y', fontsize=14)

	plt.grid(True)

	plt.show()


#Funcion que abre un histograma guardado en un archivo txt y devuelve un vector con el contenido del histograma

def AbreHistogramaTXT(file):
	with open(file) as f:
  		histo = np.loadtxt(f, delimiter="\n", dtype='double', comments="#", skiprows=0, usecols=None)
  	return histo

# -------------------------------------------------------------------------- MAIN ------------------------------------------------------------------

# HISTOGRAMAS COSENOS DIRECTORES PARA DIFERENTES ENERGIAS


# Hace un histograma de coseno director en X y los guarda  #Recordar que hay que cambiar CalculaHistogramaCosenos_Pack segun X o Y

# Histograma1 = CalculaHistogramaTotalCosenos('C:\Users\Install\Documents\EGSnrc\Phase Space/', 6.0, 20, 0.46)
# Histograma2 = CalculaHistogramaTotalCosenos('C:\Users\Install\Documents\EGSnrc\Phase Space/', 6.0, 20, 2.5)
# Histograma3 = CalculaHistogramaTotalCosenos('C:\Users\Install\Documents\EGSnrc\Phase Space/', 6.0, 20, 5.0)

# GuardaHistogramaTotalCosenos(Histograma1, 'Histograma_CosDir_X_0.46MeV.txt')
# GuardaHistogramaTotalCosenos(Histograma2, 'Histograma_CosDir_X_2.5MeV.txt')
# GuardaHistogramaTotalCosenos(Histograma3, 'Histograma_CosDir_X_5.0MeV.txt')

# Hace un histograma de coseno director en Y y los guarda   #Recordar que hay que cambiar CalculaHistogramaCosenos_Pack segun X o Y

# Histograma1 = CalculaHistogramaTotalCosenos('C:\Users\Install\Documents\EGSnrc\Phase Space/', 6.0, 20, 0.46)
# Histograma2 = CalculaHistogramaTotalCosenos('C:\Users\Install\Documents\EGSnrc\Phase Space/', 6.0, 20, 2.5)
# Histograma3 = CalculaHistogramaTotalCosenos('C:\Users\Install\Documents\EGSnrc\Phase Space/', 6.0, 20, 5.0)

# GuardaHistogramaTotalCosenos(Histograma1, 'Histograma_CosDir_Y_0.46MeV.txt')
# GuardaHistogramaTotalCosenos(Histograma2, 'Histograma_CosDir_Y_2.5MeV.txt')
# GuardaHistogramaTotalCosenos(Histograma3, 'Histograma_CosDir_Y_5.0MeV.txt')


#Abre histogramas de cosenos directores guardados en carpeta 

# Histograma_X_0_46MeV = AbreHistogramaTXT('Histograma_CosDir_X_0.46MeV.txt')
# Histograma_X_2_5MeV = AbreHistogramaTXT('Histograma_CosDir_X_2.5MeV.txt')
# Histograma_X_5_0MeV = AbreHistogramaTXT('Histograma_CosDir_X_5.0MeV.txt')

# Histograma_Y_0_46MeV = AbreHistogramaTXT('Histograma_CosDir_Y_0.46MeV.txt')
# Histograma_Y_2_5MeV = AbreHistogramaTXT('Histograma_CosDir_Y_2.5MeV.txt')
# Histograma_Y_5_0MeV = AbreHistogramaTXT('Histograma_CosDir_Y_5.0MeV.txt')


#Hace histograma de cosenos directores en X e Y para diferentes energias a partir de los archivos txt

#GraficaHistogramaTotalCosenos_dif_energias(Histograma_Y_0_46MeV, Histograma_Y_2_5MeV, Histograma_Y_5_0MeV, 20)			#Cambiar el titulo del grafico a Coseno Director segun sea en X o Y



#HISTOGRAMAS COSENOS DIRECTORES PARA DIFERENTES REGIONES


#Hace un histograma de coseno director en X EN TODA LA REGION (de -30 cm a 30 cm en X e Y) y los guarda  #Recordar que hay que cambiar CalculaHistogramaCosenos_Pack_por_region segun X o Y

# HistogramaRegionTotal_X = CalculaHistogramaTotalCosenos_por_region('C:\Users\Install\Documents\EGSnrc\Phase Space/', 6.0, 40)  #Para hacer esto, sacar las restriccion espacial en isElectron
# GuardaHistogramaTotalCosenos(HistogramaRegionTotal_X, 'Histograma_CosDir_X_Reg_Total.txt')

#Hace un histograma de coseno director en Y EN TODA LA REGION (de -30 cm a 30 cm en X e Y) y los guarda  #Recordar que hay que cambiar CalculaHistogramaCosenos_Pack_por_region segun X o Y

# HistogramaRegionTotal_Y = CalculaHistogramaTotalCosenos_por_region('C:\Users\Install\Documents\EGSnrc\Phase Space/', 6.0, 40)  #Para hacer esto, sacar las restriccion espacial en isElectron
# GuardaHistogramaTotalCosenos(HistogramaRegionTotal_Y, 'Histograma_CosDir_Y_Reg_Total.txt')


#Hace un histograma de coseno director en X EN LAS DISTINTAS REGIONES (1, 2 y 3) y los guarda  #Recordar que hay que cambiar CalculaHistogramaCosenos_Pack_por_region segun X o Y

# HistogramaReg1_X = CalculaHistogramaTotalCosenos_por_region('C:\Users\Install\Documents\EGSnrc\Phase Space/', 6.0, 40)  #Para hacer esto, sacar las restriccion espacial en isElectron
# GuardaHistogramaTotalCosenos(HistogramaReg1_X, 'Histograma_CosDir_X_Reg1.txt')

# HistogramaReg2_X = CalculaHistogramaTotalCosenos_por_region('C:\Users\Install\Documents\EGSnrc\Phase Space/', 6.0, 40)  #Para hacer esto, sacar las restriccion espacial en isElectron
# GuardaHistogramaTotalCosenos(HistogramaReg2_X, 'Histograma_CosDir_X_Reg2.txt')

# HistogramaReg3_X = CalculaHistogramaTotalCosenos_por_region('C:\Users\Install\Documents\EGSnrc\Phase Space/', 6.0, 40)  #Para hacer esto, sacar las restriccion espacial en isElectron
# GuardaHistogramaTotalCosenos(HistogramaReg3_X, 'Histograma_CosDir_X_Reg3.txt')

#Hace un histograma de coseno director en Y EN LAS DISTINTAS REGIONES (1, 2 y 3) y los guarda  #Recordar que hay que cambiar CalculaHistogramaCosenos_Pack segun X o Y

# HistogramaReg1_Y = CalculaHistogramaTotalCosenos_por_region('C:\Users\Install\Documents\EGSnrc\Phase Space/', 6.0, 40)  #Para hacer esto, sacar las restriccion espacial en isElectron
# GuardaHistogramaTotalCosenos(HistogramaReg1_Y, 'Histograma_CosDir_Y_Reg1.txt')

# HistogramaReg2_Y = CalculaHistogramaTotalCosenos_por_region('C:\Users\Install\Documents\EGSnrc\Phase Space/', 6.0, 40)  #Para hacer esto, sacar las restriccion espacial en isElectron
# GuardaHistogramaTotalCosenos(HistogramaReg2_Y, 'Histograma_CosDir_Y_Reg2.txt')

# HistogramaReg3_Y = CalculaHistogramaTotalCosenos_por_region('C:\Users\Install\Documents\EGSnrc\Phase Space/', 6.0, 40)  #Para hacer esto, sacar las restriccion espacial en isElectron
# GuardaHistogramaTotalCosenos(HistogramaReg3_Y, 'Histograma_CosDir_Y_Reg3.txt')


#Abre histogramas de cosenos directores guardados en carpeta 

HistogramaRegTotal_X = AbreHistogramaTXT('C:\Users\Install\Documents\EGSnrc\Phase Space\Histograma_CosDir_X_Reg_Total.txt')
HistogramaReg1_X = AbreHistogramaTXT('C:\Users\Install\Documents\EGSnrc\Phase Space\Histograma_CosDir_X_Reg1.txt')
HistogramaReg2_X = AbreHistogramaTXT('C:\Users\Install\Documents\EGSnrc\Phase Space\Histograma_CosDir_X_Reg2.txt')
HistogramaReg3_X = AbreHistogramaTXT('C:\Users\Install\Documents\EGSnrc\Phase Space\Histograma_CosDir_X_Reg3.txt')

# HistogramaRegTotal_Y = AbreHistogramaTXT('C:\Users\Install\Documents\EGSnrc\Phase Space\Cosenos directores por region\Histograma_CosDir_Y_RegTotal.txt')
# HistogramaReg1_Y = AbreHistogramaTXT('C:\Users\Install\Documents\EGSnrc\Phase Space\Cosenos directores por region\Histograma_CosDir_Y_Reg1.txt')
# HistogramaReg2_Y = AbreHistogramaTXT('C:\Users\Install\Documents\EGSnrc\Phase Space\Cosenos directores por region\Histograma_CosDir_Y_Reg2.txt')
# HistogramaReg3_Y = AbreHistogramaTXT('C:\Users\Install\Documents\EGSnrc\Phase Space\Cosenos directores por region\Histograma_CosDir_Y_Reg3.txt')

#Hace histograma de cosenos directores en X e Y para diferentes energias a partir de los archivos txt

GraficaHistogramaTotalCosenos_por_region(HistogramaReg1_X, HistogramaRegTotal_X, 40)		#Cambiar el titulo del grafico a Coseno Director segun sea en X o Y
