import numpy 		as np
from glob           import glob 
from os.path        import splitext, split, dirname
import matplotlib.pyplot as plt


def findFiles(phsp_path):                                                                           # FINDS FILES INSIDE FOLDER WITH IAEAphsp EXTENSION
    find_names  = phsp_path + '/*.IAEAphsp'                                                       # Names are assumed to be .......... .IAEAheader
    all_names   = glob(find_names)                                                                  # Found files with the searched names
    
    phsp_files  = []                                                                                # List to put correct names
    for name in all_names: phsp_files.append(splitext(name)[0])                                     # Appends name without extension to files
    
    print 'Espacios de fases encontrados:'

    print '\n'.join(phsp_files)

    phsp_files = np.delete(phsp_files,[1,2,3,4,5,6])

    print 'Abriendo:'

    print '\n'.join(phsp_files)

    return phsp_files

def parseFile(phsp_file):                                                                           # PARSES IAEA FILE FOR PARTICLES AND PACK FORMAT
    phsp_header = phsp_file + '.IAEAheader'                                                         # Header file with full information

    if (phsp_file == 'C:\Users\Install\Documents\EGSnrc\Phase Space\VarianClinaciX_6MV_20x20_w1' or 'C:\Users\Install\Documents\EGSnrc\Phase Space\VarianClinaciX_6MV_20x20_Homogeneous'):
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

#Calcula histograma de cosenos directores u y v de un tipo de particulas para determinada energia de referencia

def CalculaHistogramaCosenos_Pack(data, type, maxEnergy, bins, energy_ref):                         # CALCULATES FULL SPECTRUM FROM PACK DATA
    histo       = np.zeros(bins)                                                                       # Espectro de coseno director de ELECTRONES
    
    energy      = np.abs(data[:, 0])
    weight      = data[:, 5]

    epsilon = 0.3 																			 # Weight of every photon in Phase plane of current unpack

    #Coseno director de X:

    xCos = data[:,3]
    isElectron    = 1*(type == 2) * (energy < maxEnergy) * ( abs(energy - energy_ref) < epsilon) * (abs(xCos) < 1.0)								
    eneInd      = isElectron * ((xCos+1.0) * bins /2.0).astype(int)                               # Index for spectrum sumation (el 2 es porque el coseno director puede variar de -1 a 1)
    sumation    = isElectron * weight

    #Coseno director de Y:

    # yCos        = data[:, 4]                                                                    
    # isElectron    = 1*(type == 2) * (energy < maxEnergy) * ( abs(energy - energy_ref) < epsilon) * (abs(yCos) < 1.0)  
    # eneInd      = isElectron * ((yCos + 1.0) * bins / 2.0).astype(int)							 # Index for spectrum sumation (el 2 es porque el coseno director puede variar de -1 a 1)                                 
    # sumation    = isElectron * weight

    #Coseno director en Z:

    # xCos        = data[:, 3]
    # yCos        = data[:, 4]
    # zCos        = (1.0 - (data[:, 3])**2 - (data[:, 4])**2)**.5                                                                 
    # isElectron    = 1*(type == 2) * (energy < maxEnergy) 												# Only photons inside calculation limit are proccessed
    # eneInd      = isElectron * ( zCos * bins ).astype(int)                               			 # Index for spectrum sumation (el 2 es porque el coseno director puede variar de -1 a 1)
    # sumation    = isElectron * weight

    for ind, value in enumerate(eneInd): 
      	histo[value] += sumation[ind]                        									   # Spectrum photon countingt 

    return histo

def CalculaHistogramaTotalCosenos(data_path, maxEnergy, bins, energy_ref):                          # OPENS IAEA FILE AND CALCULATES PARTICLE SPECTRUM  IN small DISC
    histo     = np.zeros(bins)                                                                      # Energy spectrum (E)
    phsp_files   = findFiles(data_path)                                                             # Finds IAEA files to work with
    for phsp in phsp_files:
        particles, packFormat = parseFile(phsp)                                                     # Specifies format and amount of particles to read
        
        with open(phsp + '.IAEAphsp', 'rb') as IAEAphsp:                                            # Opens binary file to read data
            for x in range(particles/100000):
                pack        = np.fromfile(IAEAphsp, dtype=packFormat, count=100000)                 # COUNT = NUMBER OF ITEMS TO READ!!!
                histo       += CalculaHistogramaCosenos_Pack(pack['Data'], pack['Type'], maxEnergy, bins, energy_ref)                                                              # Calculates spectrum part from one file

    return histo  																		# Spectrum tiene los valores del eje Y del histograma en funcion del numero de bin 


#Calcula el espectro de cosenos directores u y v de un tipo de particulas EN UNA DETERMINADA REGION

def CalculaHistogramaCosenos_Pack_por_region(data, type, maxEnergy, bins):                                          # CALCULATES FULL SPECTRUM FROM PACK DATA
    histo       = np.zeros(bins)                                                                       # Espectro de coseno director de ELECTRONES
    
    energy      = np.abs(data[:, 0])
    PosX        = data[:, 1] 
    PosY        = data[:, 2]
    weight      = data[:, 5] 																			 # Weight of every photon in Phase plane of current unpack

    #Comentar y descomentar esto segun si uno quiere coseno director en la direccion X, Y o Z

    #Coseno director de X:

    xCos = data[:,3]
    isElectron    = 1*(type == 2) * (energy < maxEnergy) #* (-25.0<PosX) * (PosX<-15.0) * (15.0<PosY) * (PosY<25.0)									
    cosInd      = isElectron * ((xCos+1.0) * bins /2.0).astype(int)                               # Index for spectrum sumation (el 2 es porque el coseno director puede variar de -1 a 1)
    sumation    = isElectron * weight

    #Coseno director de Y:

    # yCos        = data[:, 4]                                                                    
    # isElectron    = 1*(type == 2) * (energy < maxEnergy) #* (-25.0<PosX) * (PosX<-15.0) * (15.0<PosY) * (PosY<25.0)
    # cosInd      = isElectron * ((yCos + 1.0) * bins / 2.0).astype(int)							 # Index for spectrum sumation (el 2 es porque el coseno director puede variar de -1 a 1)                                 
    # sumation    = isElectron * weight




    #Angulo director de X:

    # xCos = data[:,3]
    # isElectron    = 1*(type == 2) * (energy < maxEnergy) #* (-25.0<PosX) * (PosX<-15.0) * (15.0<PosY) * (PosY<25.0)                                 
    # cosInd      = isElectron * (np.arccos(xCos) * bins /(np.pi)).astype(int)                               # Index for spectrum sumation (el 2 es porque el coseno director puede variar de -1 a 1)
    # sumation    = isElectron * weight

    #Angulo director de Y:

    # yCos        = data[:, 4]                                                                    
    # isElectron    = 1*(type == 2) * (energy < maxEnergy) #* (-25.0<PosX) * (PosX<-15.0) * (15.0<PosY) * (PosY<25.0)
    # cosInd      = isElectron * (np.arccos(yCos) * bins /np.pi).astype(int)                               # Index for spectrum sumation (el 2 es porque el coseno director puede variar de -1 a 1)                              
    # sumation    = isElectron * weight

    cosInd = cosInd * (cosInd<100)

    # print cosInd[cosInd>99]


    for ind, value in enumerate(cosInd): 
    	histo[value] += sumation[ind]                        									   # Spectrum photon countingt 

    return histo

def CalculaHistogramaTotalCosenos_por_region(data_path, maxEnergy, bins):                                    # OPENS IAEA FILE AND CALCULATES PARTICLE SPECTRUM  IN small DISC
    histo    = np.zeros(bins)                                                                      # Energy spectrum (E)
    phsp_files   = findFiles(data_path)                                                  # Finds IAEA files to work with
    for phsp in phsp_files:
        particles, packFormat = parseFile(phsp)                                                     # Specifies format and amount of particles to read
        
        with open(phsp + '.IAEAphsp', 'rb') as IAEAphsp:                                            # Opens binary file to read data
            for x in range(particles/100000):
                pack        = np.fromfile(IAEAphsp, dtype=packFormat, count=100000)                   # COUNT = NUMBER OF ITEMS TO READ!!!
                histo   += CalculaHistogramaCosenos_Pack_por_region(pack['Data'], pack['Type'], maxEnergy, bins)            # Calculates spectrum part from one file

    return histo  																		# Spectrum tiene los valores del eje Y del histograma en funcion del numero de bin 


#Funciones que guardan y grafican histogramas

def GuardaHistogramaTotalCosenos(data, file):
    np.savetxt(file, data)

def GraficaHistogramaTotalCosenos_por_region(data, data_region_total, bins):																# CAMBIAR ACA LA LEYENDA SEGUN X E Y

    path = 'C:\Users\Install\Desktop\Tesis Roy\cColSec\Cosenos directores de electrones\Por region/'

    opacity = 0.8

    abcisas = np.arange(bins)

    abcisas2 = abcisas * np.pi / float(bins-1)

    abcisas = abcisas *2.0 /float(bins-1) - 1.0



    #Normalizacion:
    data /= data.sum()
    data_region_total /= data_region_total.sum()


    abcisas = np.ravel(zip(abcisas,abcisas + 2.0 /float(bins-1)))

    abcisas = abcisas -1.0 /float(bins-1)

    abcisas2 = np.ravel(zip(abcisas2,abcisas2 + np.pi /float(bins-1)))

    abcisas2 = abcisas2 -np.pi/2.0 /float(bins-1)

    data = np.ravel(zip(data,data))
    data_region_total = np.ravel(zip(data_region_total,data_region_total))

    fig, ax = plt.subplots()

    rects1 = plt.plot(abcisas, data, alpha=opacity, color='b', label='Homog cosenos', linewidth = 2)
 
    rects2 = plt.plot(abcisas2, data, alpha=opacity, color='g', label='Homog angulos', linewidth = 2)

    plt.legend(fontsize=11)

    # plt.axis([-1.1, 3.5, 0, 0.02])

    plt.xlabel('Coseno director en Y', fontsize=14)

    plt.grid(True)

    # fig.savefig(path + 'Coseno_dir_en_y_reg3.png', dpi = 200, bbox_inches='tight')

    plt.show()

def GraficaHistogramaTotalCosenos_dif_energias(data1, data2, data3, bins):																# CAMBIAR ACA LA LEYENDA SEGUN X E Y

    path = 'C:\Users\Install\Desktop\Tesis Roy\cColSec\Cosenos directores de electrones\Para diferentes energias/'

    opacity = 0.8

    abcisas = np.arange(bins)

    abcisas = abcisas *2.0 /float(bins-1) - 1.0

    #Para zCos:
    #abcisas = abcisas /float(bins-1)

    #Normalizacion:

    data1 /= data1.sum()
    data2 /= data2.sum()
    data3 /= data3.sum()

    abcisas = np.ravel(zip(abcisas,abcisas + 2.0 /float(bins-1)))

    abcisas = abcisas -1.0 /float(bins-1)

    data1 = np.ravel(zip(data1,data1))
    data2 = np.ravel(zip(data2,data2))
    data3 = np.ravel(zip(data3,data3))

    fig, ax = plt.subplots()

    rects1 = plt.plot(abcisas, data1, alpha=opacity, color='b', label='Energia = 0.46 MeV', linewidth = 2)
 
    rects2 = plt.plot(abcisas, data2, alpha=opacity, color='g', label='Energia = 2.5 MeV', linewidth = 2)

    rects3 = plt.plot(abcisas, data3, alpha=opacity, color='r', label='Energia = 5 MeV', linewidth = 2)

    plt.legend(fontsize=11)

    plt.axis([-1.1, 1.1, 0, 0.16])

    plt.xlabel('Coseno director en Y', fontsize=14)

    plt.grid(True)

    fig.savefig(path + 'Coseno_dir_en_y.png', dpi = 200, bbox_inches='tight')

    plt.show()


#Funcion que abre un histograma guardado en un archivo txt y devuelve un vector con el contenido del histograma

def AbreHistogramaTXT(file):
    with open(file) as f:
        histo = np.loadtxt(f, delimiter="\n", dtype='double', comments="#", skiprows=0, usecols=None)
    return histo

#Funcion que calcula histogramas de cosenos directores en X e Y en regiones del espacio de fase


# -------------------------------------------------------------------------- MAIN ------------------------------------------------------------------

# HISTOGRAMAS COSENOS DIRECTORES PARA DIFERENTES ENERGIAS


# Hace un histograma de coseno director en X y los guarda  #Recordar que hay que cambiar CalculaHistogramaCosenos_Pack segun X o Y

# Histograma1 = CalculaHistogramaTotalCosenos('C:\Users\Install\Documents\EGSnrc\Phase Space/', 6.0, 41, 0.46)
# Histograma2 = CalculaHistogramaTotalCosenos('C:\Users\Install\Documents\EGSnrc\Phase Space/', 6.0, 41, 2.5)
# Histograma3 = CalculaHistogramaTotalCosenos('C:\Users\Install\Documents\EGSnrc\Phase Space/', 6.0, 41, 5.0)

# GuardaHistogramaTotalCosenos(Histograma1, 'C:\Users\Install\Documents\EGSnrc\Phase Space\VarianClinaciX_6MV_20x20_w1\Cosenos directores por energia/Histograma_CosDir_X_0.46MeV.txt')
# GuardaHistogramaTotalCosenos(Histograma2, 'C:\Users\Install\Documents\EGSnrc\Phase Space\VarianClinaciX_6MV_20x20_w1\Cosenos directores por energia/Histograma_CosDir_X_2.5MeV.txt')
# GuardaHistogramaTotalCosenos(Histograma3, 'C:\Users\Install\Documents\EGSnrc\Phase Space\VarianClinaciX_6MV_20x20_w1\Cosenos directores por energia/Histograma_CosDir_X_5.0MeV.txt')

# Hace un histograma de coseno director en Y y los guarda   #Recordar que hay que cambiar CalculaHistogramaCosenos_Pack segun X o Y

# Histograma1 = CalculaHistogramaTotalCosenos('C:\Users\Install\Documents\EGSnrc\Phase Space/', 6.0, 41, 0.46)
# Histograma2 = CalculaHistogramaTotalCosenos('C:\Users\Install\Documents\EGSnrc\Phase Space/', 6.0, 41, 2.5)
# Histograma3 = CalculaHistogramaTotalCosenos('C:\Users\Install\Documents\EGSnrc\Phase Space/', 6.0, 41, 5.0)

# GuardaHistogramaTotalCosenos(Histograma1, 'C:\Users\Install\Documents\EGSnrc\Phase Space\VarianClinaciX_6MV_20x20_w1\Cosenos directores por energia/Histograma_CosDir_Y_0.46MeV.txt')
# GuardaHistogramaTotalCosenos(Histograma2, 'C:\Users\Install\Documents\EGSnrc\Phase Space\VarianClinaciX_6MV_20x20_w1\Cosenos directores por energia/Histograma_CosDir_Y_2.5MeV.txt')
# GuardaHistogramaTotalCosenos(Histograma3, 'C:\Users\Install\Documents\EGSnrc\Phase Space\VarianClinaciX_6MV_20x20_w1\Cosenos directores por energia/Histograma_CosDir_Y_5.0MeV.txt')


#Abre histogramas de cosenos directores guardados en carpeta 

# Histograma_X_0_46MeV = AbreHistogramaTXT('C:\Users\Install\Documents\EGSnrc\Phase Space\VarianClinaciX_6MV_20x20_w1\Cosenos directores por energia/Histograma_CosDir_X_0.46MeV.txt')
# Histograma_X_2_5MeV = AbreHistogramaTXT('C:\Users\Install\Documents\EGSnrc\Phase Space\VarianClinaciX_6MV_20x20_w1\Cosenos directores por energia/Histograma_CosDir_X_2.5MeV.txt')
# Histograma_X_5_0MeV = AbreHistogramaTXT('C:\Users\Install\Documents\EGSnrc\Phase Space\VarianClinaciX_6MV_20x20_w1\Cosenos directores por energia/Histograma_CosDir_X_5.0MeV.txt')

# Histograma_Y_0_46MeV = AbreHistogramaTXT('C:\Users\Install\Documents\EGSnrc\Phase Space\VarianClinaciX_6MV_20x20_w1\Cosenos directores por energia/Histograma_CosDir_Y_0.46MeV.txt')
# Histograma_Y_2_5MeV = AbreHistogramaTXT('C:\Users\Install\Documents\EGSnrc\Phase Space\VarianClinaciX_6MV_20x20_w1\Cosenos directores por energia/Histograma_CosDir_Y_2.5MeV.txt')
# Histograma_Y_5_0MeV = AbreHistogramaTXT('C:\Users\Install\Documents\EGSnrc\Phase Space\VarianClinaciX_6MV_20x20_w1\Cosenos directores por energia/Histograma_CosDir_Y_5.0MeV.txt')


# # #Hace histograma de cosenos directores en X e Y para diferentes energias a partir de los archivos txt

# GraficaHistogramaTotalCosenos_dif_energias(Histograma_Y_0_46MeV, Histograma_Y_2_5MeV, Histograma_Y_5_0MeV,41)			#Cambiar el titulo del grafico a Coseno Director segun sea en X o Y



#HISTOGRAMAS COSENOS DIRECTORES PARA DIFERENTES REGIONES


#Hace un histograma de coseno director en X EN TODA LA REGION (de -30 cm a 30 cm en X e Y) y los guarda  #Recordar que hay que cambiar CalculaHistogramaCosenos_Pack_por_region segun X o Y

# HistogramaRegionTotal_X = CalculaHistogramaTotalCosenos_por_region('C:\Users\Install\Documents\EGSnrc\Phase Space/', 6.0, 40)  #Para hacer esto, sacar las restriccion espacial en isElectron
# GuardaHistogramaTotalCosenos(HistogramaRegionTotal_X, 'C:\Users\Install\Documents\EGSnrc\Phase Space\Varian_without_ColSec_2.1\Cosenos directores por region/Histograma_CosDir_X_Reg_Total_sColSec.txt')

#Hace un histograma de coseno director en Y EN TODA LA REGION (de -30 cm a 30 cm en X e Y) y los guarda  #Recordar que hay que cambiar CalculaHistogramaCosenos_Pack_por_region segun X o Y

# HistogramaRegionTotal_Y = CalculaHistogramaTotalCosenos_por_region('C:\Users\Install\Documents\EGSnrc\Phase Space/', 6.0, 40)  #Para hacer esto, sacar las restriccion espacial en isElectron
# GuardaHistogramaTotalCosenos(HistogramaRegionTotal_Y, 'C:\Users\Install\Documents\EGSnrc\Phase Space\Varian_without_ColSec_2.1\Cosenos directores por region/Histograma_CosDir_Y_Reg_Total_sColSec.txt')



#Hace un histograma de coseno director en X EN LAS DISTINTAS REGIONES (1, 2 y 3) y los guarda  #Recordar que hay que cambiar CalculaHistogramaCosenos_Pack_por_region segun X o Y

# HistogramaReg1_X = CalculaHistogramaTotalCosenos_por_region('C:\Users\Install\Documents\EGSnrc\Phase Space/', 6.0, 40)  #Para hacer esto, sacar las restriccion espacial en isElectron
# GuardaHistogramaTotalCosenos(HistogramaReg1_X, 'C:\Users\Install\Documents\EGSnrc\Phase Space\Varian_without_ColSec_2.1\Cosenos directores por region/Histograma_CosDir_X_Reg1.txt')

# HistogramaReg2_X = CalculaHistogramaTotalCosenos_por_region('C:\Users\Install\Documents\EGSnrc\Phase Space/', 6.0, 40)  #Para hacer esto, sacar las restriccion espacial en isElectron
# GuardaHistogramaTotalCosenos(HistogramaReg2_X, 'C:\Users\Install\Documents\EGSnrc\Phase Space\Varian_without_ColSec_2.1\Cosenos directores por region/Histograma_CosDir_X_Reg2.txt')

# HistogramaReg3_X = CalculaHistogramaTotalCosenos_por_region('C:\Users\Install\Documents\EGSnrc\Phase Space/', 6.0, 40)  #Para hacer esto, sacar las restriccion espacial en isElectron
# GuardaHistogramaTotalCosenos(HistogramaReg3_X, 'C:\Users\Install\Documents\EGSnrc\Phase Space\Varian_without_ColSec_2.1\Cosenos directores por region/Histograma_CosDir_X_Reg3.txt')



#Hace un histograma de coseno director en Y EN LAS DISTINTAS REGIONES (1, 2 y 3) y los guarda  #Recordar que hay que cambiar CalculaHistogramaCosenos_Pack segun X o Y

# HistogramaReg1_Y = CalculaHistogramaTotalCosenos_por_region('C:\Users\Install\Documents\EGSnrc\Phase Space/', 6.0, 40)  #Para hacer esto, sacar las restriccion espacial en isElectron
# GuardaHistogramaTotalCosenos(HistogramaReg1_Y, 'C:\Users\Install\Documents\EGSnrc\Phase Space\Varian_without_ColSec_2.1\Cosenos directores por region/Histograma_CosDir_Y_Reg1.txt')

# HistogramaReg2_Y = CalculaHistogramaTotalCosenos_por_region('C:\Users\Install\Documents\EGSnrc\Phase Space/', 6.0, 40)  #Para hacer esto, sacar las restriccion espacial en isElectron
# GuardaHistogramaTotalCosenos(HistogramaReg2_Y, 'C:\Users\Install\Documents\EGSnrc\Phase Space\Varian_without_ColSec_2.1\Cosenos directores por region/Histograma_CosDir_Y_Reg2.txt')

# HistogramaReg3_Y = CalculaHistogramaTotalCosenos_por_region('C:\Users\Install\Documents\EGSnrc\Phase Space/', 6.0, 40)  #Para hacer esto, sacar las restriccion espacial en isElectron
# GuardaHistogramaTotalCosenos(HistogramaReg3_Y, 'C:\Users\Install\Documents\EGSnrc\Phase Space\Varian_without_ColSec_2.1\Cosenos directores por region/Histograma_CosDir_Y_Reg3.txt')


#Abre histogramas de cosenos directores guardados en carpeta 

# HistogramaRegTotal_X = AbreHistogramaTXT('C:\Users\Install\Documents\EGSnrc\Phase Space\VarianClinaciX_6MV_20x20_w1\Cosenos directores por region/Histograma_CosDir_X_Reg_Total.txt')
# HistogramaReg1_X = AbreHistogramaTXT('C:\Users\Install\Documents\EGSnrc\Phase Space\VarianClinaciX_6MV_20x20_w1\Cosenos directores por region/Histograma_CosDir_X_Reg1.txt')
# HistogramaReg2_X = AbreHistogramaTXT('C:\Users\Install\Documents\EGSnrc\Phase Space\VarianClinaciX_6MV_20x20_w1\Cosenos directores por region/Histograma_CosDir_X_Reg2.txt')
# HistogramaReg3_X = AbreHistogramaTXT('C:\Users\Install\Documents\EGSnrc\Phase Space\VarianClinaciX_6MV_20x20_w1\Cosenos directores por region/Histograma_CosDir_X_Reg3.txt')

# HistogramaRegTotal_Y = AbreHistogramaTXT('C:\Users\Install\Documents\EGSnrc\Phase Space\VarianClinaciX_6MV_20x20_w1\Cosenos directores por region/Histograma_CosDir_Y_Reg_Total.txt')
# HistogramaReg1_Y = AbreHistogramaTXT('C:\Users\Install\Documents\EGSnrc\Phase Space\VarianClinaciX_6MV_20x20_w1\Cosenos directores por region/Histograma_CosDir_Y_Reg1.txt')
# HistogramaReg2_Y = AbreHistogramaTXT('C:\Users\Install\Documents\EGSnrc\Phase Space\VarianClinaciX_6MV_20x20_w1\Cosenos directores por region/Histograma_CosDir_Y_Reg2.txt')
# HistogramaReg3_Y = AbreHistogramaTXT('C:\Users\Install\Documents\EGSnrc\Phase Space\VarianClinaciX_6MV_20x20_w1\Cosenos directores por region/Histograma_CosDir_Y_Reg3.txt')

#Hace histograma de cosenos directores en X e Y para diferentes energias a partir de los archivos txt

# GraficaHistogramaTotalCosenos_por_region(HistogramaReg3_Y, HistogramaRegTotal_Y, 40)		#Cambiar el titulo del grafico a Coseno Director segun sea en X o Y y el nombre de la region
     

#-------------------------------------------------------------------------- 25/11 ------------------------------------------------------------------

path = 'C:\Users\Install\Documents\EGSnrc\Phase Space/'

filename_homog_cos_x = 'pruebacosdir_homogeneo_x.txt'

filename_homog_cos_y = 'pruebacosdir_homogeneo_y.txt'

filename_homog_ang_x = 'pruebaangdir_homogeneo_x.txt'

filename_homog_ang_y = 'pruebaangdir_homogeneo_y.txt'


filename_no_homog_ang_x = 'pruebaangdir_no_homogeneo_x.txt'

filename_no_homog_cos_x = 'pruebacosdir_no_homogeneo_x.txt'

findFiles(path)

bins = 100

# Histograma = CalculaHistogramaTotalCosenos_por_region(path, 6.0, bins)

# GuardaHistogramaTotalCosenos(Histograma, path + filename_homog_cos_x)


Histograma_homogeneo_ang = AbreHistogramaTXT(path + filename_homog_ang_x)

Histograma_homogeneo_cos = AbreHistogramaTXT(path + filename_homog_cos_x)

GraficaHistogramaTotalCosenos_por_region(Histograma_homogeneo_ang, Histograma_homogeneo_cos, bins)
