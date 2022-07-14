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

# Calcula espectro de energia por intervalo de energia:

def calculatePack(data, type, maxEnergy, bins):                                                     # CALCULATES FULL SPECTRUM FROM PACK DATA
    spectrum    = np.zeros(bins)                                                                       # Energy spectrum (E)
        
    #radial      = (data[:, 1]**2 + data[:, 2]**2)**.5
    energy      = np.abs(data[:, 0])
    zCos        = (1 - (data[:, 3])**2 - (data[:, 4])**2)**.5                                       # Cosin angle in Z direction of every photon in Phase plane of current unpack
    weight      = data[:, 5]                                                                        # Weight                     of every photon in Phase plane of current unpack
    
    isPhoton    = 1*(type == 1) * (energy < maxEnergy)					 							# Only photons inside calculation limit are proccessed
    eneInd      = isPhoton * (energy * bins / maxEnergy).astype(int)                                # Index for spectrum sumation
    sumation    = isPhoton * energy * weight / zCos

    #isElectron = 1*(type == ) * (energy < maxEnergy)
    #eneInd      = isElectron * (energy * bins / maxEnergy).astype(int)                                # Index for spectrum sumation
    #sumation    = isElectron * energy * weight / zCos                                                  # Sumation value for spectrum

    for ind, value in enumerate(eneInd): 
    	spectrum[value] += sumation[ind]                           # Spectrum photon countingt 
                                                                   
    return spectrum

def kernelWeight(data_path, maxEnergy, bins=32):                                    # OPENS IAEA FILE AND CALCULATES PARTICLE SPECTRUM  IN small DISC
    spectrum     = np.zeros(bins)                                                                      # Energy spectrum (E)
    phsp_files   = findFiles(data_path)                                                  # Finds IAEA files to work with

    for phsp in phsp_files:
        particles, packFormat = parseFile(phsp)                                                     # Specifies format and amount of particles to read
        
        with open(phsp + '.IAEAphsp', 'rb') as IAEAphsp:                                            # Opens binary file to read data
            for x in range(1):
                pack        = np.fromfile(IAEAphsp, dtype=packFormat, count=10000000)                   # Gets full pack data from one file
                spectrum   += calculatePack(pack['Data'], pack['Type'], maxEnergy, bins)            # Calculates spectrum part from one file

    spectrum /= spectrum.sum()

    print spectrum


#Calcula el numero de un tipo de particula en el espaciod e fases: 

def calculasuma(data, type, maxEnergy):                                                     # SUMA PARTICULAS EN EL PACK DAT
    suma = 0
    energy      = np.abs(data[:, 0])                                       
    weight      = data[:, 5]                                                                       
    
    #isPhoton    = 1*(type == 1) * (energy < maxEnergy)					 							# Only photons inside calculation limit are proccessed

    isElectron = 1*(type == 2) * (energy < maxEnergy) 

    #isPositron = 1*(type == 3) * (energy < maxEnergy)                                                         												   
                                                                    
    return isElectron.sum()

def sumaparticulas (data_path, maxEnergy):                                                                    
	phsp_files   = findFiles(data_path)
	suma = 0
	for phsp in phsp_files:
        
		particles, packFormat = parseFile(phsp)                                                     # Specifies format and amount of particles to read
        
        with open(phsp + '.IAEAphsp', 'rb') as IAEAphsp:                                            # Opens binary file to read data
            for x in range(particles):
                pack  = np.fromfile(IAEAphsp, dtype=packFormat, count=-1)            		   # Gets full pack data from one file
                suma  += calculasuma(pack['Data'], pack['Type'], maxEnergy, bins)            		# Suma las constribuciones de particulas del pack del archivo 

	print suma


#Calcula el espectro de energia de los electrones

def calculatePackEnergiaCinetica(data, type, maxEnergy, bins):                                                     # CALCULATES FULL SPECTRUM FROM PACK DATA
    spectrum    = np.zeros(bins)                                                                       # Energy spectrum (E)
        
    #radial      = (data[:, 1]**2 + data[:, 2]**2)**.5
    energy      = np.abs(data[:, 0])
    #zCos        = (1 - (data[:, 3])**2 - (data[:, 4])**2)**.5                                       # Cosin angle in Z direction of every photon in Phase plane of current unpack
    weight      = data[:, 5]                                                                        # Weight                     of every photon in Phase plane of current unpack
    
    isElectron    = 1*(type == 2) * (energy < maxEnergy)					 							# Only photons inside calculation limit are proccessed
    eneInd      = isElectron * (energy * bins / maxEnergy).astype(int)                                # Index for spectrum sumation
    sumation    = isElectron * weight

    # print energy, weight

    for ind, value in enumerate(eneInd): 
    	spectrum[value] += sumation[ind]                           # Spectrum photon countingt 
                                                                  
    return spectrum

def CalculaEspectroEnergiaCinetica(data_path, maxEnergy, bins):                                    # OPENS IAEA FILE AND CALCULATES PARTICLE SPECTRUM  IN small DISC
    spectrum     = np.zeros(bins)                                                                      # Energy spectrum (E)
    phsp_files   = findFiles(data_path)                                                  # Finds IAEA files to work with

    for phsp in phsp_files:
        particles, packFormat = parseFile(phsp)                                                     # Specifies format and amount of particles to read
        
        with open(phsp + '.IAEAphsp', 'rb') as IAEAphsp:                                            # Opens binary file to read data
            for x in range(particles/100000):
                pack        = np.fromfile(IAEAphsp, dtype=packFormat, count=100000)                   # Gets full pack data from one file
                spectrum   += calculatePackEnergiaCinetica(pack['Data'], pack['Type'], maxEnergy, bins)            # Calculates spectrum part from one file

    # spectrum /= spectrum.sum()

    return spectrum

def GraficaHistogramaEnergiaCinetica(data, maxEnergy, bins, path):									
    
    masa_en_reposo_electron = 0.511  #MeV

    #normalizacion
    data /= ( data.sum() * (maxEnergy - masa_en_reposo_electron)/bins )

    abcisas = np.arange(bins)

    abcisas = abcisas /float(bins-1) * maxEnergy - masa_en_reposo_electron

    fig = plt.figure() 

    abcisas = np.ravel(zip(abcisas,abcisas+maxEnergy/bins))
    data= np.ravel(zip(data,data))

    plt.plot(abcisas, data, "r")

    plt.axis([0, 6.1, 0.0, 0.81])

    plt.xlabel('Energia Cinetica[MeV]', fontsize=14)

    plt.grid(True)

    plt.show()

    fig.savefig(path + 'espectro.png', dpi = 200, bbox_inches='tight')

    return abcisas

def GuardaEspectroEnergiaCinetica(archivo,data): np.savetxt(archivo, data)

def AbreEspectroTXT(file):
    with open(file) as f:
        spectrum = np.loadtxt(f, delimiter="\n", dtype='double', comments="#", skiprows=0, usecols=None)
    return spectrum

# -------------------------------------------------------------------------- MAIN ------------------------------------------------------------------

#Se ejecutan funciones para calcular espectro de energia por intervalo de energia (contribucion de energia en el eje Z en funcion de la energia de la particula)

#kernelWeight('C:\Users\Install\Documents\EGSnrc\Phase Space/', 6.0, 32)


#Cuenta la cantidad de particulas de cierto tipo (e-, p+, foton)

#sumaparticulas('C:\Users\Install\Documents\EGSnrc\Phase Space/', 6.0)


#Hace un espectro de energia de los e-

path = 'C:\Users\Install\Documents\EGSnrc\Phase Space/'

filename =  'Espectro_w1.txt'

Espectro = CalculaEspectroEnergiaCinetica(path, 6.0, 100)

GuardaEspectroEnergiaCinetica(path + filename,Espectro)

# Espectro = AbreEspectroTXT(path + filename)

absisas = GraficaHistogramaEnergiaCinetica(Espectro,6.0,100,path)



np.savez('espectrum_patronus', energia = absisas, pesos_espectro = Espectro)