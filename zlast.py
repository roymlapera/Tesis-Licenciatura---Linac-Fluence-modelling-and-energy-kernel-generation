import numpy 		as np
from glob           import glob 
from os.path        import splitext, split, dirname
import matplotlib.pyplot as plt
from mpl_toolkits.axes_grid1.inset_locator import inset_axes
from mpl_toolkits.axes_grid1.inset_locator import mark_inset


def findFiles(phsp_path):                                                                           # FINDS FILES INSIDE FOLDER WITH IAEAphsp EXTENSION
    find_names  = phsp_path + '/*.IAEAphsp'                                                       # Names are assumed to be .......... .IAEAheader
    all_names   = glob(find_names)                                                                  # Found files with the searched names
    
    phsp_files  = []                                                                                # List to put correct names
    for name in all_names: phsp_files.append(splitext(name)[0])                                     # Appends name without extension to files
    
    print 'Espacios de fases encontrados:'

    print '\n'.join(phsp_files)

    phsp_files = np.delete(phsp_files,[0,2,3,4,5])

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

#Funciones que calculan histograma de Zlast de electrones

def CalculaHistogramaZlast_Pack(data, type, maxEnergy, bins):                                          # CALCULATES FULL SPECTRUM FROM PACK DATA
    histo    = np.zeros(bins)                                                                       # Espectro de coseno director de ELECTRONES
    
    maxZ = 100.0

    energy      = np.abs(data[:, 0])
    weight      = data[:, 5]                                                                             # Weight of every photon in Phase plane of current unpack
    Zlast       = data[:, 6] 
                                                                
    isElectron    = 1*(type == 2) * (energy < maxEnergy) * (Zlast < maxZ)
    zInd      = isElectron * (Zlast * bins / maxZ).astype(int)                           # Index for spectrum sumation (el 2 es porque el coseno director puede variar de -1 a 1)                                 
    sumation    = isElectron * weight

    for ind, value in enumerate(zInd): 
        histo[value] += sumation[ind]                                                               # Spectrum photon countingt 

    return histo

def CalculaHistogramaTotalZlast(data_path, maxEnergy, bins):                                    # OPENS IAEA FILE AND CALCULATES PARTICLE SPECTRUM  IN small DISC
    histo     = np.zeros(bins)                                                                      # Energy spectrum (E)
    phsp_files   = findFiles(data_path)                                                  # Finds IAEA files to work with
    for phsp in phsp_files:
        particles, packFormat = parseFile(phsp)                                                     # Specifies format and amount of particles to read
        
        with open(phsp + '.IAEAphsp', 'rb') as IAEAphsp:                                            # Opens binary file to read data
            for x in range(particles/100000):
                pack        = np.fromfile(IAEAphsp, dtype=packFormat, count=100000)                   # COUNT = NUMBER OF ITEMS TO READ!!!
                histo       += CalculaHistogramaZlast_Pack(pack['Data'], pack['Type'], maxEnergy, bins)            # Calculates spectrum part from one file

    histo /= histo.sum()

    return histo                                                                         # Spectrum tiene los valores del eje Y del histograma en funcion del numero de bin 

def GuardaHistogramaTotalZlast(data, file):
    np.savetxt(file, data)

def AbreHistogramaTXT(file):
    with open(file) as f:
        histo = np.loadtxt(f, delimiter="\n", dtype='double', comments="#", skiprows=0, usecols=None)
    return histo

def GraficaHistogramaTotalZlast(data, bins, path):                                                              # CAMBIAR ACA LA LEYENDA SEGUN X E Y

    maxZ = 100.0
    opacity = 0.5

    abcisas = np.arange(bins)
    abcisas = abcisas * maxZ /float(bins-1)

    #normalizacion
    data /= data.sum() 

    fig, ax = plt.subplots() # create a new figure with a default 111 subplot

    abcisas = np.ravel(zip(abcisas,abcisas+100.0/bins))
    data= np.ravel(zip(data,data))

    ax.plot(abcisas, data, "g")

    plt.xticks( np.arange(min(abcisas), max(abcisas)+1.0, 5.0) )

    plt.axis([-1.0, 101.0, 0.0, 0.24])

    plt.xlabel('Z last [cm]', fontsize=14)

    plt.grid(True)


    # axins = zoomed_inset_axes(ax, 3, loc=1) # zoom-factor: 2.5, location: upper-left
    axins = inset_axes(ax, 4,2.5 , loc=1,bbox_to_anchor=(0.88, 0.88),bbox_transform=ax.figure.transFigure) # no zoom
    axins.plot(abcisas, data, "g")

    x1, x2, y1, y2 = 27.0, 100.5, 0, 0.006 # specify the limits
    axins.set_xlim(x1, x2) # apply the x-limits
    axins.set_ylim(y1, y2) # apply the y-limits
    plt.tick_params(labelsize=10)
    plt.grid(True)

    mark_inset(ax, axins, loc1=2, loc2=1, fc="none", ec="0.5")

    fig.savefig(path + 'zlast.png', dpi = 200, bbox_inches='tight')

    plt.show()

# -------------------------------------------------------------------------------- MAIN ----------------------------------------------------------------------------------

path = 'C:\Users\Install\Documents\EGSnrc\Phase Space/'

#Guardamos el histograma en un txt:

# HistogramaZlast = CalculaHistogramaTotalZlast(path, 6.0, 500)

# GuardaHistogramaTotalZlast(HistogramaZlast,path + 'HistogramaZlast_w1.txt')


# #Abrimos la informacion y la usamos

HistogramaZlast = AbreHistogramaTXT(path + 'HistogramaZlast_w1.txt')

GraficaHistogramaTotalZlast(HistogramaZlast, 500, path)

# findFiles('C:\Users\Install\Documents\EGSnrc\Phase Space/')