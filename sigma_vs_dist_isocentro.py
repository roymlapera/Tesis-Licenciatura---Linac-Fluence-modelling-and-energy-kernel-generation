import numpy 		as np
from glob           import glob 
from os.path        import splitext, split, dirname
import matplotlib.pyplot as plt
import pylab                
from scipy.optimize import curve_fit


def findFiles(phsp_path):                                                                           # FINDS FILES INSIDE FOLDER WITH IAEAphsp EXTENSION
    find_names  = phsp_path + '/*.IAEAphsp'                                                       # Names are assumed to be .......... .IAEAheader
    all_names   = glob(find_names)                                                                  # Found files with the searched names
    
    phsp_files  = []                                                                                # List to put correct names
    for name in all_names: phsp_files.append(splitext(name)[0])                                     # Appends name without extension to files
    
    print 'Archivos encontrados:'

    print phsp_files

    # Espacio de fase (chico: 4GB) con colimador secundario --> 0

    # Espacio de fase (30GB) sin colimador secundario ---> 1

    # Espacio de fase (30GB) sin colimador secundario (ESTEPE de 0.20)---> 1

    phsp_files = np.delete(phsp_files,[0,2])

    print 'Abriendo:'

    print phsp_files

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

#Calcula el espectro de cosenos directores u y v de un tipo de particulas EN UNA DETERMINADA REGION

def CalculaHistogramaCosenos_Pack_por_region(data, type, maxEnergy, bins, cos_dir, x_inf, x_sup, y_inf, y_sup):            
    histo       = np.zeros(bins)                                                                       # Espectro de coseno director de ELECTRONES
    
    energy      = np.abs(data[:, 0])
    PosX        = data[:, 1] 
    PosY        = data[:, 2]
    weight      = data[:, 5]                                                                             # Weight of every photon in Phase plane of current unpack

    if(cos_dir == 'x'):
        xCos = data[:,3]
        isElectron    = 1*(type == 2) * (energy < maxEnergy) * (x_inf<PosX) * (PosX<x_sup) * (y_inf<PosY) * (PosY<y_sup)                                    
        cosInd      = isElectron * ((xCos+1.0) * bins /2.0).astype(int)                             # (el 2 es porque el coseno director puede variar de -1 a 1)
        sumation    = isElectron * weight
    elif(cos_dir == 'y'):
        yCos        = data[:, 4]                                                                    
        isElectron    = 1*(type == 2) * (energy < maxEnergy) * (PosX<x_sup) * (y_inf<PosY) * (PosY<y_sup)  
        cosInd      = isElectron * ((yCos + 1.0) * bins / 2.0).astype(int)                          #  (el 2 es porque el coseno director puede variar de -1 a 1)                                 
        sumation    = isElectron * weight
    else:
        print 'Error en el tipo de coseno director'

    for ind, value in enumerate(cosInd): 
        histo[value] += sumation[ind]                                                               # Spectrum photon countingt 

    return histo

def CalculaHistogramaCosenos_por_region(data_path, maxEnergy, bins, cos_dir, x_inf, x_sup, y_inf, y_sup):                      
    histo    = np.zeros(bins)                                                                       # Energy spectrum (E)
    phsp_files   = findFiles(data_path)                                                             # Finds IAEA files to work with
    for phsp in phsp_files:
        particles, packFormat = parseFile(phsp)                                                     # Specifies format and amount of particles to read
        
        with open(phsp + '.IAEAphsp', 'rb') as IAEAphsp:                                            # Opens binary file to read data
            for x in range(particles/100000):
                pack        = np.fromfile(IAEAphsp, dtype=packFormat, count=100000)                   # COUNT = NUMBER OF ITEMS TO READ!!!
                histo   += CalculaHistogramaCosenos_Pack_por_region(pack['Data'], pack['Type'], maxEnergy, bins, cos_dir, x_inf, x_sup, y_inf, y_sup)          

    return histo                                                                        # Spectrum tiene los valores del eje Y del histograma en funcion del numero de bin 


def GuardaHistogramaTotalCosenos(data, file): np.savetxt(file, data)

def AbreHistogramaTXT(file):
    with open(file) as f:
        histo = np.loadtxt(f, delimiter="\n", dtype='double', comments="#", skiprows=0, usecols=None)
    return histo


def GraficaHistogramaTotalCosenos_por_region(data, path, bins, cos_dir, norm):                                                                # CAMBIAR ACA LA LEYENDA SEGUN X E Y

    opacity = 0.8

    abcisas = np.arange(bins)

    abcisas = abcisas *2.0 /float(bins-1) - 1.0

    #Normalizacion:

    if(norm == 1): data /= data.sum()
    elif(norm != 0): print 'Valor no esperado de norm. Valor esperado: "1" o "0".'

    abcisas = np.ravel(zip(abcisas,abcisas + 2.0 /float(bins-1)))

    abcisas = abcisas -1.0 /float(bins-1)

    data = np.ravel(zip(data,data))


    #fig, ax = plt.subplots()


    rects = plt.plot(abcisas, data, alpha=opacity, linewidth = 2)

    # plt.axis([-1.1, 1.1, 0, 0.15])

    plt.grid(True)


    if(cos_dir == 'x'):
        plt.xlabel('Coseno director en X', fontsize=14)
        #fig.savefig(path + 'Coseno_dir_en_x.png', dpi = 200, bbox_inches='tight')
    elif(cos_dir == 'y'):
        plt.xlabel('Coseno director en Y', fontsize=14)
        #fig.savefig(path + 'Coseno_dir_en_y_reg3.png', dpi = 200, bbox_inches='tight')
    else:
        print 'Error en el tipo de coseno director'


#Calculo de sigma y media por region

def Devuelve_Coseno_dir_por_region(data, type, maxEnergy, cos_dir, x_inf, x_sup, y_inf, y_sup):            
    
    energy      = np.abs(data[:, 0])
    PosX        = data[:, 1] 
    PosY        = data[:, 2]
    weight      = data[:, 5]                                                                             # Weight of every photon in Phase plane of current unpack

    # print type(cos_dir)


    if cos_dir is 'x':
        print 'Hola entre'
        xCos                = data[:,3]
        isElectron          = 1*(type == 2) * (energy < maxEnergy) * (x_inf<PosX) * (PosX<x_sup) * (y_inf<PosY) * (PosY<y_sup)                                    
        cos_director        = xCos * isElectron 
    elif cos_dir is 'y':
        print 'Hola entre aca'
        yCos                = data[:, 4]                                                                    
        isElectron          = 1*(type == 2) * (energy < maxEnergy) * (x_inf<PosX) * (PosX<x_sup) * (y_inf<PosY) * (PosY<y_sup)                                    
        cos_director        = yCos * isElectron 
    else:
        print 'Error en el tipo de coseno director!!!!!!!!!!!!!!!!!'


    weight *= isElectron

    #asumo que precindir de las particulas cuyo coseno director es 0, no afecta al calculo

    cos_director = np.delete(cos_director,np.where(cos_director == 0)[0])
    weight = np.delete(weight,np.where(cos_director == 0)[0])

    return cos_director, weight

def Calcula_prom_var_por_region(data_path, maxEnergy, cos_dir, x_inf, x_sup, y_inf, y_sup):                      
    phsp_files   = findFiles(data_path)                                                          
    for phsp in phsp_files:
        particles, packFormat = parseFile(phsp)                                                  
        
        with open(phsp + '.IAEAphsp', 'rb') as IAEAphsp:                                         
            for x in range(1):
                if (x==0):
                    pack                  = np.fromfile(IAEAphsp, dtype=packFormat, count=100000)                
                    cos_director, weight       = Devuelve_Coseno_dir_por_region(pack['Data'], pack['Type'], maxEnergy, cos_dir, x_inf, x_sup, y_inf, y_sup)           
                    pack                      = np.fromfile(IAEAphsp, dtype=packFormat, count=100000)               

                cos_director_pack, weight_pack = Devuelve_Coseno_dir_por_region(pack['Data'], pack['Type'], maxEnergy, cos_dir, x_inf, x_sup, y_inf, y_sup)
                cos_director = np.concatenate( (cos_director,cos_director_pack) )  
                weight = np.concatenate( (weight,weight_pack) )            
    
    print cos_director.shape
    print weight.shape
    promedio = np.average(cos_director, weights=weight)

    varianza = np.average((cos_director-promedio)**2, weights=weight)

    return promedio,varianza 

##############################################  REGIONES DISPUESTAS HORIZONTALMENTE  ##########################################

#Genero la data de histogramas en zona horizontal y las guardo, las regiones son de 4cm x 4cm

def Genera_histos_txt_horizontal(phsp_path, folder_path, maxEnergy,bins):
    for i in range(8):
        x_inf = -30.0 + 4.0*i 
        x_sup = -30.0 + 4.0 + 4.0*i
        y_inf = -4.0
        y_sup = 4.0

        cos_dir = 'x'
        Histograma_x = CalculaHistogramaCosenos_por_region(phsp_path, maxEnergy, bins, cos_dir, x_inf, x_sup, y_inf, y_sup)
        GuardaHistogramaTotalCosenos(Histograma_x, folder_path + 'Histograma_horiz_' + cos_dir + '_' + str(i) + '.txt')

        cos_dir = 'y'
        Histograma_y = CalculaHistogramaCosenos_por_region(phsp_path, maxEnergy, bins, cos_dir, x_inf, x_sup, y_inf, y_sup)
        GuardaHistogramaTotalCosenos(Histograma_y, folder_path + 'Histograma_horiz_' + cos_dir + '_' + str(i) + '.txt')

#Grafico los histogramas guardados en txt

def Grafica_histo_en_txt_horizontal(phsp_path, folder_path, maxEnergy,bins):
    for i in range(8):
        # cos_dir = 'x'

        # histo = AbreHistogramaTXT(folder_path + 'Histograma_horiz_' + cos_dir + '_' + str(i) + '.txt')
        # GraficaHistogramaTotalCosenos_por_region(histo, folder_path, bins, cos_dir, 1)

        cos_dir = 'y'

        histo = AbreHistogramaTXT(folder_path + 'Histograma_horiz_' + cos_dir + '_' + str(i) + '.txt')
        GraficaHistogramaTotalCosenos_por_region(histo, folder_path, bins, cos_dir, 1)
                                                                                                                                                                                    
    plt.show()  

#Calcula el promedio y la varianza de particulas por region    SIN BINES

def Calcula_sigma_vs_dist_iso_horizontal(phsp_path, folder_path, maxEnergy):
    dist_isocentro = np.zeros(8)

    prom_x = np.zeros(8)
    var_x = np.zeros(8)

    prom_y = np.zeros(8)
    var_y = np.zeros(8)


    for i in range(8):
        x_inf = -30.0 + 4.0*i 
        x_sup = -30.0 + 4.0 + 4.0*i
        y_inf = -4.0
        y_sup = 4.0

        dist_isocentro[i] = 28.0 - 4.0 * i
        cos_dir = 'x'
        # print cos_dir
        prom_x[i], var_x[i] = Calcula_prom_var_por_region(phsp_path, maxEnergy, cos_dir, x_inf, x_sup, y_inf, y_sup)

        cos_dir = 'y'
        prom_y[i], var_y[i] = Calcula_prom_var_por_region(phsp_path, maxEnergy, cos_dir, x_inf, x_sup, y_inf, y_sup)

    print prom_x, var_x, prom_y, var_y

    np.savetxt(folder_path + 'datos_considerando_peso.txt', np.c_[prom_x,var_x,prom_y,var_y])

def Grafica_varianza_vs_dist_iso_horizontal(phsp_path, folder_path, maxEnergy,data):
    varianza_x = data[:,1].astype(np.float64) 
    varianza_y = data[:,3].astype(np.float64) 

    # print varianza_x,varianza_y

    dist_isocentro = np.zeros(8)

    # for i in range(8):
    #     dist_isocentro[i] = 28.0 - 4.0 * i


    fig = plt.figure()

    plt.plot(dist_isocentro, varianza_x,'b', label = 'Coseno director en x')

    plt.plot(dist_isocentro, varianza_y,'g', label = 'Coseno director en y')

    plt.xlabel('Distancia al isocentro [cm]', fontsize=14)
    plt.ylabel('Varianza del coseno director', fontsize=14)
    plt.axis([0, 30, -1, 1])

    plt.legend(loc=2,prop={'size':14})

    plt.show()

def Grafica_promedio_vs_dist_iso_horizontal(phsp_path, folder_path, maxEnergy,data):
    promedio_x = data[:,0].astype(np.float64) 
    promedio_y = data[:,2].astype(np.float64) 

    # print varianza_x,varianza_y

    dist_isocentro = np.zeros(8)

    for i in range(8):
        dist_isocentro[i] = 28.0 - 4.0 * i

    fig = plt.figure()

    plt.plot(dist_isocentro, promedio_x,'b*', label = 'Coseno director en x')

    plt.plot(dist_isocentro, promedio_y,'g*', label = 'Coseno director en y')

    plt.xlabel('Distancia al isocentro [cm]', fontsize=14)
    plt.ylabel('Promedio del coseno director', fontsize=14)
    plt.axis([0, 30, -0.5, 0.5])

    plt.legend(loc=2,prop={'size':14})

    plt.show()

##############################################  REGIONES DISPUESTAS EN DIAGONAL  ##########################################

#Genero la data de histogramas en zona diagonal y las guardo, las regiones son de 4cm x 4cm HAY QUE ARREGLAR ESTO; LA VARIANZA ESTA MAL

def Genera_histos_txt_diagonal(phsp_path, folder_path, maxEnergy,bins):
    for i in range(8):
        x_inf = -30.0 + 4.0*i 
        x_sup = -30.0 + 4.0 + 4.0*i
        y_inf = 30.0 - 4.0 - 4.0*i
        y_sup = 30.0 - 4.0*i

        cos_dir = 'x'
        Histograma_x = CalculaHistogramaCosenos_por_region(phsp_path, maxEnergy, bins, cos_dir, x_inf, x_sup, y_inf, y_sup)
        GuardaHistogramaTotalCosenos(Histograma_x, folder_path + 'Histograma_horiz_' + cos_dir + '_' + str(i) + '.txt')

        cos_dir = 'y'
        Histograma_y = CalculaHistogramaCosenos_por_region(phsp_path, maxEnergy, bins, cos_dir, x_inf, x_sup, y_inf, y_sup)
        GuardaHistogramaTotalCosenos(Histograma_y, folder_path + 'Histograma_horiz_' + cos_dir + '_' + str(i) + '.txt')

def Grafica_histo_en_txt_diagonal(phsp_path, folder_path, maxEnergy,bins):
    for i in range(8):
        cos_dir = 'x'

        histo = AbreHistogramaTXT(folder_path + 'Histograma_horiz_' + cos_dir + '_' + str(i) + '.txt')
        GraficaHistogramaTotalCosenos_por_region(histo, folder_path, bins, cos_dir, 1)

        # cos_dir = 'y'

        # histo = AbreHistogramaTXT(folder_path + 'Histograma_horiz_' + cos_dir + '_' + str(i) + '.txt')
        # GraficaHistogramaTotalCosenos_por_region(histo, folder_path, bins, cos_dir, 1)
                                                                                                                                                                                    
    plt.show() 

def Calcula_sigma_vs_dist_iso_diagonal(phsp_path, folder_path, maxEnergy,data):
    dist_isocentro = np.zeros(8)

    prom_x = np.zeros(8)
    var_x = np.zeros(8)

    prom_y = np.zeros(8)
    var_y = np.zeros(8)


    for i in range(8):
        x_inf = -30.0 + 4.0*i 
        x_sup = -30.0 + 4.0 + 4.0*i
        y_inf = 30.0 - 4.0 - 4.0*i
        y_sup = 30.0 - 4.0*i

        dist_isocentro[i] = 2.0*np.sqrt(2.0) * (7.0-i)

        cos_dir = 'x'
        prom_x[i], var_x[i] = Calcula_prom_var_por_region(phsp_path, maxEnergy, cos_dir, x_inf, x_sup, y_inf, y_sup)

        cos_dir = 'y'
        prom_y[i], var_y[i] = Calcula_prom_var_por_region(phsp_path, maxEnergy, cos_dir, x_inf, x_sup, y_inf, y_sup)

    print prom_x, var_x, prom_y, var_y

    np.savetxt(folder_path + 'datos_s_bins.txt', np.c_[prom_x,var_x,prom_y,var_y])

#-------------------------------------------------------------------- MAIN -----------------------------------------------------------------------------


phsp_path = 'C:\Users\Install\Documents\EGSnrc\Phase Space/'

folder_path = 'C:\Users\Install\Desktop\Tesis Roy\sigma_vs_dist_isocentro\sColSec\Horizontal/'

maxEnergy = 6.0

bins = 100

Calcula_sigma_vs_dist_iso_horizontal(phsp_path, folder_path, maxEnergy)

# print range(1000000000000/100000)

# data = np.loadtxt('C:\Users\Install\Desktop\Tesis Roy\sigma_vs_dist_isocentro\cColSec\Horizontal/datos_considerando_peso.txt')

# Grafica_promedio_vs_dist_iso_horizontal(phsp_path, folder_path, maxEnergy,data)

# Grafica_histo_en_txt_horizontal(phsp_path, folder_path, maxEnergy,bins)










