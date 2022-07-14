import numpy 		as np
from glob           import glob 
from os.path        import splitext, split, dirname
import matplotlib.pyplot as plt
from scipy.optimize import curve_fit

def gauss(x, *p):
    A, mu, sigma = p
    return A*np.exp(-(x-mu)**2/(2.*sigma**2))

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

##################################################################################################################################################################

def DevuelveCosenosOAngulos(data, type, maxEnergy, bins, variable, axis):                                          # CALCULATES FULL SPECTRUM FROM PACK DATA
    energy      = np.abs(data[:, 0])
    PosX        = data[:, 1] 
    PosY        = data[:, 2]
    weight      = data[:, 5] 																			 # Weight of every photon in Phase plane of current unpack

    if (variable == 'cos' and axis == 'x'):
	    xCos = data[:,3]
	    isElectron    = 1*(type == 2) * (energy < maxEnergy) 							
	    weight = weight*isElectron
	    xCos = xCos[weight>0]
	    weight = weight[weight>0]
	    return xCos, weight

    if (variable == 'cos' and axis == 'y'):
	    yCos = data[:,4]
	    isElectron    = 1*(type == 2) * (energy < maxEnergy) 							
	    weight = weight*isElectron
	    yCos = yCos[weight>0]
	    weight = weight[weight>0]
	    return yCos, weight

    if (variable == 'ang' and axis == 'x'):
    	xCos = data[:,3]
    	isElectron    = 1*(type == 2) * (energy < maxEnergy) 
    	ang = np.arccos(xCos)							
    	weight = weight*isElectron
    	ang = ang[weight>0]
    	weight = weight[weight>0]
    	return ang, weight

    if (variable == 'ang' and axis == 'y'):
	    yCos = data[:,4]
	    isElectron    = 1*(type == 2) * (energy < maxEnergy) 
	    ang = np.arccos(yCos)							
	    weight = weight*isElectron
	    ang = ang[weight>0]
	    weight = weight[weight>0]
	    return ang, weight

def DevuelveCosenosOAngulosDeArchivo(data_path, maxEnergy, bins, variable, axis):                                    # OPENS IAEA FILE AND CALCULATES PARTICLE SPECTRUM  IN small DISC
    phsp_files   = findFiles(data_path)                                                  # Finds IAEA files to work with
    for phsp in phsp_files:
        particles, packFormat = parseFile(phsp)                                                     # Specifies format and amount of particles to read
        
        with open(phsp + '.IAEAphsp', 'rb') as IAEAphsp:                                            # Opens binary file to read data
            pack       = np.fromfile(IAEAphsp, dtype=packFormat, count=-1)                   # COUNT = NUMBER OF ITEMS TO READ!!!
            data, weight   = DevuelveCosenosOAngulos(pack['Data'], pack['Type'], maxEnergy, bins, variable, axis)            # Calculates spectrum part from one file
    return data, weight  	

def GraficaHistogramaCosenosAngulos(data1, data2, weight1, weight2, bins):																# CAMBIAR ACA LA LEYENDA SEGUN X E Y

    path = 'C:\Users\Install\Desktop\Tesis Roy\cColSec\Cosenos directores de electrones\Por region/'

    fig, ax = plt.subplots()

    plt.hist(data1, bins=bins, histtype='step', normed=True, color='b', label='Coseno director respecto al eje X', weights=weight1)
 
    plt.hist(data2, bins=bins, histtype='step', normed=True, color='g', label='Coseno director respecto al eje Y', weights=weight2)

    plt.legend(fontsize=11)

    plt.axis([-1.1, 1.1, 0, 1.7])

    plt.xlabel('Coseno director', fontsize=14)

    plt.grid(True)

    fig.savefig(path + 'Coseno_dir_globales_w1.png', dpi = 200, bbox_inches='tight')

    plt.show()

##################################################################################################################################################################

def DevuelveCosenosOAngulos_PorRegion(data, type, maxEnergy, bins, variable, axis, x_inf, x_sup, y_inf, y_sup): 
    energy      = np.abs(data[:, 0])
    PosX        = data[:, 1] 
    PosY        = data[:, 2]
    weight      = data[:, 5]                                                                             # Weight of every photon in Phase plane of current unpack

    if (variable == 'cos' and axis == 'x'):
        xCos = data[:,3]
        isElectron    = 1*(type == 2) * (energy < maxEnergy)  * (x_inf<PosX) * (PosX<x_sup) * (y_inf<PosY) * (PosY<y_sup)                            
        weight = weight*isElectron
        xCos = xCos[weight>0]
        weight = weight[weight>0]
        return xCos, weight

    if (variable == 'cos' and axis == 'y'):
        yCos = data[:,4]
        isElectron    = 1*(type == 2) * (energy < maxEnergy)  * (x_inf<PosX) * (PosX<x_sup) * (y_inf<PosY) * (PosY<y_sup)                           
        weight = weight*isElectron
        yCos = yCos[weight>0]
        weight = weight[weight>0]
        return yCos, weight

def DevuelveCosenosOAngulosDeArchivo_PorRegion(data_path, maxEnergy, bins, variable, axis, x_inf, x_sup, y_inf, y_sup):   
    phsp_files   = findFiles(data_path)                                
    for phsp in phsp_files:
        particles, packFormat = parseFile(phsp)                                                   
        
        with open(phsp + '.IAEAphsp', 'rb') as IAEAphsp:                                        
            pack       = np.fromfile(IAEAphsp, dtype=packFormat, count=-1)             
            data, weight   = DevuelveCosenosOAngulos_PorRegion(pack['Data'], pack['Type'], maxEnergy, bins, variable, axis, x_inf, x_sup, y_inf, y_sup)
    return data, weight     

def GraficaHistogramaCosenosAngulos_PorRegion(data1, data2, weight1, weight2, bins):                                                              # CAMBIAR ACA LA LEYENDA SEGUN X E Y

    path = 'C:\Users\Install\Desktop\Tesis Roy/'

    # fig, ax = plt.subplots()

    # plt.hist(data1, bins=bins, histtype='step', normed=True, color='b', label='Region 1', weights=weight1)
 
    # plt.hist(data2, bins=bins, histtype='step', normed=True, color='r', label='Region 3', weights=weight2)

    # plt.legend(fontsize=11)

    # # plt.axis([-1.1, 1.1, 0, 1.7])

    # plt.xlabel('Coseno director respecto al eje Y', fontsize=14)

    # plt.grid(True)

    # # fig.savefig(path + 'Coseno_dir_reg3_y.png', dpi = 200, bbox_inches='tight')

    # plt.show()

    histo1, bin_edges = np.histogram(data1, bins=bins, weights=weight1)
    histo2, bin_edges = np.histogram(data2, bins=bins, weights=weight2)

    bin_centres = (bin_edges[:-1] + bin_edges[1:])/2
    p0 = [1., 0., 1.]

    coeff1, var_matrix1 = curve_fit(gauss, bin_centres, histo1, p0=p0)
    coeff2, var_matrix2 = curve_fit(gauss, bin_centres, histo2, p0=p0)

    print coeff1,coeff2


##################################################################################################################################################################

def DevuelveCosenosOAngulos_PorEnergia(data, type, maxEnergy, bins, variable, axis, energy_ref, porcentaje_energy_ref): 
    energy      = np.abs(data[:, 0])
    PosX        = data[:, 1] 
    PosY        = data[:, 2]
    weight      = data[:, 5]

    delta_E = porcentaje_energy_ref * energy_ref

    if (variable == 'cos' and axis == 'x'):
        xCos = data[:,3]
        isElectron    = 1*(type == 2) * (energy < maxEnergy) * ( abs(energy - energy_ref) < delta_E) * (abs(xCos) < 1.0)                            
        weight = weight*isElectron
        xCos = xCos[weight>0]
        weight = weight[weight>0]
        return xCos, weight

    if (variable == 'cos' and axis == 'y'):
        yCos = data[:,4]
        isElectron    = 1*(type == 2) * (energy < maxEnergy) * ( abs(energy - energy_ref) < delta_E) * (abs(yCos) < 1.0)                           
        weight = weight*isElectron
        yCos = yCos[weight>0]
        weight = weight[weight>0]
        return yCos, weight

def DevuelveCosenosOAngulosDeArchivo_PorEnergia(data_path, maxEnergy, bins, variable, axis, energy_ref, porcentaje_energy_ref):   
    phsp_files   = findFiles(data_path)                                
    for phsp in phsp_files:
        particles, packFormat = parseFile(phsp)                                                   
        
        with open(phsp + '.IAEAphsp', 'rb') as IAEAphsp:                                        
            pack       = np.fromfile(IAEAphsp, dtype=packFormat, count=-1)             
            data, weight   = DevuelveCosenosOAngulos_PorEnergia(pack['Data'], pack['Type'], maxEnergy, bins, variable, axis, energy_ref, porcentaje_energy_ref)
    return data, weight     
    
def GraficaHistogramaCosenosAngulos_PorEnergia(data1, data2, data3, weight1, weight2, weight3, bins):                                                              # CAMBIAR ACA LA LEYENDA SEGUN X E Y

    path = 'C:\Users\Install\Desktop\Tesis Roy/'

    fig, ax = plt.subplots()

    plt.hist(data1, bins=bins, histtype='step', normed=True, color='b', label='E = 0.46 MeV', weights=weight1)
 
    plt.hist(data2, bins=bins, histtype='step', normed=True, color='g', label='E = 2.50 MeV', weights=weight2)

    plt.hist(data3, bins=bins, histtype='step', normed=True, color='r', label='E = 5.00 MeV', weights=weight3)

    plt.legend(fontsize=11)

    # plt.axis([-1.1, 1.1, 0, 1.7])

    plt.xlabel('Coseno director respecto al eje Y', fontsize=14)

    plt.grid(True)

    fig.savefig(path + 'Coseno_dir_energias_y.png', dpi = 200, bbox_inches='tight')

    plt.show()

##################################################################################################################################################################

def Genera_cosenos_txt_horizontal(data_path, folder_path, maxEnergy,bins):
    for i in range(8):
        x_inf = -30.0 + 4.0*i 
        x_sup = -30.0 + 4.0 + 4.0*i
        y_inf = -4.0
        y_sup = 4.0

        variable = 'cos'

        axis = 'x'

        data_x, weight_x = DevuelveCosenosOAngulosDeArchivo_PorRegion(data_path, maxEnergy, bins, variable, axis, x_inf, x_sup, y_inf, y_sup)
        GuardaDataCosenos(folder_path + 'Histograma_horiz_' + axis + '_' + str(i) + '.txt', np.c_[data_x,weight_x])

        axis = 'y'

        data_y, weight_y = DevuelveCosenosOAngulosDeArchivo_PorRegion(phsp_path, maxEnergy, bins, variable, axis, x_inf, x_sup, y_inf, y_sup)
        GuardaDataCosenos(folder_path + 'Histograma_horiz_' + axis + '_' + str(i) + '.txt', np.c_[data_y,weight_y])

def Genera_cosenos_txt_diagonal(data_path, folder_path, maxEnergy,bins):
    for i in range(8):
        x_inf = -30.0 + 4.0*i 
        x_sup = -30.0 + 4.0 + 4.0*i
        y_inf = 30.0 - 4.0 - 4.0*i
        y_sup = 30.0 - 4.0*i

        variable = 'cos'

        axis = 'x'

        data_x, weight_x = DevuelveCosenosOAngulosDeArchivo_PorRegion(data_path, maxEnergy, bins, variable, axis, x_inf, x_sup, y_inf, y_sup)
        GuardaDataCosenos(folder_path + 'Histograma_diagon_' + axis + '_' + str(i) + '.txt', np.c_[data_x,weight_x])

        axis = 'y'

        data_y, weight_y = DevuelveCosenosOAngulosDeArchivo_PorRegion(phsp_path, maxEnergy, bins, variable, axis, x_inf, x_sup, y_inf, y_sup)
        GuardaDataCosenos(folder_path + 'Histograma_diagon_' + axis + '_' + str(i) + '.txt', np.c_[data_y,weight_y])

def GraficaHistogramasSegunDistanciaIsocentro(data1, data2, data3, data4, data5, data6, data7, data8, 
                                                weight1, weight2, weight3, weight4, weight5, weight6, weight7, weight8, 
                                                bins):                                     

    path = 'C:\Users\Install\Desktop\Tesis Roy/'

    fig, ax = plt.subplots()

    plt.hist(data1, bins=bins, histtype='step', normed=True, label='28.07 cm', weights=weight1)
 
    plt.hist(data2+0.5, bins=bins, histtype='step', normed=True, label='24.08 cm', weights=weight2)

    plt.hist(data3+1, bins=bins, histtype='step', normed=True, label='20.10 cm', weights=weight3)

    plt.hist(data4+1.5, bins=bins, histtype='step', normed=True, label='16.12 cm', weights=weight4)
 
    plt.hist(data5+2, bins=bins, histtype='step', normed=True, label='12.17 cm', weights=weight5)

    plt.hist(data6+2.5, bins=bins, histtype='step', normed=True, label='8.25 cm', weights=weight6)

    plt.hist(data7+3, bins=bins, histtype='step', normed=True, label='4.47 cm', weights=weight7)
 
    plt.hist(data8+3.5, bins=bins, histtype='step', normed=True, label='0 cm', weights=weight8)

    plt.legend(fontsize=11)

    # plt.axis([-1.1, 9, 0, 3])

    plt.xlabel('Coseno director respecto al eje X', fontsize=14)

    plt.grid(True)

    fig.savefig(path + 'Coseno_dir_dist_isocentro_diagon_x.png', dpi = 200, bbox_inches='tight')

    plt.show()

##################################################################################################################################################################

def AbreDataTXT(file):
    with open(file) as f:
        data = np.loadtxt(f, dtype='double', comments="#", skiprows=0)
    return data

def GuardaDataCosenos(file, data): np.savetxt(file, data)

# -------------------------------------------------------------------------- MAIN ------------------------------------------------------------------

#-------------------------------------------------------------------------- 25/11 ------------------------------------------------------------------

# path = 'C:\Users\Install\Documents\EGSnrc\Phase Space/'

# filename_homog_cos_x = 'pruebacosdir_homogeneo_x.txt'

# filename_homog_cos_y = 'pruebacosdir_homogeneo_y.txt'

# filename_homog_ang_x = 'pruebaangdir_homogeneo_x.txt'

# filename_homog_ang_y = 'pruebaangdir_homogeneo_y.txt'


# filename_no_homog_ang_x = 'pruebaangdir_no_homogeneo_x.txt'

# filename_no_homog_cos_x = 'pruebacosdir_no_homogeneo_x.txt'

# findFiles(path)

# bins = 100

# axis = 'x'

# variable = 'cos'

# data1, weight1 = DevuelveCosenosOAngulosDeArchivo(path, 6.0, bins, variable, axis)

# variable = 'ang'

# data2, weight2 = DevuelveCosenosOAngulosDeArchivo(path, 6.0, bins, variable, axis)

# GraficaHistogramaCosenosAngulos(data1, data2, weight1, weight2, bins)

#-------------------------------------------------------------------------- 28/11 ------------------------------------------------------------------

#cosenos directores globales

# path = 'C:\Users\Install\Documents\EGSnrc\Phase Space/'

# bins = 100

# variable = 'cos'

# axis = 'x'

# data1, weight1 = DevuelveCosenosOAngulosDeArchivo(path, 6.0, bins, variable, axis)

# axis = 'y'

# data2, weight2 = DevuelveCosenosOAngulosDeArchivo(path, 6.0, bins, variable, axis)

# GraficaHistogramaCosenosAngulos(data1, data2, weight1, weight2, bins)

#--------------------------------------------

#cosenos directores por region

path = 'C:\Users\Install\Documents\EGSnrc\Phase Space/'

bins = 100

variable = 'cos'

axis = 'y'

#Region1

x_inf = -5.0
x_sup = 5.0
y_inf = -5.0
y_sup = 5.0

data1, weight1 = DevuelveCosenosOAngulosDeArchivo_PorRegion(path, 6.0, bins, variable, axis, x_inf, x_sup, y_inf, y_sup)

#Region2,3

x_inf = -25.0
x_sup = -15.0
y_inf = 15.0
y_sup = 25.0

data2, weight2 = DevuelveCosenosOAngulosDeArchivo_PorRegion(path, 6.0, bins, variable, axis, x_inf, x_sup, y_inf, y_sup)

GraficaHistogramaCosenosAngulos_PorRegion(data1, data2, weight1, weight2, bins)


#--------------------------------------------

#cosenos directores por energia

path = 'C:\Users\Install\Documents\EGSnrc\Phase Space/'

bins = 100

variable = 'cos'

axis = 'y'

#0.46MeV

energy_ref = 0.46

porcentaje_energy_ref = 0.6

data1, weight1 = DevuelveCosenosOAngulosDeArchivo_PorEnergia(path, 6.0, bins, variable, axis, energy_ref, porcentaje_energy_ref)

#2.5 MeV

energy_ref = 2.5

porcentaje_energy_ref = 0.2

data2, weight2 = DevuelveCosenosOAngulosDeArchivo_PorEnergia(path, 6.0, bins, variable, axis, energy_ref, porcentaje_energy_ref)

#5 MeV

energy_ref = 5.0

porcentaje_energy_ref = 0.1

data3, weight3 = DevuelveCosenosOAngulosDeArchivo_PorEnergia(path, 6.0, bins, variable, axis, energy_ref, porcentaje_energy_ref)

GraficaHistogramaCosenosAngulos_PorEnergia(data1, data2, data3, weight1, weight2, weight3, bins)

#--------------------------------------------

#cosenos directores en funcion de la distancia al isocentro

# phsp_path = 'C:\Users\Install\Documents\EGSnrc\Phase Space/'

# folder_path = 'C:\Users\Install\Documents\EGSnrc\Phase Space\dist_isocentro_nuevo/'

# bins = 100

# # Genera_cosenos_txt_horizontal(phsp_path, folder_path, 6.0,bins)

# # Genera_cosenos_txt_diagonal(phsp_path, folder_path, 6.0,bins)

# # datos1 = AbreDataTXT(folder_path + 'Histograma_horiz_y_0.txt')
# # datos2 = AbreDataTXT(folder_path + 'Histograma_horiz_y_1.txt')
# # datos3 = AbreDataTXT(folder_path + 'Histograma_horiz_y_2.txt')
# # datos4 = AbreDataTXT(folder_path + 'Histograma_horiz_y_3.txt')
# # datos5 = AbreDataTXT(folder_path + 'Histograma_horiz_y_4.txt')
# # datos6 = AbreDataTXT(folder_path + 'Histograma_horiz_y_5.txt')
# # datos7 = AbreDataTXT(folder_path + 'Histograma_horiz_y_6.txt')
# # datos8 = AbreDataTXT(folder_path + 'Histograma_horiz_y_7.txt')

# datos1 = AbreDataTXT(folder_path + 'Histograma_diagon_x_0.txt')
# datos2 = AbreDataTXT(folder_path + 'Histograma_diagon_x_1.txt')
# datos3 = AbreDataTXT(folder_path + 'Histograma_diagon_x_2.txt')
# datos4 = AbreDataTXT(folder_path + 'Histograma_diagon_x_3.txt')
# datos5 = AbreDataTXT(folder_path + 'Histograma_diagon_x_4.txt')
# datos6 = AbreDataTXT(folder_path + 'Histograma_diagon_x_5.txt')
# datos7 = AbreDataTXT(folder_path + 'Histograma_diagon_x_6.txt')
# datos8 = AbreDataTXT(folder_path + 'Histograma_diagon_x_7.txt')

# data1 = datos1[:,0].astype(np.float64)
# weight1 = datos1[:,1].astype(np.float64)
# data2 = datos2[:,0].astype(np.float64)
# weight2 = datos2[:,1].astype(np.float64)
# data3 = datos3[:,0].astype(np.float64)
# weight3 = datos3[:,1].astype(np.float64)
# data4 = datos4[:,0].astype(np.float64)
# weight4 = datos4[:,1].astype(np.float64)
# data5 = datos5[:,0].astype(np.float64)
# weight5 = datos5[:,1].astype(np.float64)
# data6 = datos6[:,0].astype(np.float64)
# weight6 = datos6[:,1].astype(np.float64)
# data7 = datos7[:,0].astype(np.float64)
# weight7 = datos7[:,1].astype(np.float64)
# data8 = datos8[:,0].astype(np.float64)
# weight8 = datos8[:,1].astype(np.float64)

# GraficaHistogramasSegunDistanciaIsocentro(data1, data2, data3, data4, data5, data6, data7, data8, weight1, weight2, weight3, weight4, weight5, weight6, weight7, weight8, bins)

