import numpy 		as np
from glob           import glob 
from os.path        import splitext, split, dirname
import sys

directorio_scripts = 'C:/Users/OWner/Documents/Facu/Bariloche 2/Viejos Scripts/'

sys.path.insert(0, directorio_scripts)

import funciones_con_histogramas

# ------------------------------------------------------------------------------------------------------------------------------------------------------------------------

def findFiles(phsp_path,):                                                                           # FINDS FILES INSIDE FOLDER WITH IAEAphsp EXTENSION
    find_names  = phsp_path + '/*.IAEAphsp'                                                       # Names are assumed to be .......... .IAEAheader
    all_names   = glob(find_names)                                                                  # Found files with the searched names
    
    phsp_files  = []                                                                                # List to put correct names
    for name in all_names: phsp_files.append(splitext(name)[0])                                     # Appends name without extension to files
    
    print ('Espacios de fases encontrados:')

    print ('\n'.join(phsp_files))

    return phsp_files

def parseFile(phsp_file):                                                                           # PARSES IAEA FILE FOR PARTICLES AND PACK FORMAT
    phsp_header = phsp_file + '.IAEAheader'                                                         # Header file with full information

    with open(phsp_header, 'r') as fop: lines = fop.readlines()                                     # Opens file to operate. IT will closed automatically
    for index, line in enumerate(lines):                                                            # Search in every line for keywords
        if ('$BYTE_ORDER' in line) and ('1234' in lines[index + 1]):                                # Determins the type of unpack format
            packFormat  = np.dtype([('Type', '<i1'), ('Data', '<f4', (7)), ('Extra', '<i4', (2))])   # Assuming IAEA format with ZLAST
        if ('$BYTE_ORDER' in line) and ('4321' in lines[index + 1]):
            packFormat  = np.dtype([('Type', '>i1'), ('Data', '>f4', (7)), ('Extra', '>i4', (2))])   # Assuming IAEA format with ZLAST
        if '$PARTICLES' in line:  particles = int(lines[index + 1])                                 # Total number of particles in file

    # different_EF_w1 = 'C:/Users/FisMedCom 01/Documents/Tesis Maestria Roy/Tesis primer cuatri/Espacios de fases\VarianClinaciX_6MV_20x20_w1'

    # if ( (phsp_file == different_EF_w1) ):
    #     with open(phsp_header, 'r') as fop: lines = fop.readlines()                                     # Opens file to operate. IT will closed automatically
    #     for index, line in enumerate(lines):                                                            # Search in every line for keywords
    #         if ('$BYTE_ORDER' in line) and ('1234' in lines[index + 1]):                                # Determins the type of unpack format
    #             packFormat  = np.dtype([('Type', '<i1'), ('Data', '<f4', (7)), ('Extra', '<i4', (2))])   # Assuming IAEA format with ZLAST
    #         if ('$BYTE_ORDER' in line) and ('4321' in lines[index + 1]):
    #             packFormat  = np.dtype([('Type', '>i1'), ('Data', '>f4', (7)), ('Extra', '>i4', (2))])   # Assuming IAEA format with ZLAST
    #         if '$PARTICLES' in line:  particles = int(lines[index + 1])                                 # Total number of particles in file
    # else:
    #     with open(phsp_header, 'r') as fop: lines = fop.readlines()                                     # Opens file to operate. IT will closed automatically
    #     for index, line in enumerate(lines):                                                            # Search in every line for keywords

    #         print ('aaaaaaaa')

    #         if ('$BYTE_ORDER' in line) and ('1234' in lines[index + 1]):                                # Determins the type of unpack format
    #             packFormat  = np.dtype([('Type', '<i1'), ('Data', '<f4', (7)), ('Extra', '<i4', (1))])   # Assuming IAEA format with ZLAST
    #         if ('$BYTE_ORDER' in line) and ('4321' in lines[index + 1]):
    #             packFormat  = np.dtype([('Type', '>i1'), ('Data', '>f4', (7)), ('Extra', '>i4', (1))])   # Assuming IAEA format with ZLAST
    #         if '$PARTICLES' in line:  particles = int(lines[index + 1])                                 # Total number of particles in file
    
    return particles, packFormat

# ------------------------------------------------------------------------------------------------------------------------------------------------------------------------

def DevuelveAngulos(data, type_part, maxEnergy, xinf, xsup, yinf, ysup):              # Devuelve las energias de los electrones dentro de [xinf,xsup]x[yinf,ysup] y su peso correspondiente
    energy      = np.abs(data[:, 0])
    PosX        = data[:, 1] 
    PosY        = data[:, 2]
    CosDirX     = data[:, 3]
    CosDirY     = data[:, 4]
    weight      = data[:, 5] 

    isElectron    = 1* (type_part == 2) * (energy < maxEnergy) * (xinf<PosX) * (PosX<xsup) * (yinf<PosY) * (PosY<ysup)

    CosDirZ     = np.sqrt( 1 - CosDirX**2 - CosDirY**2 )

    weight = weight*isElectron*energy/CosDirZ          # para pasar de fluencia a fluencia energ, hay que multiplicar por energy

    angulos = np.arccos(CosDirZ)

    angulos = angulos[weight>0]
    weight = weight[weight>0]

    return angulos, weight

def DevuelveAngulosDeArchivo(data_path, maxEnergy, tamano_pack_particulas, xinf, xsup, yinf, ysup,nro_EF):      # HAce lo mismo que DevuelveEnergias pero para cada pack y junta todos los datos
    phsp_files   = findFiles(data_path) 

    phsp_files = [phsp_files[nro_EF]]                                                       #ACA ESTOY ELIGIENDO QUE EF USAR

    print ('Abriendo:')
    print ('\n'.join(phsp_files))

    for phsp in phsp_files:
        particles, packFormat = parseFile(phsp)  

        if(tamano_pack_particulas==-1):
            with open(phsp + '.IAEAphsp', 'rb') as IAEAphsp:                                        
                pack       = np.fromfile(IAEAphsp, dtype=packFormat, count=-1)    #count me dice la cantidad de items que leo                
                angulos, weight   = DevuelveAngulos(pack['Data'], pack['Type'], maxEnergy, xinf, xsup, yinf, ysup)              
        else:
            angulos = []
            weight = []  
            for x in range(particles//tamano_pack_particulas):
                with open(phsp + '.IAEAphsp', 'rb') as IAEAphsp:            
                    pack       = np.fromfile(IAEAphsp, dtype=packFormat, count=tamano_pack_particulas)    #count me dice la cantidad de items que leo                
                    a =  DevuelveAngulos(pack['Data'], pack['Type'], maxEnergy, xinf, xsup, yinf, ysup) 
                    angulos += list(a[0])
                    weight += list(a[1])
                    sys.stdout.write("\r{.2f}".format((float(x)/(particles//tamano_pack_particulas))*100) + "%")
                    sys.stdout.flush()  
                
    return angulos, weight 

def Genera_PesoAngular(path, filename, bins, maxEnergy):                       # Toma el path y el nombre de archivo, lo abre, grafica y devuelve los pesos estadisticos ya normalizados
    peso_angular = AbreHistograma(path + filename)

    fig = plt.figure() 

    plt.plot(peso_angular, "r.")

    plt.show()

    print(peso_angular.sum() * np.pi/180.0)
    return peso_angular

# ------------------------------------------------------------------------------------------------------------------------------------------------------------------------

def DevuelveFluenciaEnergeticaYPosXY(data, type_part, maxEnergy):                                                     #  POS X Y DE LAS PARTICULAS DE PACK DATA
    PosXparticula    = data[:, 1]                                                                  
    PosYparticula    = data[:, 2]   
    energy           = np.abs(data[:, 0])
    CosDirX     = data[:, 3]
    CosDirY     = data[:, 4]
    weight           = data[:, 5] 

    isElectron  = 1*(type_part == 2) * (energy < maxEnergy)

    CosDirZ     = np.sqrt( 1 - CosDirX**2 - CosDirY**2 )
    
    PosX = PosXparticula[np.nonzero(isElectron)]                                            # El array PosX contiene la posicion X de las particulas que cumplen isElectron == 1
    PosY = PosYparticula[np.nonzero(isElectron)]

    Fenergetica = weight[np.nonzero(isElectron)] * energy[np.nonzero(isElectron)] / CosDirZ[np.nonzero(isElectron)]           #Fluencia energetica

    return Fenergetica, PosX, PosY

def PuntosXY_Fluenciaenergetica_particulas(data_path, maxEnergy,nro_EF):                                    # ABRE EL ARCHIVO IAEA Y DEVUELVE DOS VECTORES CON POS X E Y DE LAS PARTICULAS  IN small DISC
    VectorX        = [0]
    VectorY        = [0]      
    FE             = [0]                                                                 
    phsp_files   = findFiles(data_path)                                                  # Finds IAEA files to work with

    phsp_files = [phsp_files[nro_EF]]                                                       #ACA ESTOY ELIGIENDO QUE EF USAR

    print ('Abriendo:')
    print ('\n'.join(phsp_files))

    for phsp in phsp_files:
        particles, packFormat = parseFile(phsp)                                 
        
        with open(phsp + '.IAEAphsp', 'rb') as IAEAphsp:                                        
                pack        = np.fromfile(IAEAphsp, dtype=packFormat, count=-1)                # Gets full pack data from one file Cuanto mas chico el count menos electrones considera
                Fenergetica, PosX,PosY   = DevuelveFluenciaEnergeticaYPosXY(pack['Data'], pack['Type'], maxEnergy) 
                VectorX = np.concatenate([VectorX,PosX])
                VectorY = np.concatenate([VectorY,PosY])
                FE = np.concatenate([FE, Fenergetica]) 
    return FE, VectorX, VectorY

# ------------------------------------------------------------------------------------------------------------------------------------------------------------------------

def DevuelveEnergias(data, type_part, maxEnergy, xinf, xsup, yinf, ysup):              # Devuelve las energias de los electrones dentro de [xinf,xsup]x[yinf,ysup] y su peso correspondiente
    energy      = np.abs(data[:, 0])
    PosX        = data[:, 1] 
    PosY        = data[:, 2]
    CosDirX     = data[:, 3]
    CosDirY     = data[:, 4]
    weight      = data[:, 5] 

    isElectron    = 1* (type_part == 2) * (energy < maxEnergy) * (xinf<PosX) * (PosX<xsup) * (yinf<PosY) * (PosY<ysup)

    CosDirZ     = np.sqrt( 1 - CosDirX**2 - CosDirY**2 )

    weight = weight*isElectron#/CosDirZ#*energy                                      # para pasar de fluencia a fluencia energ, hay que multiplicar por energy

    energy = energy[weight>0]
    weight = weight[weight>0]

    return energy, weight

def DevuelveEnergiasDeArchivo(data_path, maxEnergy, tamano_pack_particulas, xinf, xsup, yinf, ysup,nro_EF):      # HAce lo mismo que DevuelveEnergias pero para cada pack y junta todos los datos
    phsp_files   = findFiles(data_path) 

    phsp_files = [phsp_files[nro_EF]]                                                       #ACA ESTOY ELIGIENDO QUE EF USAR

    print ('Abriendo:')
    print ('\n'.join(phsp_files))

    for phsp in phsp_files:
        particles, packFormat = parseFile(phsp)  

        if(tamano_pack_particulas==-1):
            with open(phsp + '.IAEAphsp', 'rb') as IAEAphsp:                                        
                pack       = np.fromfile(IAEAphsp, dtype=packFormat, count=-1)    #count me dice la cantidad de items que leo                
                energy, weight   = DevuelveEnergias(pack['Data'], pack['Type'], maxEnergy, xinf, xsup, yinf, ysup)              
        else:
            energy = []
            weight = []  
            for x in range(particles//tamano_pack_particulas):
                with open(phsp + '.IAEAphsp', 'rb') as IAEAphsp:            
                    pack       = np.fromfile(IAEAphsp, dtype=packFormat, count=tamano_pack_particulas)    #count me dice la cantidad de items que leo                
                    a =  DevuelveEnergias(pack['Data'], pack['Type'], maxEnergy, xinf, xsup, yinf, ysup) 
                    energy += list(a[0])
                    weight += list(a[1])
                    sys.stdout.write("\r{.2f}".format((float(x)/(particles//tamano_pack_particulas))*100) + "%")
                    sys.stdout.flush()  
                
    return energy, weight     

def Guarda_data(archivo,data,weights): np.savetxt(archivo, np.c_[data,weights])

def Abre_energias(file):
    with open(file) as f:
        energias, weights = np.loadtxt(f, dtype='double', comments="#", unpack=True)
    return energias, weights

# ----------------------------------------------------------------------------------------------------------------------------------------------------------------------

def Edita_pack(data, type, maxEnergy, xinf, xsup, yinf, ysup):                                                        # Devuelve los angulos directores respecto a z y los pesos estadisticos de los e-  

    electron_en_reg_central = 1*(type == 2) * (xinf<data[:,1]) * (data[:,1]<xsup) * (yinf<data[:,2]) * (data[:,2]<ysup)

    electron_fuera_reg_central = 1*(type == 2) * ((xinf>data[:,1]) + (data[:,1]>xsup) + (yinf>data[:,2]) + (data[:,2]>ysup))

    data[electron_en_reg_central>0,1] = 0
    data[electron_en_reg_central>0,2] = 0

    data[electron_fuera_reg_central>0,0] = 0.511

    return data

def Crea_Nuevo_EF(data_path, maxEnergy, tamano_pack_particulas, xinf, xsup, yinf, ysup):                                                  # HAce lo mismo que DevuelveAngulos pero para cada pack y junta todos los datos
    phsp_files   = findFiles(data_path)
    file = 'C:/Users/FisMedCom 01/Documents/Tesis Maestria Roy/Tesis segundo cuatri/Nuevos Scripts/Creacion del EF para el kernel de comparacion\VarianClinaciX_6MV_20x20_w1'

    if(file in phsp_files):
        phsp_files = [file]
        print ('Abriendo:')
        print ('\n'.join(phsp_files))
    else:
        phsp_files = []
        print('No esta el archivo.')
    
    for phsp in phsp_files:
        particles, packFormat = parseFile(phsp)                                             
        
        with open(phsp + '.IAEAphsp', 'rb') as IAEAphsp:                                        
            pack       = np.fromfile(IAEAphsp, dtype=packFormat, count=-1)    #count me dice la cantidad de items que leo  
            pack['Data']   = Edita_pack(pack['Data'], pack['Type'], maxEnergy, xinf, xsup, yinf, ysup)
        with open(phsp + '_part_colapsadas.IAEAphsp', 'wb') as nuevoIAEAphsp: 
            pack.tofile(nuevoIAEAphsp, sep='')   

# ------------------------------------------------------------------------------------------------------------------------------------------------------------------------

def DevuelveFluenciaEnergeticaYEnergiaYAngulo(data, type_part, maxEnergy):                                                     #  POS X Y DE LAS PARTICULAS DE PACK DATA  
    energy      = np.abs(data[:, 0])
    CosDirX     = data[:, 3]
    CosDirY     = data[:, 4]
    weight      = data[:, 5] 

    isElectron  = 1*(type_part == 2) * (energy < maxEnergy)

    CosDirZ     = np.sqrt( 1 - CosDirX**2 - CosDirY**2 )
    
    Energia = energy[np.nonzero(isElectron)]                                            # El array PosX contiene la posicion X de las particulas que cumplen isElectron == 1
    CosDirZ = CosDirZ[np.nonzero(isElectron)]

    Angulo = np.arccos(CosDirZ)

    Fenergetica = weight[np.nonzero(isElectron)] * energy[np.nonzero(isElectron)]           #Fluencia energetica

    return Fenergetica, Energia, Angulo

def PuntosEnergiaAngulo_Fluenciaenergetica_particulas(data_path, maxEnergy,nro_EF):                                    # ABRE EL ARCHIVO IAEA Y DEVUELVE DOS VECTORES CON POS X E Y DE LAS PARTICULAS  IN small DISC                                                               
    phsp_files   = findFiles(data_path)                                                  # Finds IAEA files to work with

    phsp_files = [phsp_files[nro_EF]]                                                       #ACA ESTOY ELIGIENDO QUE EF USAR

    print ('Abriendo:')
    print ('\n'.join(phsp_files))

    for phsp in phsp_files:
        particles, packFormat = parseFile(phsp)                                 
        
        with open(phsp + '.IAEAphsp', 'rb') as IAEAphsp:                                        
                pack        = np.fromfile(IAEAphsp, dtype=packFormat, count=-1)                # Gets full pack data from one file Cuanto mas chico el count menos electrones considera
                Fenergetica, Energia,Angulo   = DevuelveFluenciaEnergeticaYEnergiaYAngulo(pack['Data'], pack['Type'], maxEnergy)

    print(Energia.min(),Energia.max())
    print(Angulo.min(),Angulo.max())
    print(Fenergetica.min(),Fenergetica.max())

    return Fenergetica, Energia,Angulo








