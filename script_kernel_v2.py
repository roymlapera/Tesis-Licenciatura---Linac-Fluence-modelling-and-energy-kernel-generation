import numpy as np


def genera_regiones_mitad_agua_mitad_aire(nc_water,nr):
    
    nc_total        = nc_water + 2
    
    start_region    = np.zeros(nr*2)
    stop_region     = np.zeros(nr*2)
    for i in range(nr):
        start_region[2*i  ] = nc_total * i + 2
        start_region[2*i+1] = nc_total * i + 2 + nc_water
        stop_region[2*i   ] = nc_total * i + 2 + nc_water - 1
        stop_region[2*i+1]  = nc_total * i + 2 + nc_water + 1

    start_region = [ str(int(x)) for x in start_region ]
    stop_region = [ str(int(x)) for x in stop_region ]
    
    # print(start_region)
    # print(stop_region)

    # print( ',\n'.join([', '.join(i) for i in zip(*[iter(start_region)]*10)]) )
    print( ',\n'.join([', '.join(i) for i in zip(*[iter(stop_region)]*10)]) )

def genera_regiones_mitad_agua_mitad_aire_2(nc_water,nr):                           #Devuelve regiones de agua en una configuracion de regiones con mitad agua mitad aire
    
    nc_total        = nc_water + 2
    
    start_region    = np.zeros(nr)
    stop_region     = np.zeros(nr)

    for i in range(nr):
        start_region[i] = nc_total * i + 2
        stop_region[ i] = nc_total * i + 2 + nc_water - 1

    start_region = [ str(int(x)) for x in start_region ]
    stop_region = [ str(int(x)) for x in stop_region ]

    # print( ',\n'.join([', '.join(i) for i in zip(*[iter(start_region)]*10)]) )
    print( ',\n'.join([', '.join(i) for i in zip(*[iter(stop_region)]*10)]) )

# genera_regiones_mitad_agua_mitad_aire(90, 210)
# genera_regiones_mitad_agua_mitad_aire_2(90, 210)

a = np.tile([2,1], 210)
print len(a)
a = [ str(int(x)) for x in a ]
print( ',\n'.join([', '.join(i) for i in zip(*[iter(a)]*10)]) )
