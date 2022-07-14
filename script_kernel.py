import numpy 		as np
from glob           import glob 
from os.path        import splitext, split, dirname
import matplotlib.pyplot as plt
import re
import fileinput
import sys
import matplotlib.cm as cm
from matplotlib.pyplot import figure, show, rc

#Funcion que cambia el tipo de particula del archivo egsinp

def Cambia_particula_de_egsinp(egsinp_file, old_particle, new_particle):
	old_string = 'INCIDENT PARTICLE= ' + old_particle + '          # electron,photon,positron'

	new_string = 'INCIDENT PARTICLE= ' + new_particle + '          # electron,photon,positron'

	for line in fileinput.input(egsinp_file, inplace=1):
	        if old_string in line:
	            line = line.replace(old_string,new_string)
	        sys.stdout.write(line)

#Funcion que asigna material aire y agua en incidencia normal




#Funcion que saca los caracteres '[', ']' y '%' del txt. Ojo esto esta hecho solo para archivo egslst la parte del kernel

def Saca_porcentaje_y_corchetes_de_egslst(fname):
	string = open(fname).read()
	new_str = re.sub(r'[^a-zA-Z0-9\n\.\-\+]', ' ', string)
	open(fname, 'w').write(new_str)



def Genera_regiones_mitad_agua_mitad_aire(nc,nr):
	start_region = np.arange(2.0,nc*nr+2,nc/2)
	stop_region = np.zeros(2*nr)

	for i in range(nr*2):
		stop_region[i] = start_region[i] + nc/2 - 1
		start_region[i] = start_region[i]
		stop_region[i] = int ( stop_region[i] )

	start_region = [ int(x) for x in start_region ]
	stop_region = [ int(x) for x in stop_region ]

	return start_region,stop_region

def Cambia_string_de_egsinp(egsinp_file, old_string, new_string):
	for line in fileinput.input(egsinp_file, inplace=1):
	        if old_string in line:
	            line = line.replace(old_string,new_string)
	        sys.stdout.write(line)

#------------------------------------------------------------  MAIN  -----------------------------------------------------------------------------

# egsinp_file = 'C:\EGSnrc\EGS_HOME\edknrc/edknrc_electron_3MeV.egsinp'

# Cambia_particula(egsinp_file,'photon','electron')

# data = np.loadtxt(fname, dtype=np.float32, comments='#', delimiter=None, converters=None, skiprows=5, usecols=None, unpack=False, ndmin=0)



# nc = 48
# nr = 24

# start_region,stop_region =  Genera_regiones_mitad_agua_mitad_aire(nc,nr)

# string_start = ','.join(map(str, start_region))

# string_stop = ','.join(map(str, stop_region))

# print start_region,stop_region


# Cambia_string_de_egsinp(egsinp_file,'START REGION=','START REGION=' + string_start + ' #')

# Cambia_string_de_egsinp(egsinp_file,'STOP REGION=','STOP REGION=' + string_stop + ' #')

txt_file = 'C:\Users\Install\Desktop\edk_nrc_electron_3MeV_50Mhistorias.txt'

# Saca_porcentaje_y_corchetes_de_egslst(txt_file)

data = np.loadtxt(txt_file, dtype=np.float32, comments='#', delimiter=None, converters=None, skiprows=5, usecols=None, unpack=False, ndmin=0)

theta = data[:,0].astype(np.float64)
radii = data[:,1].astype(np.float64)
kernel = data[:,2].astype(np.float64)
theta *= np.pi/180  #de grados a radianes

nc= 48
nr= 24


fig = figure(figsize=(8,8))
ax = fig.add_axes([0.1, 0.1, 0.8, 0.8], polar=True)

#prueba
# theta = np.tile( np.arange(np.pi/nc, np.pi, np.pi/nc) ,nr)
# radii = np.concatenate( ([7.5]*nc,[15.0]*nc,[22.5]*nc,[30.0]*nc))
# kernel = np.arange(0.0,1.0,1.0/(nc*nr))

theta -= np.pi/nc

radii = radii[::-1]
# aux = np.split(radii,nr)
# aux = np.fliplr(aux)
# radii = np.concatenate(aux)

kernel = kernel[::-1]
# aux = np.split(kernel,nr)
# aux = np.fliplr(aux)
# kernel = np.concatenate(aux)

theta_2 = theta + np.pi - np.pi/nc

kernel_2 = kernel
kernel_2 = kernel[::-1]
aux = np.split(kernel_2,nr)
aux = np.fliplr(aux)
kernel_2 = np.concatenate(aux)


width = np.pi/nc
bars = ax.bar(np.concatenate((theta,theta_2)), np.concatenate((radii,radii)), width=width, bottom=0.0)
for k,bar in zip(np.concatenate((kernel,kernel_2)), bars):
    k = k.clip(min=0)
    bar.set_facecolor( cm.hot( 100*k/0.01 ) )
    # print k/0.01
    bar.set_alpha(1)

show()

# fig = figure()

# plt.hist(kernel)

# plt.show()