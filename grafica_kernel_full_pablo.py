# -*- coding: utf-8 -*-
#### kernel check


import numpy as np
import matplotlib.pyplot as plt
import matplotlib.cm as cm
import matplotlib.colors as colors
import scipy.ndimage
import random

############################
#### Funcion para graficar
############################
def plotPolarContour_full(kernel2draw, rAxis2draw, thetaAxis2draw, Rmax):

	
	flippedkernel2draw = np.fliplr(kernel2draw)
	kernel2draw = np.concatenate((kernel2draw, flippedkernel2draw) ,axis=1)
	thetaAxis2draw = np.concatenate((thetaAxis2draw ,180.0 + thetaAxis2draw),axis=0)
	# print theta_axis
	thetaAxis2draw_ext = np.append(-thetaAxis2draw[0], thetaAxis2draw)  

	# print theta_axis_ext

	kernel2draw_ext = np.zeros([kernel2draw.shape[0],kernel2draw.shape[1]+1])
	kernel2draw_ext[:,0] = kernel2draw[:,-1]
	kernel2draw_ext[:,1:] = kernel2draw

	deltaR = np.zeros(len(rAxis2draw))
	deltaR[0] = rAxis2draw[0]
	k=0
	while k < len(rAxis2draw)-1:
		deltaR[k+1] = rAxis2draw[k+1]-rAxis2draw[k]
		k+=1

	rAxis2draw=rAxis2draw-deltaR/2
	deltaTheta = np.zeros(len(thetaAxis2draw_ext))
	deltaTheta[0] = thetaAxis2draw_ext[0] 
	k=0

	while k < len(thetaAxis2draw)-1:
		deltaTheta[k+1] = thetaAxis2draw_ext[k+1]-thetaAxis2draw_ext[k]
		k+=1

	thetaAxis2draw_ext=thetaAxis2draw_ext-deltaTheta/2
	thetaAxis2draw_ext = np.radians(thetaAxis2draw_ext)

	thetaAxisMesh, rAxisMesh = np.meshgrid(thetaAxis2draw_ext, rAxis2draw)

	# kernel2draw = scipy.ndimage.zoom(kernel2draw, 3)
	# rAxisMesh = scipy.ndimage.zoom(rAxisMesh, 3)
	# thetaAxisMesh = scipy.ndimage.zoom(thetaAxisMesh, 3)

	# fig = plt.figure()
	# ax = fig.add_subplot(111, polar = True)

	# polarImg = ax.pcolor(thetaAxisMesh, rAxisMesh, kernel2draw, cmap=cm.jet, norm=colors.LogNorm()) #
	# polarImg = ax.pcolormesh(thetaAxisMesh, rAxisMesh, kernel2draw, cmap=cm.jet, norm=colors.LogNorm()) #
	# polarImg = ax.pcolormesh(thetaAxis2draw, rAxis2draw, kernel2draw, cmap=cm.jet, norm=colors.LogNorm())
	# kernel2draw.max()*0.01
	# L = [kernel2draw.max()*0.001, kernel2draw.max()*0.01]

	L_min = -12
	L_max = 4
	L = np.logspace(L_min, L_max, 1000)

	l_exp = np.arange(L_min, L_max+1)

	l = 10.0**l_exp

	# print l
	# polarImg = ax.contour(thetaAxisMesh, rAxisMesh, kernel2draw,  cmap=cm.jet, levels = L)
	#plt.title('kernel')

	fig, ax = plt.subplots(subplot_kw=dict(projection='polar'))
	polarImg = ax.contourf(thetaAxisMesh, rAxisMesh, kernel2draw_ext,L, cmap=cm.jet, norm=colors.LogNorm())

	ax.set_theta_zero_location("S")
	#ax.set_theta_direction(-1)
	ax.set_rmax(Rmax)
	ax.grid(True)
	polarImg.set_clim(np.amin(L), np.amax(L))
	cb = plt.colorbar(polarImg) 
	cb.set_ticks(l)
	# path = 'C:\EGSnrc\EGS_HOME\edknrc/'
	# fig.savefig(path + 'fullWater_50M_0500MeV_prueba.png', dpi = 200, bbox_inches='tight')
	plt.show()


##############################################
#### Get rAxis vector from a edknrc file #####
##############################################
simulatedTheta = 180 #cantidad de angulos equiespaciados en la generacion MC del kernel
simulatedRadii = 210

header_rows = 216
footer_rows = 8

monoEDK = np.genfromtxt('electron_fullWater_2MeV.egslst', skip_header=header_rows, skip_footer=footer_rows, usecols=(0, 1, 2), unpack=True) #, dtype=np.float32

# print monoEDK.shape

kernel = np.zeros([simulatedRadii, simulatedTheta])

rAxis = monoEDK[1,0:simulatedRadii]
thetaAxis = np.zeros(simulatedTheta)

deltaR = np.zeros(len(rAxis))
deltaR[0] = rAxis[0]
k=0
while k < simulatedRadii-1:
    deltaR[k+1] = rAxis[k+1]-rAxis[k]
    k+=1
    
# shift rAxis to the center of the voxel

rAxis_mid = rAxis - deltaR/2

for k in range(simulatedTheta):
#    print k
    start = k*simulatedRadii
    # print start
    stop = (k+1)*simulatedRadii
    # print stop
    kernel[:,k] = monoEDK[2,start:stop]
    thetaAxis[k] = monoEDK[0,start]

    
##############################
#### Check normalization #####
##############################
# debo multiplicar cada elemento por su volumen y sumar
# construction of volume elements in cm**3



RAxis = np.concatenate((np.array([0.0]), rAxis[0:-1]) ,axis=0)
ThetaAxis = np.concatenate((np.array([0.0]), thetaAxis[0:-1]) ,axis=0)

volume = np.zeros(kernel.shape)

for k in range(simulatedTheta):
   volume[:,k] = 2.0*np.pi*( np.cos(np.radians(ThetaAxis[k])) - np.cos(np.radians(thetaAxis[k])))*(RAxis*deltaR**2+deltaR*RAxis**2+(deltaR**3)/3 )
   
kernel = kernel/volume

print kernel.sum()

# kernel = kernel[:,:18900]

# print kernel.shape

# for theta in range(7):
# 	for i in range(210):
# 		if (i>=6):
# 			kernel[i,theta] = kernel[i,random.randint(7,14)] 

# np.savez("fullWater_50M_5625MeV.npz", kernel = kernel, rAxis = rAxis, thetaAxis = thetaAxis)


#################
#### Figures ####
#################   

# print kernel.shape

# print rAxis

# plt.figure()	

# plt.plot(rAxis, kernel[:,10], "g.")
# plt.plot(rAxis, kernel[:,20], "r.")
# plt.plot(rAxis, kernel[:,30], "k.")
# plt.plot(rAxis, kernel[:,40], "y.")
# plt.plot(rAxis, kernel[:,50], "m.")
# plt.plot(rAxis, kernel[:,60], "c.")
# plt.plot(rAxis, kernel[:,70], "b.")
# plt.plot(rAxis, kernel[:,80], "m.")
# plt.plot(rAxis, kernel[:,90], "c.")
# plt.plot(rAxis, kernel[:,100], "c.")        
# plt.plot(rAxis, kernel[:,110], "g.")                                                 
# plt.plot(rAxis, kernel[:,120], "r.")
# plt.plot(rAxis, kernel[:,130], "k.")
# plt.plot(rAxis, kernel[:,140], "y.")
# plt.plot(rAxis, kernel[:,150], "m.")
# plt.plot(rAxis, kernel[:,160], "c.")
# plt.plot(rAxis, kernel[:,170], "b.")

# plt.axis([0, 10, 1e-6, 1e7])

# plt.yscale('log')

# plt.show()


    
plotPolarContour_full(kernel, rAxis, thetaAxis, 10)
# #
# #
# #
# fig = plt.figure()
# ax = fig.add_subplot(111)

# for k in range(simulatedTheta):
#     ax.plot(rAxis_mid, kernel[:,k], '.')

# ax.set_yscale('log')
# plt.show()


# fig = plt.figure()
# ax = fig.add_subplot(111)
# for k in range(simulatedRadii):
#     plt.plot(thetaAxis, kernel[k,:])

# ax.set_yscale('log')
# plt.show()


# # np.savetxt('prueba',kernel)

# #fig = plt.figure()
# #ax = fig.add_subplot(111)
# #ax.plot(rAxis_mid, iDK[:,0], 'r', label='3.75')
# #ax.plot(rAxis_mid, iDK[:,6], 'g--', label='48.75')
# #ax.plot(rAxis_mid, iDK[:,12],'b--', label = ' 93.75' )
# #ax.plot(rAxis_mid, iDK[:,18],'m--', label = ' 138.75' )
# #legend = plt.legend(loc='upper right')
# #plt.xlim([0, 6.0])
# #plt.show()
# #
# #fig = plt.figure()
# #ax = fig.add_subplot(111)
# #ax.plot(rAxis_mid, iC[:,0], 'r', label='3.75')
# #ax.plot(rAxis_mid, iC[:,6], 'g--', label='48.75')
# #ax.plot(rAxis_mid, iC[:,12],'b--', label = ' 93.75' )
# #ax.plot(rAxis_mid, iC[:,18],'m--', label = ' 138.75' )
# #ax.plot(rAxis_mid, niC[:,0], 'r*', label='3.75')
# #ax.plot(rAxis_mid, niC[:,6], 'g*', label='48.75')
# #ax.plot(rAxis_mid, niC[:,12],'b*', label = ' 93.75' )
# #ax.plot(rAxis_mid, niC[:,18],'m*', label = ' 138.75' )
# #legend = plt.legend(loc='upper right')
# #plt.xlim([0, 6.0])
# #plt.show()
# #
# #fig = plt.figure()
# #ax = fig.add_subplot(111)
# #ax.plot(rAxis_mid, iCC[:,0], 'r', label='3.75')
# #ax.plot(rAxis_mid, iCC[:,6], 'g--', label='48.75')
# #ax.plot(rAxis_mid, iCC[:,12],'b--', label = ' 93.75' )
# #ax.plot(rAxis_mid, iCC[:,18],'m--', label = ' 138.75' )
# #ax.plot(rAxis_mid, niCC[:,0], 'r*', label='3.75')
# #ax.plot(rAxis_mid, niCC[:,6], 'g*', label='48.75')
# #ax.plot(rAxis_mid, niCC[:,12],'b*', label = ' 93.75' )
# #ax.plot(rAxis_mid, niCC[:,18],'m*', label = ' 138.75' )
# #legend = plt.legend(loc='upper right')
# #plt.xlim([0, 6.0])
# #plt.show()