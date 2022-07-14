import numpy as np
import matplotlib.pyplot as plt
import matplotlib.cm as cm
import matplotlib.colors as colors
import scipy.ndimage

############################
#### Funcion para graficar
############################
def plotPolarContour(kernel2draw, rAxis2draw, thetaAxis2draw, Rmax, drawType):

	# if drawType == 'full':
	# 	flippedkernel2draw = np.fliplr(kernel2draw)
	# 	kernel2draw = np.concatenate((kernel2draw, flippedkernel2draw) ,axis=1)
	# 	thetaAxis2draw = np.concatenate((thetaAxis2draw ,180.0 + thetaAxis2draw),axis=0)

	deltaR = np.zeros(len(rAxis2draw))
	deltaR[0] = rAxis2draw[0]
	k=0
	while k < len(rAxis2draw)-1:
		deltaR[k+1] = rAxis2draw[k+1]-rAxis2draw[k]
		k+=1

	rAxis2draw=rAxis2draw-deltaR/2
	deltaTheta = np.zeros(len(thetaAxis2draw))
	deltaTheta[0] = thetaAxis2draw[0] 
	k=0

	while k < len(thetaAxis2draw)-1:
		deltaTheta[k+1] = thetaAxis2draw[k+1]-thetaAxis2draw[k]
		k+=1

	thetaAxis2draw=thetaAxis2draw-deltaTheta/2
	thetaAxis2draw = np.radians(thetaAxis2draw)

	thetaAxisMesh, rAxisMesh = np.meshgrid(thetaAxis2draw, rAxis2draw)

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

	L = np.logspace(-9, 4, 1000)
	# polarImg = ax.contour(thetaAxisMesh, rAxisMesh, kernel2draw,  cmap=cm.jet, levels = L)
	#plt.title('kernel')

	fig, ax = plt.subplots(subplot_kw=dict(projection='polar'))
	polarImg = ax.contourf(thetaAxisMesh, rAxisMesh, kernel2draw, 1000, cmap=cm.jet)#, norm=colors.LogNorm())

	ax.set_theta_zero_location("S")
	#ax.set_theta_direction(-1)
	ax.set_rmax(Rmax)
	ax.grid(True)
	# polarImg.set_clim(np.amin(L), np.amax(L))
	plt.colorbar(polarImg) 
	path = 'C:\EGSnrc\EGS_HOME\edknrc/'
	fig.savefig(path + 'kernel_error.png', dpi = 200, bbox_inches='tight')
	plt.show()



# kernel_full_water_npz = np.load("electron_fullWater_6MeV_VS_DC.npz")

# kernel_half_water_npz = np.load("electron_halfWater_6MeV_VS_DC.npz")

kernel_full_water_npz = np.load("electron_fullWater_6MeV_VS_DC.npz")

kernel_half_water_npz = np.load("electron_halfWater_6MeV_VS_DC.npz")


kernel_half = kernel_half_water_npz['kernel']
kernel_full = kernel_full_water_npz['kernel']


print kernel_half.shape
print kernel_full.shape

error = np.abs(kernel_half - kernel_full)/kernel_full

rAxis = kernel_full_water_npz['rAxis']
thetaAxis = kernel_full_water_npz['thetaAxis']

print rAxis.shape
print thetaAxis.shape




#################
#### Figures ####
#################   
   
plotPolarContour(error, rAxis, thetaAxis, 10, 'hemi')

# fig, ax = plt.subplots()
# ax.imshow(error)
# plt.show()

print error.max()

print error.shape


