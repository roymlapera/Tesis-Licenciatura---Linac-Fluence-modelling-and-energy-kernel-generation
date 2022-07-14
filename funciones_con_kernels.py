import numpy as np
import matplotlib.pyplot as plt
import matplotlib.cm as cm
import matplotlib.colors as colors
import scipy.ndimage
import math

def GeneraKernel360(kernel180):                                                                   # A partir de un kernel de 180 grados genera uno de 360 (espejando)
	flippedkernel180 = np.fliplr(kernel180)
	kernel360 = np.concatenate((kernel180, flippedkernel180) ,axis=1)
	return kernel360

def RotaKernel(kernel, nc_rotacion, NR, NC):                                                      # Rota el kernel "nc_rotacion" grados 

	kernel_aux = np.copy(kernel)

	for i in range(nc_rotacion):
		for fila in range(NR):
			kernel_aux[fila,:] = np.concatenate((kernel_aux[fila,-1:],kernel_aux[fila,:-1]))

	return kernel_aux  

def Crea_lista_kernels_rotados(kernel_poli, NR, NC, nc_rotacion):
	kernel_poli = GeneraKernel360(kernel_poli)

	kernels = []

	kernels.append(kernel_poli)

	for i in range(-90//nc_rotacion,90//nc_rotacion-1):					                         # Genero los kernels rotados a partir de rotar el kernel_poli y los meto todos en un vector de kernels llamado "kernels"
		kernel_rotado = RotaKernel(kernel_poli,nc_rotacion*i,NR,NC)
		kernels.append(kernel_rotado)

	return kernels

def plotPolarContour(kernel2draw, rAxis2draw, thetaAxis2draw, Rmax, drawType):

	if drawType == 'full':
		# flippedkernel2draw = np.fliplr(kernel2draw)
		# kernel2draw = np.concatenate((kernel2draw, flippedkernel2draw) ,axis=1)
		thetaAxis2draw = np.concatenate((thetaAxis2draw ,180.0 + thetaAxis2draw),axis=0)

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

	L_min = -12
	L_max = 7
	L = np.logspace(L_min, L_max, 1000)

	l_exp = np.arange(L_min, L_max+1)

	l = 10.0**l_exp

	# print l
	# polarImg = ax.contour(thetaAxisMesh, rAxisMesh, kernel2draw,  cmap=cm.jet, levels = L)
	#plt.title('kernel')

	fig, ax = plt.subplots(subplot_kw=dict(projection='polar'))
	polarImg = ax.contourf(thetaAxisMesh, rAxisMesh, kernel2draw,L, cmap=cm.jet, norm=colors.LogNorm())

	ax.set_theta_zero_location("S")
	#ax.set_theta_direction(-1)
	ax.set_rmax(Rmax)
	ax.grid(True)
	polarImg.set_clim(np.amin(L), np.amax(L))
	cb = plt.colorbar(polarImg) 
	cb.set_ticks(l)
	# path = 'C:\Users\Install\Desktop\Tesis Roy\Imagenes finales/'
	# fig.savefig(path + 'kernel_rotado_45grados.png', dpi = 200, bbox_inches='tight')
	plt.show()

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

	L_min = -7
	L_max = 7
	L = np.logspace(L_min, L_max, 1000)

	l_exp = np.arange(L_min, L_max+1)

	l = 10.0**l_exp

	# print l
	# polarImg = ax.contour(thetaAxisMesh, rAxisMesh, kernel2draw,  cmap=cm.jet, levels = L)
	#plt.title('kernel')

	fig, ax = plt.subplots(subplot_kw=dict(projection='polar'))
	polarImg = ax.contourf(thetaAxisMesh, rAxisMesh, kernel2draw_ext,L, cmap=cm.hot, norm=colors.LogNorm())

	ax.set_theta_zero_location("S")
	#ax.set_theta_direction(-1)
	ax.set_rmax(Rmax)
	ax.grid(True)
	polarImg.set_clim(np.amin(L), np.amax(L))
	cb = plt.colorbar(polarImg) 
	cb.set_ticks(l)
	# path = 'C:\Users\Install\Desktop\Tesis Roy\Imagenes finales/'
	# fig.savefig(path + 'kernel_poliangular.png', dpi = 200, bbox_inches='tight')
	plt.show()

def plotPolarContour_full_360(kernel2draw, rAxis2draw, thetaAxis2draw, Rmax):

	
	# flippedkernel2draw = np.fliplr(kernel2draw)
	# kernel2draw = np.concatenate((kernel2draw, flippedkernel2draw) ,axis=1)
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

	L_min = -7
	L_max = 7
	L = np.logspace(L_min, L_max, 1000)

	l_exp = np.arange(L_min, L_max+1)

	l = 10.0**l_exp

	# print l
	# polarImg = ax.contour(thetaAxisMesh, rAxisMesh, kernel2draw,  cmap=cm.jet, levels = L)
	#plt.title('kernel')

	fig, ax = plt.subplots(subplot_kw=dict(projection='polar'))
	polarImg = ax.contourf(thetaAxisMesh, rAxisMesh, kernel2draw_ext,L, cmap=cm.hot, norm=colors.LogNorm())

	ax.set_theta_zero_location("S")
	#ax.set_theta_direction(-1)
	ax.set_rmax(Rmax)
	ax.grid(True)
	polarImg.set_clim(np.amin(L), np.amax(L))
	cb = plt.colorbar(polarImg) 
	cb.set_ticks(l)
	# path = 'C:\Users\Install\Desktop\Tesis Roy\Imagenes finales/'
	# fig.savefig(path + 'kernel_rotado_90grados.png', dpi = 200, bbox_inches='tight')
	plt.show()

	# print(rAxis2draw.shape)
	# print(kernel2draw.shape)

	# fig = plt.figure()
	# ax = fig.add_subplot(2,1,1)

	# ax.plot(-rAxis2draw,kernel2draw[:,0], marker='.', color='b')

	# ax.set_yscale('log')

	# plt.show()

def CreaKernelPoliangular(kernels, pesos_angulares, NC, NR, nc_rotacion):
	kernel_poliangular = np.zeros((NR, NC*2))                                                     #Crea kernel del doble de conos

	for i in range(-90//nc_rotacion,0):							                                  # hago la suma pesada de kernels polienergeticos en diferentes incidencias
		kernel_poliangular += pesos_angulares[45+i] * kernels[45+i] * nc_rotacion*(np.pi/180)

	for i in range(0,90//nc_rotacion):							                                  # hago la suma pesada de kernels polienergeticos en diferentes incidencias
		kernel_poliangular += pesos_angulares[i] * kernels[i] * nc_rotacion*(np.pi/180)

	kernel_poliangular /= 2 #aca divido el kernel por dos porque los pesos estan normalizados de 0 a 90 grados, no de 0 a 180 grados

	return kernel_poliangular


