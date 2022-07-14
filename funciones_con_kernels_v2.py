import numpy as np
import matplotlib.pyplot as plt
import matplotlib.cm as cm
import matplotlib.colors as colors
import scipy.ndimage
import math
import sys
from scipy.interpolate import interp1d
from scipy.interpolate import interp2d

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

	L_min = -7
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
	polarImg = ax.contourf(thetaAxisMesh, rAxisMesh, kernel2draw_ext,L, cmap=cm.jet, norm=colors.LogNorm())

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

def Devuelve_coordenadas_centro_voxeles(rAxis,delta_theta,delta_phi):
    conos_salteados_th = int(delta_theta / (np.pi/180))

    conos_salteados_phi = int(delta_phi / (np.pi/180))

    # print(conos_salteados_th,conos_salteados_phi)

    deltaR = np.zeros(len(rAxis))
    deltaR[0] = rAxis[0]

    k=0
    while k < len(rAxis)-1:
        deltaR[k+1] = rAxis[k+1]-rAxis[k]
        k+=1
    r = [0] + rAxis
    r = r - deltaR/2

    th = np.arange(conos_salteados_th/2,180,conos_salteados_th) * np.pi/180 

    phi = np.arange(conos_salteados_phi/2,360,conos_salteados_phi) * np.pi/180 

    # print(r,th*180/np.pi,phi*180/np.pi)

    return r,th,phi

def Theta_prima(theta,alpha,beta):
	# print(len(theta),len(alpha),len(beta))
	Theta,Alpha,Beta = np.meshgrid(alpha,theta,beta)

	Theta2 = np.zeros(Theta.shape)

	Theta2 = np.arccos( np.cos(Theta)*np.cos(Alpha) - np.cos(Beta)*np.sin(Alpha)*np.sin(Theta) )

	# print(Theta2.min(),Theta2.max())

	return Theta2

def Crea_kernel_continuo(i, theta, kernel):
	return lambda x: interp1d(theta, kernel[i,:], kind='cubic', bounds_error=False, fill_value=(kernel[i,0],kernel[i,-1]))(x)

def CreaKernelPoliangular(rAxis, kernel, p, NC_th, NR, delta_alpha, delta_beta,alpha_max):
	kernel_poliangular = np.zeros((NR, NC_th))                                                  #Crea kernel poliangular

	delta_theta = np.pi/180
	delta_phi = np.pi/180

	rAxis = rAxis[:NR]

	rad, theta, phi = Devuelve_coordenadas_centro_voxeles(rAxis,delta_theta,delta_phi)

	irmax = len(rad)
	ithmax = len(theta)

	alpha = np.arange(0,alpha_max,delta_alpha) + delta_alpha/2
	beta = np.arange(0,2*np.pi,delta_beta) + delta_alpha/2

	# print(alpha*180/np.pi,beta*180/np.pi)

	Theta2 = Theta_prima(theta,alpha,beta)

	kernel_interpolado = np.zeros((Theta2.shape[1],Theta2.shape[2]))

	for k,r in enumerate(rad[:irmax]):                           #r
		for l,th in enumerate(theta[:ithmax]):       #theta
			for m,a in enumerate(alpha):
				for n,b in enumerate(beta):
					kernel_interpolado[m][n] = Crea_kernel_continuo(k, theta, kernel)(Theta2[l,m,n])
			
			kernel_poliangular[k,l] = np.array([ ( kernel_interpolado.sum(axis=1)*np.cos(alpha)+(kernel_interpolado*np.cos(beta)).sum(axis=1)*np.sin(alpha)/np.tan(th) ) * p(a) for a in alpha]).sum()
			
			if(kernel_poliangular[k,l]<=0):
				print('WARNING')

			sys.stdout.write("\r{:.2f}".format( (l+k*len(theta[:ithmax]))*100/(len(rad[:irmax])*len(theta[:ithmax])) ) + "%")
			sys.stdout.flush() 

	kernel_poliangular *= delta_alpha*delta_beta

	return kernel_poliangular

def Theta_prima_func(theta,alpha,beta):
	return np.arccos(np.cos(theta)*np.cos(alpha) + np.sin(alpha)*np.sin(theta)*np.cos(beta))

def RotaKernel(kernel,angulo_solido,alpha,beta,Theta_prima_func,rAxis,delta_theta,delta_phi):
	kernel_rotado = np.zeros((len(rAxis),int(np.pi/delta_theta)))

	rad, theta, phi = Devuelve_coordenadas_centro_voxeles(rAxis,delta_theta,delta_phi)

	for k in range(len(rad)):
		for l,th in enumerate(theta):
			theta2 = Theta_prima_func(th,alpha,beta)
			kernel_rotado[k][l] =  interp1d(theta, kernel[k,:], kind='cubic', bounds_error=False, fill_value='extrapolate')(theta2) #* (2*delta_phi*np.sin(theta2)*np.sin(delta_theta/2) / angulo_solido[l] )
			# sys.stdout.write("\r{:.2f}".format( (l+k*len(theta))*100/(len(rad)*len(theta)) ) + "%")
			# sys.stdout.flush()

	print('1 ')
	return kernel_rotado

def Integral_vol_del_kernel(kernel,rAxis,thetaAxis):
    '''Calcula la integral volumetrica del kernel mediante el calculo de una matriz de volumenes esfericos, la funcion imprime por pantalla el resultado de la integral
    y devuelve la matriz de volumenes con ejes segun rAxis thetaAxis y phiAxis.'''
    deltaR = np.zeros(len(rAxis))
    deltaR[0] = rAxis[0]
    k=0
    while k < (len(rAxis)-1):
        deltaR[k+1] = rAxis[k+1]-rAxis[k]
        k+=1

    RAxis = np.concatenate((np.array([0.0]), rAxis[0:-1]) ,axis=0)
    ThetaAxis = np.concatenate((np.array([0.0]), thetaAxis[0:-1]) ,axis=0) 

    volume = np.zeros(kernel.shape)
    ang_sol = np.zeros(len(thetaAxis))

    for k in range(len(thetaAxis)):
        ang_sol[k] = 2*np.pi/360 * (np.cos(ThetaAxis[k])-np.cos(thetaAxis[k]))
        volume[:,k] = ang_sol[k]*(RAxis*deltaR**2+deltaR*RAxis**2+(deltaR**3)/3)

    Frac_energ_esf = kernel*volume

    print('La integral del kernel en la esfera de radio 10 da: ' + str(360*Frac_energ_esf.sum()))

    return volume,ang_sol

def CreaKernelPoliangular_v2(kernel,Theta_prima_func,rAxis,thetaAxis,delta_theta,delta_phi,delta_alpha,delta_beta,p):
	kernel_poliangular = np.zeros(kernel.shape)

	volume, angulo_sol = Integral_vol_del_kernel(kernel,rAxis,thetaAxis)

	alpha = np.arange(0,np.pi/2,delta_alpha) + delta_alpha/2
	beta = np.arange(0,np.pi,delta_beta) + delta_beta/2

	# kernels_colapsados = []

	# for m,a in enumerate(alpha):
	# 	kernel_colapsado.append(np.array([RotaKernel(kernel[:140,:],angulo_sol[:140,:],a,b,Theta_prima_func,rAxis,delta_theta,delta_phi) for b in beta])[:,:90,:].sum(axis=2)) 

	a = 45 * np.pi/180

	kernel_colapsado = np.zeros(kernel.shape)

	for b in beta:
		aux = RotaKernel(kernel,angulo_sol,a,b,Theta_prima_func,rAxis,delta_theta,delta_phi)

		# plotPolarContour_full(aux[:140,:], rAxis[:140], thetaAxis*180/np.pi, 4)

		kernel_colapsado += aux

	kernel_colapsado *= 2

	volume,ang_sol = Integral_vol_del_kernel(kernel_colapsado,rAxis,thetaAxis)

	print((kernel_colapsado*volume).sum())

	return kernel_colapsado
