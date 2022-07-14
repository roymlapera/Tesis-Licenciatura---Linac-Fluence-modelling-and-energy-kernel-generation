import numpy as np
import matplotlib.pyplot as plt

simulatedTheta = 180 #cantidad de angulos equiespaciados en la generacion MC del kernel
simulatedRadii = 210

header_rows_1 = 216
footer_rows_1 = 8

header_rows_2 = 245
footer_rows_2 = 8

monoEDK1 = np.genfromtxt('electron_fullWater_6MeV_VS_DC.egslst', skip_header=header_rows_1, skip_footer=footer_rows_1, usecols=(0, 1, 2), unpack=True)

monoEDK2 = np.genfromtxt('electron_halfWater_6MeV_VS_DC.egslst', skip_header=header_rows_2, skip_footer=footer_rows_2, usecols=(0, 1, 2), unpack=True)

kernel1 = np.zeros([simulatedRadii, simulatedTheta])

rAxis = monoEDK1[1,0:simulatedRadii]
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
    kernel1[:,k] = monoEDK1[2,start:stop]
    thetaAxis[k] = monoEDK1[0,start]

kernel2 = np.zeros([simulatedRadii, simulatedTheta])

rAxis = monoEDK2[1,0:simulatedRadii]
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
    kernel2[:,k] = monoEDK2[2,start:stop]
    thetaAxis[k] = monoEDK2[0,start]


resta = np.abs(kernel2 - kernel1)

resta = resta[kernel2>0]/kernel2[kernel2>0]

# fig = plt.figure()
fig, ax = plt.subplots()

plt.hist(resta.ravel(), bins = 500)
# rects1=plt.hist(monoEDK1.ravel(), bins = 100)
# rects2=plt.hist(monoEDK2.ravel(), bins = 100)
# print resta.ravel()

# plt.axis([0, 2500, 0, 30])

plt.show()

# print 100.0*len(resta[resta>0.01])/len(resta)

print np.amin(resta)
print np.amax(resta)



