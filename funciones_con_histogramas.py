import numpy as np
import matplotlib.pyplot as plt
import csv

def Guarda_datos_angulos(archivo,angulos,weights): np.savetxt(archivo, np.c_[angulos,weights])

def Abre_datos_angulos(file):
    with open(file) as f:
        data = np.loadtxt(f, dtype='double', comments="#")

    angulos = data[:,0]
    weights = data[:,1]
    return angulos, weights

def GeneraHistograma_angulos(data, weight, bins):
	histo, bin_edges = np.histogram(data, bins=bins, range=[0,np.pi/2.0], weights=weight, normed=True)
	return histo, bin_edges

def Guarda_datos_angulos(archivo,angulos,weights): np.savetxt(archivo, np.c_[angulos,weights])

def Guarda_histograma(archivo, archivos_bins, histograma, bin_edges): 
	np.savetxt(archivo, histograma[None, :], delimiter='\t')
	np.savetxt(archivos_bins, bin_edges[None, :], delimiter='\t')

def Abre_archivo_pesos_angulares(file):
    with open(file) as f:
        pesos = np.loadtxt(f, dtype='double', comments="#",delimiter="\t")
    return pesos

# -------------------------------------------------------------------------------------------------------------------------------------

def Grafica_PosXY_fenergetica(x, y, FE, bins):    #ACA poner PosX PosY y FE
    fig = plt.figure()

    # gs = gridspec.GridSpec(2, 2,
    #                    width_ratios=[8,4],
    #                    height_ratios=[4,8])

    # sub1 = fig.add_subplot(gs[0])
    # sub1.set_xticks([])
    # sub1.set_yticks(np.arange(0, max(perfil_x), 0.2))
    # sub1.plot(abcisas, perfil_x, "g.")
    # sub1.tick_params(labelsize=10)
    # sub1.grid()

    # sub2 = fig.add_subplot(gs[3])
    # sub2.set_yticks([])
    # sub2.set_xticks(np.arange(0, max(perfil_y), 0.2))
    # sub2.plot(perfil_y, abcisas, "g.")
    # sub2.tick_params(labelsize=10)
    # sub2.grid()
 
    hist, xedges, yedges = np.histogram2d(x*10.0, y*10.0, bins=[bins,bins], range=[[-300, 300], [-300, 300]], normed=True, weights=FE)

    # np.savetxt('fluencia_energetica_imshow_ejes.txt', np.c_[xedges,yedges])
    # np.savetxt('fluencia_energetica_imshow_histo.txt', hist)

    # binwidth = 2
    # xymax = np.max( [np.max(np.fabs(perfil_x)), np.max(np.fabs(perfil_y))] )
    # lim = ( int(xymax/binwidth) + 1) * binwidth

    
    plt.imshow(hist, extent=[-300,300,-300,300], interpolation='None')
    # sub3.axhline(linewidth=2, color='g')
    # sub3.axvline(linewidth=2, color='g')
    # sub3.axis([-300, 300, -300, 300])
    # sub3.tick_params(labelsize=10)
    # sub3.set_aspect('equal')
    # cm3 = plt.colorbar(im3)

    plt.tight_layout()
    plt.show()
    # fig.savefig(path + 'fluencia energetica.png', dpi = 200, bbox_inches='tight')
    return hist, xedges, yedges

def Guarda_datos_fenergetica(archivo, fenergetica, PosX, PosY): np.savetxt(archivo, np.c_[fenergetica,PosX,PosY])

def Abre_datos_fenergetica(file):
    fenergetica, PosX, PosY = np.loadtxt(file, unpack=True)
    return fenergetica, PosX, PosY










