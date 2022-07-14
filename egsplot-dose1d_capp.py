#!/usr/bin/python

import sys, re, os.path
import numpy as np
import matplotlib.pyplot as plt


# usage
if (len(sys.argv) < 2):
    print '''\n usage:
    
    %s  x|y|z  a,b  filenames
    
    x|y|z       axis of the dose profile
    a,b         coordinates along the other axes
    filenames   space separated list of files to plot
    
    examples:   egsplot-dose1d.py z 0,0 *.3ddose
                (extract the central z-axis dose profile for all .3ddose files)
    
                egsplot-dose1d.py z 0,0 *.3ddose | xmgrace -settype xydy -pipe
                (pipe the output to xmgrace plotting tool)
    '''%(os.path.basename(sys.argv[0]))
    sys.exit(1)

# simplistic command-line parsing
var   = sys.argv[1]                                     # arg1 is the independent variable: x, y, or z
axis  = ['x', 'y', 'z'].index(var)                      # axis number: x=0, y=1, z=3
where = map(float,sys.argv[2].split(','))               # position of the axis in x,y form
files = sys.argv[3:]                                    # array of all files to process

# print xmgrace header
print '''\
@   title "dose distribution"
@   xaxis label "%s (cm)"
@   yaxis label "dose (Gy)"
@   type xydy''' % var

# loop over all files
count = 0
for file in files:

    # open 3ddose file
    dosefile = open(file, 'r')

    # get voxel counts on first line
    nx, ny, nz = map(int,dosefile.readline().split())    # number of voxels along x, y, z
    Ng = (nx+1) + (ny+1) + (nz+1)                        # total number of voxel grid values (voxels+1) 
    Nd = nx*ny*nz                                        # total number of data points

    # get voxel grid, dose and relative errors
    data  = map(float,dosefile.read().split())           # read the rest of the file
    xgrid = data[:nx+1]                                  # voxel boundaries in x (nx+1 values, 0 to nx)
    ygrid = data[nx+1:nx+1+ny+1]                         # voxel boundaries in y (ny+1 values, nx+1 to nx+1+ny)
    zgrid = data[nx+1+ny+1:Ng]                           # voxel boundaries in z (nz+1 values, rest up to Ng-1)
    dose  = data[Ng:Nd+Ng]                               # flat array of Nd = nx*ny*nz dose values
    # print dose
    errs  = data[Nd+Ng:]                                 # flat array of Nd = nx*ny*nz relative error values
    del data                                             # make sure we don't refer to data by mistake from now on

    # close 3ddose file
    dosefile.close()

    # setup variables for current axis
    grid   = [xgrid, ygrid, zgrid]                       # voxel boundaries in x, y, z
    step   = [1, nx, nx*ny]                              # step between values along x, y, z
    jump   = [nx, nx*ny, nx*ny*nz]                       # distance between start and end voxels
    mygrid = grid[axis]                                  # grid for plot axis
    mystep = step[axis]                                  # step for plot axis
    del grid[axis]                                       # remove plot axis from grid
    del step[axis]                                       # remove plot axis from step

    # get voxel indices for location (along two remaining axes)
    # (say you are plotting along z, then index will hold the indices for the requested x,y positions
    index = []
    for g,w in zip(grid,where):
        if (w<g[0] or w>g[-1]):
            print "ERROR: location", where, "outside of data grid!\n"
            sys.exit(1)
        i = len(g)-1
        while (w < g[i]): i -= 1
        index.append(i)

    # get the actual dose profiles
    start  = index[0]*step[0] + index[1]*step[1]         # starting voxel index
    end    = start + jump[axis]                          # "end" voxel index
    mydose = dose[start:end:mystep]                      # dose slice
    myerrs = errs[start:end:mystep]                      # relative error slice

    # print xmgrace format commands for this set
    print '''\
@   s%d errorbar length 0
@   s%d legend "%s" ''' % (count, count, file)

    ##### just for plotting #####
    doseArray = np.zeros([len(mydose),3])
    
    archivo = open('pdd_phsp_homog.txt','a')
    # print dose profile with absolute errors
    print "#\n# %10s %12s %12s" % (var, "dose", "err")
    for i in range(len(mydose)):
        t = (mygrid[i]+mygrid[i+1])/2.0
        doseArray[i,0] = t
        doseArray[i,1] = mydose[i]
        doseArray[i,2] = myerrs[i]*mydose[i]
        print "%12.6f %12g %12g" % (t, mydose[i], myerrs[i]*mydose[i])
        archivo.write('%12.6f\t%12g\t%12g\n' % (t, mydose[i], myerrs[i]*mydose[i]))
    print
    archivo.close()
    
    count += 1

np.save('PDD_Simulado_15M_a',doseArray)    
plt.plot(doseArray[:,0], doseArray[:,1], "b.")
plt.show()

