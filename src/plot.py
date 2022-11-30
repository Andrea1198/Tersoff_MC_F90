import matplotlib.pyplot as plt
from numpy import loadtxt
import numpy as np
from mpl_toolkits.mplot3d import Axes3D
dir = './files/output/'
img_dir = './images/'
filenames = [ "gdr_initial.out", "gdr_crystal.out", "gdr_liquid.out", "gdr_glass.out"]
# filenames = ['crystal.csv', 'liquid.csv', 'glass.csv', 'conf_crystal.xyz', 'conf_liquid.xyz', 'conf_glass.xyz']
# show      = [True, True, True, True, True, True]
# save      = [True, True, True, True, True, True]
show = [False] * len(filenames)
save = [False] * len(filenames)
save = [True] * len(filenames)
show = [True] * len(filenames)


for i, file in enumerate(filenames):
    if file[-3:] == "csv":
        r = loadtxt(dir+file, delimiter=',', skiprows=1)
        x = r[:, 0]
        y = r[:, 1]
        z = r[:, 2]
        xo = x[r[:, 3] == 0]
        yo = y[r[:, 3] == 0]
        zo = z[r[:, 3] == 0]
        xSi = x[r[:, 3] == 1]
        ySi = y[r[:, 3] == 1]
        zSi = z[r[:, 3] == 1]
        fig = plt.figure()
        ax = fig.add_subplot(111, projection='3d', proj_type = 'ortho')
        ax.scatter(xo, yo, zo, s=1.5)
        ax.scatter(xSi, ySi, zSi, s=3)
        plt.axis('off')
        if show[i]: plt.show()
        if save[i]: fig.savefig(img_dir+file[:-3]+".pdf")
        
    elif file[-3:] == "xyz":
        list_string = []
        with open(dir+file) as f:
            lines = f.readlines()
            for line in lines:
                ll = line.split(" ")
                if ll[0] == "Si":
                    list_string.append(0)
                elif ll[0] == "O":
                    list_string.append(1)
                
        data = np.genfromtxt(dir+file, skip_header=1)
        x = data[:,1]
        y = data[:,2]
        z = data[:,3]
        typ = np.array(list_string)
        xo = x[typ==0]
        yo = y[typ==0]
        zo = z[typ==0]
        xSi = x[typ==1]
        ySi = y[typ==1]
        zSi = z[typ==1]

        fig = plt.figure()
        ax = fig.add_subplot(111, projection='3d', proj_type='ortho')
        ax.scatter(xo, yo, zo)
        ax.scatter(xSi, ySi, zSi)

        if show[i]: plt.show()
        if save[i]: plt.savefig(img_dir+file[:-3]+".pdf")
    elif file[-3:] == "out":
            
        fil = loadtxt(dir+file, skiprows=1)
        fil = fil[:400,:]
        r = fil[:, 0]
        Si_Si   = fil[:, 1]
        Si_O    = fil[:, 2]
        O_o     = fil[:, 3]
        fig = plt.figure()
        ax = plt.axes()
        ax.scatter(r, Si_Si, s=0.3)
        ax.scatter(r, Si_O, s=0.3)
        ax.scatter(r, O_o, s=0.3)
        if show[i]: plt.show()
        if save[i]: fig.savefig(img_dir+file[:-3]+".pdf")
