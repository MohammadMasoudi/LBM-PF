# -*- coding: utf-8 -*-
"""
Created on Wed Oct  7 10:28:08 2020

@author: masoudi
"""

# -*- coding: utf-8 -*-
"""
Created on Mon Jun 29 13:22:13 2020

@author: masoudi
"""
import pickle
import matplotlib.pyplot as plt
import numpy as np
# from pyevtk.hl import gridToVTK

#read the files
A= [float(i) for i in range(22)]
for i in A:
    with open('C:/Users/masoudi/Documents/GitHub/LBM-PF/2-Results/2-Rayleigh-Taylor/LDR/results_%s.pkl' %i, 'rb') as f:
        results = pickle.load(f)
    
    ny = len(results['phase_field'])                # domain size along y axis
    nx = len(results['phase_field'][0])             # domain size along x axis
    Z = results['phase_field']
    x = np.arange(0, nx, 1)
    y = np.arange(0, ny, 1)
    
    fig, ax = plt.subplots()
    im=ax.pcolormesh(x, y, Z)
    fig.colorbar(im)
# plt.show()

#
#with open('C:/Users/masoudi/Desktop/privat/Papers/Nature/Results/results_2/RN.pkl', 'rb') as f:
#    RNkn = pickle.load(f)


#     NY=np.shape(struc)[0]
#     NX=np.shape(struc)[1]
#     NZ=1
#     nx, ny,nz = NX,NY,NZ
    
#     struc = np.reshape(struc,[NX,NY,NZ])
#     struc =struc .astype(int)
#     # Coordinates 
#     x = np.linspace(0,NX, (nx+1)) 
#     y = np.linspace(0,NY, (ny+1))
#     z = np.linspace(0,NZ, (nz+1))
    
#     X,Y ,Z = np.meshgrid(x,y,z)


# #for i in range (1,201):
# #    
# #    with open('C:/Users/masoudi/Desktop/privat/Papers/clogging model/Results_70in70/Results/RESULTS_1/results_1/results_%s.pkl'%(i), 'rb') as f:
# #        dat = pickle.load(f)
# #    
# #    Precipitant_3D = dat[1:NZ-1][:,1:NY-1][:,:,1:NX-1]
# #    final=initialStruc_3D.astype(int)
# #    final [Precipitant_3D]=2

#     gridToVTK("C:/Users/masoudi/Desktop/privat/Papers/Nature/Results/movie/50_NUC_run2_%s"%(i), X, Y, Z, cellData={"PRC": struc})
