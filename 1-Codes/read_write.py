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
# import matplotlib.pyplot as plt
import numpy as np
from pyevtk.hl import gridToVTK

#read the files
A= [float(i) for i in range(61)]
for i in A:
    with open('C:/Users/masoudi/Documents/GitHub/LBM-PF/1-Codes/RESULTS/P2/results_%s.pkl' %(i), 'rb') as f:
    # with open('Z:/Multiphase LBM/Results/Benchmarks/2-RT instability/LDR/P2/results_%s.pkl' %(i), 'rb') as f:
        results = pickle.load(f)
    
    NY = len(results['phase_field'])                # domain size along y axis
    NX = len(results['phase_field'][0])             # domain size along x axis
    NZ=1
    nx, ny,nz = NX,NY,NZ
    phi = results['phase_field']
    phi = np.reshape(phi,[NY,NX,NZ])

    x = np.linspace(0,NX, (nx+1)) 
    y = np.linspace(0,NY, (ny+1))
    z = np.linspace(0,NZ, (nz+1))

    X,Y ,Z = np.meshgrid(x,y,z)
    
    
    # gridToVTK("Z:/Multiphase LBM/Results/Benchmarks/2-RT instability/LDR/P2/LDR_P2_%s" %(int(i)), X, Y, Z, cellData={"Phi": phi})
    gridToVTK("C:/Users/masoudi/Documents/GitHub/LBM-PF/1-Codes/RESULTS/P2/LDR_P2_%s" %(int(i)), X, Y, Z, cellData={"Phi": phi})
