# -*- coding: utf-8 -*-
"""
Created on Fri Oct 30 13:44:25 2020

@author: masoudi
"""
import matplotlib.pyplot as plt
import numpy as np
from copy import deepcopy


X = np.arange(0, 52, 1)
Y = np.arange(0, 64, 1)
U=deepcopy(u1)
V=deepcopy(v1)
U, V = np.meshgrid(X, Y)

fig, ax = plt.subplots()
q = ax.quiver(X, Y, U, V)


plt.show()


data_x = [10.0, 8.0, 13.0, 9.0, 11.0, 14.0, 6.0, 4.0, 12.0, 7.0, 5.0]
data_y = [8.04, 6.95, 7.58, 8.81, 8.33, 9.96, 7.24, 4.26, 10.84, 4.82, 5.68]

fig, ax = plt.subplots()

ax.scatter(x=u1[:,5]/np.amax(u1), y=[i for i in range(64)], c="#E69F00")