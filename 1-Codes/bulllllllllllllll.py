# -*- coding: utf-8 -*-
"""
Created on Tue Oct 27 14:06:29 2020

@author: masoudi
"""


#u_in_average = q_in / (ny-2)/dx    # lu_x/lu_t
# implement a parabolic velocity profile at the inlet (like Poiseuille):
#U_max = 3/2*u_in_average
#U_inlet= -U_max /((ny-2)/2)* array([i for i in range(ny-2)]) * (array([i for i in range(1,ny-1)])-(ny-2))


# velocity pre-definition
u_aparent = zeros((ny,nx))
v_aparent = zeros((ny,nx))


# Setting velocity of inlet and outlet boundary condition
#u1[struct==0]=0
#v1[struct==0]=0
#u1[1:ny-1,0] = U_inlet 
#v1[:,0] = 0



# eq.(11) - Mobility = thau_phi*Cs2*dt
mobility = double(0.01) # --> this value is proposed in text Sect. 3.2 and Fakhari et al. (2017a)
thau_phi = mobility/(Cs2*dt)




# force pre-definition
Fb_x = zeros((ny,nx), dtype=double)
Fb_y = zeros((ny,nx), dtype=double)

Fs_x = zeros((ny,nx), dtype=double)
Fs_y = zeros((ny,nx), dtype=double)
force = zeros((ny,nx,n_pop) , dtype=double)
# Surface tension (F_s) calculation --> in text after eq.(7)
Fs_x = chemical_potential * phi_x_1st_derivative 
Fs_y = chemical_potential * phi_y_1st_derivative 
g_b= 1e-6
Fb_x= rho * g_b     

g_pc = zeros((ny,nx,n_pop), dtype=double )
h_pc = zeros((ny,nx,n_pop), dtype=double) # phase field distribuition function



    

g_pc = deepcopy(g_eq)
h_pc = deepcopy(h_eq)

#phi = sum(h, axis=2)

# u calculation --> #eq.(21a)
#g_right = g[:,:,[1,5,8]]
#g_left = g[:,:,[3,6,7]]
#rho_u = sum(g_right, axis=2) - sum(g_left, axis=2)
#u_aparent = rho_u/rho
#u1 = u_aparent/Cs2 + dt*Fs_x/(2*rho)     # eq.(21a)
#    
## v calculation --> #eq.(21a)
#g_up = g[:,:,[2,5,6]]
#g_down = g[:,:,[4,7,8]]
#rho_v = sum(g_up, axis=2) - sum(g_down, axis=2) 
#v_aparent = rho_v/rho  
#v1 = v_aparent/Cs2 + dt*Fs_y/(2*rho)     #eq.(21a)

u1_old = deepcopy(u1)


#writer = ExcelWriter(r"D:/PythonLBM/Free-Energy/force.xlsx" , engine='xlsxwriter')