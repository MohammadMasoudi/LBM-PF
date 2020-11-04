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




#writer = ExcelWriter(r"D:/PythonLBM/Free-Energy/force.xlsx" , engine='xlsxwriter')
# =============================================================================
#     #Streaming (when there is no streamming from solid boundary nodes)
#     for i in range(nx):
#         for j in range(ny):
#             for k in range(n_pop):
#                 if struct[j,i]:
#                     newx = i + cx[k]
#                     newy = j - cy[k]
#                     if (newx!=-1) and (newx!=nx) and (newy!=-1) and (newy!=ny) :
#                         if struct[newy,newx]:
#                             g[newy,newx,k] = g_pc[j,i,k]
#                             h[newy,newx,k] = h_pc[j,i,k]
#                         elif not struct[newy,newx]:
#                             g[j,i,BounceBackD2Q9(k)] = g_pc[j,i,k]
#                             h[j,i,BounceBackD2Q9(k)] = h_pc[j,i,k]
# =============================================================================
      
                      
# =============================================================================
#     # Streamming (when there is streaming from solid boundary nodes)
#     # source node ---> destination node
#     for i in range(nx):
#         for j in range(ny):
#             for k in range(n_pop):
#                 
#                 if struct[j,i]: # if the source node is fluid
#                     newx = i + cx[k]
#                     newy = j - cy[k]
#                     if (newx!=-1) and (newx!=nx) and (newy!=-1) and (newy!=ny) :
#                         if struct[newy,newx]: # if the destination is fluid node
#                             g[newy,newx,k] = g_pc[j,i,k]
#                             h[newy,newx,k] = h_pc[j,i,k]
#                         elif not struct[newy,newx]: # if the destination is solid node
#                             # using eq.(27)
#                             newy_fluid_node = newy - solid_normal_vectors_2d[newy,newx,0]
#                             newx_fluid_node = newx + solid_normal_vectors_2d[newy,newx,1]
#                             g[j,i,BounceBackD2Q9(k)] = g_pc[newy_fluid_node,newx_fluid_node,k]
#                             h[j,i,BounceBackD2Q9(k)] = h_pc[newy_fluid_node,newx_fluid_node,k]
#                 elif (not struct[j,i]) and (solid_normal_vectors_2d[j,i,0] != -100): # identifying a solid boundary node --> then performing collision on it using eq.(27)         
#                     y_fluid_node = j - solid_normal_vectors_2d[j,i,0]
#                     x_fluid_node = i + solid_normal_vectors_2d[j,i,1]
#                     if (x_fluid_node!=-1) and (x_fluid_node!=nx) and (y_fluid_node!=-1) and (y_fluid_node!=ny) :
#                         g_pc[j,i,BounceBackD2Q9(k)] = g_pc[y_fluid_node,x_fluid_node,k]
#                         h_pc[j,i,BounceBackD2Q9(k)] = h_pc[y_fluid_node,x_fluid_node,k]
# =============================================================================
                    
         
    # inlet velocity B.C --> just incomming populations 1,5,8:
    for j in range(1,ny-1):
#        g[j,0,1] = g[j,0,3] + 2*rho[j,0]*u1[j,0]/3 
#        h[j,0,1] = h[j,0,3] + 2*phi[j,0]*u1[j,0]/3 
#        
#        g[j,0,5] = g[j,0,7] + (g[j,0,4]-g[j,0,2])/2  + (u1[j,0]/6 + v1[j,0]/2)*rho[j,0]
#        h[j,0,5] = h[j,0,7] + (h[j,0,4]-h[j,0,2])/2  + (u1[j,0]/6 + v1[j,0]/2)*phi[j,0]
#        
#        g[j,0,8] = g[j,0,6] + (g[j,0,2]-g[j,0,4])/2 + (u1[j,0]/6 - v1[j,0]/2)*rho[j,0]
#        h[j,0,8] = h[j,0,6] + (h[j,0,2]-h[j,0,4])/2 + (u1[j,0]/6 - v1[j,0]/2)*phi[j,0]

#        g[j,0,1] = g[j,0,3] + 2*Density_l*u1[j,0]/3 
        
        
#        g[j,0,5] = g[j,0,7] + (g[j,0,4]-g[j,0,2])/2  + (u1[j,0]/6)*Density_l
        
        
#        g[j,0,8] = g[j,0,6] + (g[j,0,2]-g[j,0,4])/2 + (u1[j,0]/6)*Density_l
        

#        g[j,0,1] = g[j,0,3] + g_eq[j,0,1] - g_eq[j,0,3] #eq. 29
        h[j,0,1] = h[j,0,3] + h_eq[j,0,1] - h_eq[j,0,3]
#        
#        g[j,0,5] = g[j,0,7] + g_eq[j,0,5] - g_eq[j,0,7]
        h[j,0,5] = h[j,0,7] + h_eq[j,0,5] - h_eq[j,0,7]
#        
#        g[j,0,8] = g[j,0,6] + g_eq[j,0,8] - g_eq[j,0,6]
        h[j,0,8] = h[j,0,6] + h_eq[j,0,8] - h_eq[j,0,6]
        
    # outlet convective B.C --> just incomming populations 3,6,7
    u_convective = max(u1[1:ny-1,nx-2]) # getting convective velocity from the node before the last node --> Fakhari eq.(31) and Lou et el.(2013) eq.(14)
    for j in range(1,ny-1):
#        g[j,nx-1,3] = (g[j,nx-1,3] - u_convective*g[j,nx-2,3])/(1+u_convective)
        h[j,nx-1,3] = (h[j,nx-2,3]) #- u_convective*h[j,nx-2,3])/(1+u_convective)
        
#        g[j,nx-1,6] = (g[j,nx-1,6] - u_convective*g[j,nx-2,6])/(1+u_convective)
        h[j,nx-1,6] = (h[j,nx-2,6]) #- u_convective*h[j,nx-2,6])/(1+u_convective)
        
#        g[j,nx-1,7] = (g[j,nx-1,7] - u_convective*g[j,nx-2,6])/(1+u_convective)
        h[j,nx-1,7] = (h[j,nx-2,7]) #- u_convective*h[j,nx-2,6])/(1+u_convective)
    
    



