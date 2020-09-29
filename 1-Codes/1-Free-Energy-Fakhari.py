# Besm Allah Al-Rahman Al-Rahim
# Ref. : https://doi.org/10.1016/j.advwatres.2018.02.005
from numpy import logical_and,array, ones, zeros, empty, diag, sum, sqrt, cos, pi, mean, mod, double, matmul
from math import tanh
from tkinter import filedialog
from pandas import read_excel, DataFrame, ExcelWriter
from FakhariAdditionalFunctions import normal_vector, BounceBackD2Q9, phi_1st_2nd_derivative, phi_1st_2nd_SimplewithSolid, phi_1st_2nd_isotropic
from time import time
from copy import deepcopy

    
x = time()
# rading style
#Adress = filedialog.askopenfilename()
Adress = 'C:/Users/masoudi/Documents/GitHub/LBM-PF/1-Codes/FakhariMedia.xlsx'
excel_file = read_excel(Adress, sheet_name='new300_400' , header=None)
Style = array(excel_file.values)
Style = Style.astype(int)

# writer for debugging
#writer_g = ExcelWriter(r"D:/PythonLBM/Free-Energy/Fakhari-14-Aban/g.xlsx" , engine='xlsxwriter')
#writer_h = ExcelWriter(r"D:/PythonLBM/Free-Energy/Fakhari-14-Aban/h.xlsx" , engine='xlsxwriter')
#writer_force = ExcelWriter(r"D:/PythonLBM/Free-Energy/Fakhari-14-Aban/force.xlsx" , engine='xlsxwriter')


# general parameters
dx = 1
dt = 1
c = dx/dt 
Cs2 = double(1/3)
Cs4 = Cs2**2
ny = len(Style)
nx = len(Style[0])
n_steps = 1000

# D2Q9 Parameters
w1 = double(1/9)
w2 = double(1/36)
w3 = double(4/9)
w = array([w3,w1,w1,w1,w1,w2,w2,w2,w2] , dtype=double)
cx = array([0,1,0,-1,0,1,-1,-1,1], dtype=int)
cy = array([0,0,1,0,-1,1,1,-1,-1], dtype=int)
n_pop= len(w);

# real physical parameters
rho_real_1 = double(1001.60)                           # kg/m3         # Heavier component --> 1
rho_real_0 = double(818.55)                           # kg/m3         # lighter component --> 0
visc_real_1 = double(0.000975)                         # kg/(m.s) = rho*m2/s
visc_real_0 = double(0.000074)                     # kg/(m.s)
surface_tension_real = double(0.0294)                 # N/m = kg/s2 = rho*m3/s2
nu_real_1 = visc_real_1/rho_real_1          # m2/s
nu_real_0 = visc_real_0/rho_real_0          # m2/s
contact_angle = double(45)                          # in degrees (for user input)
contact_angle = double(contact_angle * pi / 180)    # converting to radian to use in simulation
q_in_real = double(0.0005) / 60000000                 # m3/s = ml/min / 60000000

# conversion parameters
rho_conversion = 1000            # (kg/m3)/lu  = rho           
x_conversion = double(10) * 10**(-6)     # m/lu = dx
t_conversion = double(10) * 10**(-6)     # s/lu = dt

# physical parameters in lattice units
rho_1 = rho_real_1 / rho_conversion     # Heavier component --> 1
rho_0 = rho_real_0 / rho_conversion     # lighter component --> 0
visc_1 = visc_real_1 * t_conversion / rho_conversion / x_conversion**2
visc_0 = visc_real_0 * t_conversion / rho_conversion / x_conversion**2
nu_1 = nu_real_1 * t_conversion / x_conversion**2
nu_0 = nu_real_0 * t_conversion / x_conversion**2
surface_tension = surface_tension_real * t_conversion**2 / (rho_conversion) / (x_conversion**3) 
q_in = q_in_real * t_conversion / x_conversion**3     # lu_x**3/lu_t
u_in_average = q_in / (ny-2)/dx    # lu_x/lu_t

# neccessary constants for calculation of chemical potential --> in text - above eq.(4)
interface_thickness = 3 #lu --> this value is proposed in text Sect. 3.2 and Fakhari et al. (2017a) (https://doi.org/10.1016/j.jcp.2017.03.062)
betha = 12*surface_tension/interface_thickness
kappa = 3*surface_tension*interface_thickness/2

# Phase field parameter pre-definition 
#   setting initial phase field as 0 shows --> initail present phase is the 'lighter phase'
phi = zeros((ny,nx))        # phi=0 for lighter fluid --> CO2 
                            # phi=1 for heavier fluid --> Water
phi_initialization = zeros((1,nx))
for i in range(nx): # using eq.(28)
    phi_initialization[0,i] = 0.5 * (1 + tanh(((i)+1 - interface_thickness)/(interface_thickness/2)))
phi[:,:] = deepcopy(phi_initialization)

# Identifing boundery nodes in form of 2D & 3D arrays + finding normal vector for all boundary solid nodes
    # coords_solid_boundaries[:,:,0] -->  "j"s of solid boundary nodes - if an element is equal to -100 it is not a solid boundary node
    # coords_solid_boundaries[:,:,1] -->  "i"s of solid boundary nodes - if an element is equal to -100 it is not a solid boundary node
    # neighbour_arrangement[j,i,:] --> contains all neighbours of the [j,i] nodes in direction of ":" - if all elements in "k" direction is zeros, the node is not a solid boundary node
    # solid_normal_vectors_2d[:,:,0] --> contains "j"s of all solid noundary nodes, if an element is equal to -100 it isn't a solid boundary node
    # solid_normal_vectors_2d[:,:,1] --> contains "i"s of all solid noundary nodes, if an element is equal to -100 it isn't a solid boundary node 
coords_solid_boundaries,neighbour_arrangement_3d = normal_vector().matrix_boundary_finder(Style)
solid_normal_vectors_2d = normal_vector().normal_bank_2d(neighbour_arrangement_3d)

# phase field for solid boundaries--> eq.(25) and eq.(26) 
# in matrix form
for j in range(ny):
    for i in range(nx):
        if coords_solid_boundaries[j,i,0] != -100 : # if the node is a solid boundary node
            newy = j - solid_normal_vectors_2d[j,i,0]
            newx = i + solid_normal_vectors_2d[j,i,1]
            normal_half_length = sqrt(solid_normal_vectors_2d[j,i,0]**2 + solid_normal_vectors_2d[j,i,1]**2) / 2
            if contact_angle == 0:
                #phi[y_solid_boundary[i],x_solid_boundary[i]] = phi[newy,newx]
                phi[j,i] = phi[newy,newx] #phi_s=phi_f
            else: 
                a1 = double(-1 * normal_half_length * sqrt(2*betha/kappa) * cos(contact_angle))
                phi[j,i] = 1/a1 * (1 + a1 - sqrt((1+a1)**2-4*a1*phi[newy,newx] )) -  phi[newy,newx]
                if (1+a1)**2-4*a1*phi[newy,newx] < 0:
                    print('j= {0} , i= {1} \n' .format(j,i))

# density pre-definition
rho = zeros((ny,nx))
rho = rho_0 + phi*(rho_1 - rho_0)

# viscosity pre-definition
visc = visc_0 + phi*(visc_1 - visc_0)

# velocity pre-definition
u_aparent = zeros((ny,nx))
v_aparent = zeros((ny,nx))
u1 = zeros((ny,nx))
v1 = zeros((ny,nx))
#P = ones((ny,nx)) * rho * Cs2    ############### very uncertain about this 
#P = sum(g, axis=2) + (dt/2)*(rho_1-rho_0)*Cs2*(u1*phi_x_1st_derivative + v1*phi_y_1st_derivative)   # using eq. (21b)
P = zeros((ny,nx))

# Setting velocity of inlet and outlet boundary condition
u1[1:ny-1,0] = u_in_average 
v1[1:ny-1,0] = 0
u1_old = deepcopy(u1)

###############################################################################################################################################################
# Derivatives: pre-definition and calculations ################################################################################################################
###############################################################################################################################################################
phi_x_1st_derivative = empty((ny,nx))
phi_y_1st_derivative = empty((ny,nx))
phi_x_2nd_derivative = zeros((ny,nx))    # = 'second derivatuve' of phi in x-dir
phi_y_2nd_derivative = zeros((ny,nx))    # = 'second derivatuve' of phi in y-dir
normal_x = empty((ny,nx))   # x component of normal vector to interface
normal_y = empty((ny,nx))   # y cimponent of normal vector to interface
laplacian_phi = zeros((ny,nx))
epsilon = double(10**(-32))     # in text - under eq.(2)
# =============================================================================
# # calculating first and second derivative: (simple central difference + ignoring solid boundary nodes)
# phi_x_1st_derivative,phi_y_1st_derivative,phi_x_2nd_derivative,phi_y_2nd_derivative = phi_1st_2nd_derivative(phi,Style)
# =============================================================================
# =============================================================================
# #calculating first and second derivatives (simple central difference + considering solid boundary nodes)
# phi_x_1st_derivative,phi_y_1st_derivative,phi_x_2nd_derivative,phi_y_2nd_derivative = phi_1st_2nd_SimplewithSolid(phi,Style)
# =============================================================================
# calculating first derivative of x,y and laplacian ( isotrpic centered diffrence)
phi_x_1st_derivative,phi_y_1st_derivative,laplacian_phi = phi_1st_2nd_isotropic(phi,Style)
# normal vector --> eq.(2)
magnitude_of_grad_phi = sqrt(phi_x_1st_derivative**2 + phi_y_1st_derivative**2) + epsilon
normal_x = phi_x_1st_derivative / (magnitude_of_grad_phi)
normal_y = phi_y_1st_derivative / (magnitude_of_grad_phi)
# laplacian --> laplacian = 2ndOrderDerivative_x + 2ndOrderDerivative_y
#laplacian_phi = phi_x_2nd_derivative + phi_y_2nd_derivative 
###############################################################################################################################################################
###############################################################################################################################################################
###############################################################################################################################################################


# eq.(11) - Mobility = thau_phi*Cs2*dt
mobility = double(0.01) # --> this value is proposed in text Sect. 3.2 and Fakhari et al. (2017a)
thau_phi = mobility/(Cs2*dt)

# chemical potential & free-energy pre-definition
free_energy = (phi**2) * (1-phi)**2     # eq.(4)
chemical_potential = 4*betha*phi*(phi-1)*(phi-0.5) - kappa*laplacian_phi  + 0.  #eq.(5)



# force pre-definition
Fs_x = zeros((ny,nx), dtype=double)
Fs_y = zeros((ny,nx), dtype=double)
force = zeros((ny,nx,n_pop) , dtype=double)
# Surface tension (F_s) calculation --> in text after eq.(7)
Fs_x = chemical_potential * phi_x_1st_derivative 
Fs_y = chemical_potential * phi_y_1st_derivative 
    
# initialization of distribution function
g_eq = zeros((ny,nx,n_pop), dtype=double)
h_eq = zeros((ny,nx,n_pop), dtype=double)
g_pc = zeros((ny,nx,n_pop), dtype=double )
h_pc = zeros((ny,nx,n_pop), dtype=double) # phase field distribuition function
collision_operator = zeros((ny,nx,n_pop), dtype= double)

# equilibrium distribuion
for k in range(n_pop):
    # eq. 10 of fakhari:
    Gamma = w[k] * (1 + (cx[k]*u1 + cy[k]*v1)/Cs2 + (cx[k]*u1 + cy[k]*v1)**2/(2*Cs4) - (u1*u1+v1*v1)/(2*Cs2))
    
    # eq. 9 of fakhari :
    h_eq[:,:,k] = Gamma*phi + w[k] * mobility/Cs2 * (4/interface_thickness*phi*(1-phi)) * (cx[k]*normal_x + cy[k]*normal_y)
    
    # force calculation on each lattice node --> eq.(15) :
    force[:,:,k] = ( (Gamma-w[k])*(rho_1-rho_0)*Cs2 + Gamma*chemical_potential ) * ((cx[k]-u1)*phi_x_1st_derivative + (cy[k]-v1)*phi_y_1st_derivative)
    # g_eq --> combining eq.(18) into eq.(17) :
    g_eq[:,:,k] = ( P*w[k] + rho*Cs2*(Gamma-w[k]) ) - 0.5 * force[:,:,k]
    
    # output for debugging:
    #df_g = DataFrame(g_eq[:,:,k])
    #df_h = DataFrame(h_eq[:,:,k])
    #df_force = DataFrame(force[:,:,k])
    
    #df_g.to_excel(writer_g, sheet_name=str(k))
    #df_h.to_excel(writer_h, sheet_name=str(k))  
    #df_force.to_excel(writer_force, sheet_name=str(k))

#writer_g.save()
#writer_g.close()
#writer_h.save()
#writer_h.close()
#writer_force.save()
#writer_force.close()

    
g = deepcopy(g_eq)
h = deepcopy(h_eq)
g_pc = deepcopy(g_eq)
h_pc = deepcopy(h_eq)

phi = sum(h, axis=2)

# u calculation --> #eq.(21a)
g_right = g[:,:,[1,5,8]]
g_left = g[:,:,[3,6,7]]
rho_u = sum(g_right, axis=2) - sum(g_left, axis=2)
u_aparent = rho_u/rho
u1 = u_aparent/Cs2 + dt*Fs_x/(2*rho)     # eq.(21a)
    
# v calculation --> #eq.(21a)
g_up = g[:,:,[2,5,6]]
g_down = g[:,:,[4,7,8]]
rho_v = sum(g_up, axis=2) - sum(g_down, axis=2) 
v_aparent = rho_v/rho  
v1 = v_aparent/Cs2 + dt*Fs_y/(2*rho)     #eq.(21a)


# MRT parameters --> from another paper from Fakhari, 2013 --> "Multiple-relaxation-time lattice Boltzmann method for immiscible fluids at high Reynolds numbers"
M = array([[ 1,  1,  1,  1,  1,  1,  1,  1,  1],
           [-4, -1, -1, -1, -1,  2,  2,  2,  2],
           [ 4, -2, -2, -2, -2,  1,  1,  1,  1],
           [ 0,  1,  0, -1,  0,  1, -1, -1,  1],
           [ 0, -2,  0,  2,  0,  1, -1, -1,  1],
           [ 0,  0,  1,  0, -1,  1,  1, -1, -1],
           [ 0,  0, -2,  0,  2,  1,  1, -1, -1],
           [ 0,  1, -1,  1, -1,  0,  0,  0,  0],
           [ 0,  0,  0,  0,  0,  1, -1,  1, -1]] , dtype=double)
M_inverse = array([[4,   -4,   4,   0,	 0,	 0,	 0,	 0,	 0],
                   [4,   -1,  -2,   6,  -6,	 0,	 0,	 9,	 0],
                   [4,   -1,  -2,   0,	 0,	 6,	-6,	-9,	 0],
                   [4,   -1,  -2,  -6,	 6,	 0,	 0,	 9,	 0],
                   [4,   -1,  -2,   0,	 0,	-6,	 6,	-9,	 0],
                   [4,    2,   1,   6,   3,	 6,	 3,	 0,	 9],
                   [4,    2,   1,  -6,  -3,	 6,	 3,	 0,	-9],
                   [4,    2,   1,  -6,  -3,	-6,	-3,	 0,	 9],
                   [4,    2,   1,   6,	 3,	-6,	-3,	 0,	-9]] , dtype=double) / 36 
  # parameter           # related to:                   # role of parameter in LB:
# from Fakhari and Bolster 
omega0_rho = 0          # density               -->     physical
omega1_e = 1            # energy                -->     tunning
omega2_epsilon = 1      # enegy**2              -->     tunning
omega3_jx = 0           # mass flux ()momentum  -->     physical
omega4_qx = 1.7           # energy flux           -->     tunning
omega5_jy = 0           # mass flux ()momentum  -->     physical
omega6_qx = 1.7           # energy flux           -->     tunning
omega7_pxx = 1          # digonal stress tensor -->     tunning
omega8_pxy = 1          # off-digonal stress tensor --> tunning
S_hat = array([[ 0.,  0,   0,   0,   0,    0,  0,    0,  0],
               [ 0,   1.,  0,   0,   0,    0,  0,    0,  0],
               [ 0,   0,   1.,  0,   0,    0,  0,    0,  0],
               [ 0,   0,   0,   0.,  0,    0,  0,    0,  0],
               [ 0,   0,   0,   0,   1.7,  0,  0,    0,  0],
               [ 0,   0,   0,   0,   0,    0., 0,    0,  0],
               [ 0,   0,   0,   0,   0,    0,  1.7,  0,  0],
               [ 0,   0,   0,   0,   0,    0,  0,    1., 0],
               [ 0,   0,   0,   0,   0,    0,  0,    0,  1.]] , dtype=double)
# S_hat = diag((omega0_rho,omega1_e,omega2_epsilon,omega3_jx,omega4_qx,omega5_jy,omega6_qx,omega7_pxx,omega8_pxy))
minus_M_S_M = (-1) * matmul( matmul(M_inverse, S_hat) , M )    # = -M_inverse*S*M

#writer = ExcelWriter(r"D:/PythonLBM/Free-Energy/force.xlsx" , engine='xlsxwriter')
# main loop
for t in range(n_steps):
    print("t= {0}" .format(t))
    if mod(t,5)==0:
        zzz=1
    
    # phase field calculation --> eq.(12)
    phi = sum(h, axis=2)
    phi[:,0] = 0        # in text- above eq.(30) zero phase field at inlet
#    phi[logical_and (Style==0,(coords_solid_boundaries[:,:,0]==-100)) ]=0
          
    # phase field for solid boundaries--> eq.(25) and eq.(26) 
    # in matrix form
    for j in range(ny):
        for i in range(nx):
            if coords_solid_boundaries[j,i,0] != -100 :
                newy = j - solid_normal_vectors_2d[j,i,0]
                newx = i + solid_normal_vectors_2d[j,i,1]
                normal_half_length = sqrt(solid_normal_vectors_2d[j,i,0]**2 + solid_normal_vectors_2d[j,i,1]**2) / 2
                if contact_angle == 0:
                    #phi[y_solid_boundary[i],x_solid_boundary[i]] = phi[newy,newx]
                    phi[j,i] = phi[newy,newx]
                else: 
                    a1 = double(-1 * normal_half_length * sqrt(2*betha/kappa) * cos(contact_angle))
                    phi[j,i] = 1/a1 * (1 + a1 - sqrt((1+a1)**2-4*a1*phi[newy,newx] )) -  phi[newy,newx] 
                    if (1+a1)**2-4*a1*phi[newy,newx] < 0:
                        print('j= {0} , i= {1} , t= {2}\n' .format(j,i,t))
    
    # Density calculation --> eq.(13) which was incorrect in paper --> 
    rho = rho_0 + phi*(rho_1 - rho_0)

# =============================================================================
#     # calculating first and second derivative: (simple central difference + ignoring solid boundary nodes)
#     phi_x_1st_derivative,phi_y_1st_derivative,phi_x_2nd_derivative,phi_y_2nd_derivative = phi_1st_2nd_derivative(phi,Style)
# =============================================================================
# =============================================================================
#     #calculating first and second derivatives (simple central difference + considering solid boundary nodes)
#     phi_x_1st_derivative,phi_y_1st_derivative,phi_x_2nd_derivative,phi_y_2nd_derivative = phi_1st_2nd_SimplewithSolid(phi,Style)
#     
# =============================================================================
    # calculating first derivative of x,y and laplacian ( isotrpic centered diffrence)
    phi_x_1st_derivative,phi_y_1st_derivative,laplacian_phi = phi_1st_2nd_isotropic(phi,Style)

    # normal vector --> eq.(2)
    magnitude_of_grad_phi = sqrt(phi_x_1st_derivative**2 + phi_y_1st_derivative**2) + epsilon
    normal_x = phi_x_1st_derivative / (magnitude_of_grad_phi )
    normal_y = phi_y_1st_derivative / (magnitude_of_grad_phi )
    
# =============================================================================
#     # laplacian --> laplacian = 2ndOrderDerivative_x + 2ndOrderDerivative_y
#     laplacian_phi = phi_x_2nd_derivative + phi_y_2nd_derivative
# =============================================================================
    
    # chemiucal potential calculation --> eq.(5)
    free_energy = (phi**2) * (1-phi)**2     # eq.(4)
    chemical_potential = 4*betha*phi*(phi-1)*(phi-1/2) - kappa*laplacian_phi    #eq.(5)

    
    # Surface tension (F_s) calculation --> in text after eq.(7)
    Fs_x = chemical_potential * phi_x_1st_derivative 
    Fs_y = chemical_potential * phi_y_1st_derivative 
    
    # u calculation --> #eq.(21a)
    g_right = g[:,:,[1,5,8]]
    g_left = g[:,:,[3,6,7]]
    rho_u = sum(g_right, axis=2) - sum(g_left, axis=2)
    u_aparent = rho_u/rho 
    u1 = u_aparent/Cs2 + dt*Fs_x/(2*rho)     # eq.(21a)
    
    # v calculation --> #eq.(21a)
    g_up = g[:,:,[2,5,6]]
    g_down = g[:,:,[4,7,8]]
    rho_v = sum(g_up, axis=2) - sum(g_down, axis=2) 
    v_aparent = rho_v/rho  
    v1 = v_aparent/Cs2 + dt*Fs_y/(2*rho)     #eq.(21a)
    
    #u1[Style==0]=0
    #v1[Style==0]=0
    # Setting velocity of inlet and outlet boundary condition
    u1[1:ny-1,0] = u_in_average 
    v1[1:ny-1,0] = 0
    
    if (mod(t,50)==0)&(t!=0):
        xx = time()
        conv= (mean(sum(u1,axis=0))-mean(sum(u1_old,axis=0)))/mean(sum(u1,axis=0)) ;
        print('Conv= {0} \t\t runtime= {1}' .format(conv,xx-x))
        # Checking convergence
        if abs(conv)<10**(-10):
            break;
        else:
            u1_old=u1
            
    # pressure calculation --> eq.(21b)
    P = sum(g, axis=2) + (dt/2)*(rho_1-rho_0)*Cs2*(u1*phi_x_1st_derivative + v1*phi_y_1st_derivative)
    
    # viscosity calculation --> eq.(19)
    visc = visc_0 + phi*(visc_1 - visc_0)
    
    ################################################################
    # relaxation time --> eq.(20)
    thau_rho = visc / (rho*Cs2)         # what is this used for?
    ################################################################
    Sv = 1/(thau_rho + 0.5)     # from Fakhari and Bolster 2017 - eq.(23)
    
    # equilibrium distribuion
    for k in range(n_pop):
        # eq. 10 of fakhari:
        Gamma = w[k] * (1 + (cx[k]*u1 + cy[k]*v1)/Cs2 + (cx[k]*u1 + cy[k]*v1)**2/(2*Cs4) - (u1*u1+v1*v1)/(2*Cs2))
        
        # eq. 9 of fakhari: 
        h_eq[:,:,k] = Gamma*phi + w[k] * mobility/Cs2 * (4/interface_thickness*phi*(1-phi)) * (cx[k]*normal_x + cy[k]*normal_y)
        
        # force calculation on each lattice node --> eq.(15)
        force[:,:,k] = dt * ( (Gamma-w[k])*(rho_1-rho_0)*Cs2 + Gamma*chemical_potential ) * ((cx[k]-u1)*phi_x_1st_derivative + (cy[k]-v1)*phi_y_1st_derivative)
        #df = DataFrame(force[:,:,k])
        #df.to_excel(writer, sheet_name='k= '+ str(k))
        
        # g_eq --> combining eq.(18) into eq.(17)
        g_eq[:,:,k] = ( P*w[k] + rho*Cs2*(Gamma-w[k]) ) - 0.5*force[:,:,k]
        
    #writer.save()
    #writer.close
    
    # MRT collision operator for Hydrodynamic --> eq.(16) in another algorithmic form
    for i in range(nx):
        for j in range(ny):
            if Style[j,i]:
                S_hat[7,7] = Sv[j,i]  # fakhari and bolster eq.(22)
                S_hat[8,8] = Sv[j,i]  # fakhari and bolster eq.(22)
                minus_M_S_M = (-1) * matmul( matmul(M_inverse, S_hat) , M )
                #collision_operator[j,i,:] = (-1) * M_inverse.dot(S_hat.dot(M.dot(g[j,i,:]-g_eq[j,i,:])))
                collision_operator[j,i,:] = matmul(minus_M_S_M , (g[j,i,:]-g_eq[j,i,:]))
                #collision_operator[j,i,:] = minus_M_S_M * (g[j,i,:]-g_eq[j,i,:])
    
# =============================================================================
#     # Single Relaxation-time
#     for j in range(ny):
#         for i in range(nx):
#             collision_operator[j,i,:] = - (g[j,i,:]-g_eq[j,i,:])/(thau_rho[j,i]+0.5)
# =============================================================================
     
    
    # collision step for phase field --> eq.(8)
    h_pc = h - (h-h_eq)/(thau_phi+0.5)
    # collision step for hydrodynamics --> eq.(14)
    g_pc = g + collision_operator + force

    # replacing the distribution of solid boundary nodes
    for j in range(ny):
        for i in range(nx):
            if coords_solid_boundaries[j,i,0]!=-100 : # if this is a solid boundari node
                for k in range(n_pop):
                    newy = j - solid_normal_vectors_2d[j,i,0]
                    newx = i + solid_normal_vectors_2d[j,i,1]
                    
                    g_pc[j,i,k] = g_pc[newy,newx,BounceBackD2Q9(k)]
                    h_pc[j,i,k] = h_pc[newy,newx,BounceBackD2Q9(k)]

                    #newy = j - cy[k]
                    #newx = i + cx[k]
                    #if (newx!=-1) and (newx!=nx) and (newy!=-1) and (newy!=ny) :
                        #g_pc[j,i,k] = g_pc[newy,newx,BounceBackD2Q9(k)]
                        #h_pc[j,i,k] = h_pc[newy,newx,BounceBackD2Q9(k)]

    # streamimng the fluid and solid boundary nodes 
    for j in range(ny):
        for i in range(nx):
            if Style[j,i] or (coords_solid_boundaries[j,i,0]!=-100) :  # if this is a fluid node or a solid boundary node 
                for k in range(n_pop):
                    newy = j - cy[k] ;
                    newx = i + cx[k] ;
                    if (newx!=-1) and (newx!=nx) and (newy!=-1) and (newy!=ny) :
                        g[newy,newx,k] = g_pc[j,i,k] ;
                        h[newy,newx,k] = h_pc[j,i,k] ;
        
# =============================================================================
#     #Streaming (when there is no streamming from solid boundary nodes)
#     for i in range(nx):
#         for j in range(ny):
#             for k in range(n_pop):
#                 if Style[j,i]:
#                     newx = i + cx[k]
#                     newy = j - cy[k]
#                     if (newx!=-1) and (newx!=nx) and (newy!=-1) and (newy!=ny) :
#                         if Style[newy,newx]:
#                             g[newy,newx,k] = g_pc[j,i,k]
#                             h[newy,newx,k] = h_pc[j,i,k]
#                         elif not Style[newy,newx]:
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
#                 if Style[j,i]: # if the source node is fluid
#                     newx = i + cx[k]
#                     newy = j - cy[k]
#                     if (newx!=-1) and (newx!=nx) and (newy!=-1) and (newy!=ny) :
#                         if Style[newy,newx]: # if the destination is fluid node
#                             g[newy,newx,k] = g_pc[j,i,k]
#                             h[newy,newx,k] = h_pc[j,i,k]
#                         elif not Style[newy,newx]: # if the destination is solid node
#                             # using eq.(27)
#                             newy_fluid_node = newy - solid_normal_vectors_2d[newy,newx,0]
#                             newx_fluid_node = newx + solid_normal_vectors_2d[newy,newx,1]
#                             g[j,i,BounceBackD2Q9(k)] = g_pc[newy_fluid_node,newx_fluid_node,k]
#                             h[j,i,BounceBackD2Q9(k)] = h_pc[newy_fluid_node,newx_fluid_node,k]
#                 elif (not Style[j,i]) and (solid_normal_vectors_2d[j,i,0] != -100): # identifying a solid boundary node --> then performing collision on it using eq.(27)         
#                     y_fluid_node = j - solid_normal_vectors_2d[j,i,0]
#                     x_fluid_node = i + solid_normal_vectors_2d[j,i,1]
#                     if (x_fluid_node!=-1) and (x_fluid_node!=nx) and (y_fluid_node!=-1) and (y_fluid_node!=ny) :
#                         g_pc[j,i,BounceBackD2Q9(k)] = g_pc[y_fluid_node,x_fluid_node,k]
#                         h_pc[j,i,BounceBackD2Q9(k)] = h_pc[y_fluid_node,x_fluid_node,k]
# =============================================================================
                    
         
    # inlet velocity B.C --> just incomming populations 1,5,8:
    for j in range(1,ny-1):
        g[j,0,1] = g[j,0,3] + 2*rho[j,0]*u1[j,0]/3 
        h[j,0,1] = h[j,0,3] + 2*phi[j,0]*u1[j,0]/3 
        
        g[j,0,5] = g[j,0,7] + (g[j,0,4]-g[j,0,2])/2  + (u1[j,0]/6 + v1[j,0]/2)*rho[j,0]
        h[j,0,5] = h[j,0,7] + (h[j,0,4]-h[j,0,2])/2  + (u1[j,0]/6 + v1[j,0]/2)*phi[j,0]
        
        g[j,0,8] = g[j,0,6] + (g[j,0,2]-g[j,0,4])/2 + (u1[j,0]/6 - v1[j,0]/2)*rho[j,0]
        h[j,0,8] = h[j,0,6] + (h[j,0,2]-h[j,0,4])/2 + (u1[j,0]/6 - v1[j,0]/2)*phi[j,0]
        
    # outlet convective B.C --> just incomming populations 3,6,7
    u_convective = max(u1[1:ny-1,nx-2]) # getting convective velocity from the node before the last node --> Fakhari eq.(31) and Lou et el.(2013) eq.(14)
    for j in range(1,ny-1):
        g[j,nx-1,3] = (g[j,nx-1,3] - u_convective*g[j,nx-2,3])/(1+u_convective)
        h[j,nx-1,3] = (h[j,nx-1,3] - u_convective*h[j,nx-2,3])/(1+u_convective)
        
        g[j,nx-1,6] = (g[j,nx-1,6] - u_convective*g[j,nx-2,6])/(1+u_convective)
        h[j,nx-1,6] = (h[j,nx-1,6] - u_convective*h[j,nx-2,6])/(1+u_convective)
        
        g[j,nx-1,7] = (g[j,nx-1,7] - u_convective*g[j,nx-2,6])/(1+u_convective)
        h[j,nx-1,7] = (h[j,nx-1,7] - u_convective*h[j,nx-2,6])/(1+u_convective)
    
    
                        
    

















