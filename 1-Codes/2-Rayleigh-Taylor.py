# Ref.: 
# 1) https://doi.org/10.1016/j.advwatres.2018.02.005 : 
# 2) https://doi.org/10.1103/PhysRevE.96.053301
from numpy import identity,array, ones, zeros, sum, sqrt, cos, pi, mean, mod, double, matmul ,diag
#from math import tanh
#from tkinter import filedialog
#from pandas import read_excel, DataFrame, ExcelWriter
from AdditionalFunctions import normal_vector, BounceBackD2Q9, phi_1st_2nd_isotropic#phi_1st_2nd_derivative, phi_1st_2nd_SimplewithSolid, 
from time import time
from copy import deepcopy
import pickle
x = time()

# 1: inputs ##################################################################
## 1.1 Geometry:  (y,x)
# 0 ---> Solid node      1---> fluid node
# the code is written in a way that y=0 is the top boundary and x=0 is the left boundary

# rading the structure from an external source: 
#Adress = filedialog.askopenfilename()
#Adress = 'C:/Users/masoudi/Documents/GitHub/LBM-PF/1-Codes/FakhariMedia.xlsx'
#excel_file = read_excel(Adress, sheet_name='new300_400' , header=None)
#struct = array(excel_file.values)

# Poiseuille flow :
struct= ones ([64,30]) 
struct[(0,63),:]=0


struct = struct.astype(int) # to make sure the data type is ok
ny = len(struct)                # domain size along y axis
nx = len(struct[0])             # domain size along x axis

## 1.2 real physical parameters:
dx_ph = 10e-6                                    # grid resolution(m) 
dt_ph = 10e-6                                    # physical time step (m)
Density_real_H = double(1000)                 # kg/m3 # Heavier component
Density_real_L = double(100)                  # kg/m3 # lighter component
visc_real_H = double(0.005)                   # kg/(m.s) = rho*m2/s
visc_real_L = double(0.000005)                   # kg/(m.s)
surface_tension_real = double(0)            # N/m = kg/s2 = rho*m3/s2
contact_angle = double(90)                       # in degrees (for user input)
#q_in_real = double(0.05) / 60000000              # injction rate m3/s = ml/min / 60000000
g_b_x=1e-6                                         # body force in x direction
g_b_y=0.0
## 1.3 Model related or LB parameters:
n_steps = 10000                                 # time steps
PhaseField_l = double(0)                                 # for lighter fluid --> CO2 
PhaseField_h = double(1)                                 # for heavier fluid --> Water
PhaseField_av=(PhaseField_h+PhaseField_l)/2  
Density_h = double(1)                            # Heavier component LB density
interface_thickness =double(4)                         #lu --> this value is proposed in text Sect. 3.2 and Ref.2
mobility = double(0.01)                               # eq.(11) - Mobility = thau_phi*Cs2*dt # --> this value is proposed in text Sect. 3.2 and Fakhari et al. (2017a)

## 1.4 initial condition:



#phi
#Poiseuille:
phi = zeros((ny,nx))                                   
AA=int((((ny-2)/2)+1))
phi [AA:ny-1,:]=1

#phi_initialization = zeros((1,nx))
#for i in range(nx): # using eq.(28)
#    phi_initialization[0,i] = 0.5 * (1 + tanh(((i)+1 - interface_thickness)/(interface_thickness/2)))
#phi[:,:] = deepcopy(phi_initialization)

#P
P = zeros((ny,nx))
# velocity
u1 = zeros((ny,nx))
v1 = zeros((ny,nx))

u1_old = deepcopy(u1)
#2: functions ################################################################
def getPressure(g):
    # P here is normalized pressure ?? P* = P / (rho*Cs2)
    P = sum(g, axis=2) + (dt/2)*(Density_h-Density_l)*Cs2*(u1*phi_x_1st_derivative + v1*phi_y_1st_derivative)   # using eq. (21b)
    # P = sum(g, axis=2) # eq 32a ref. 2
    
    return P                        
    
def update_force(P,rho,phi,visc,u1,v1,g,phi_x_1st_derivative,phi_y_1st_derivative,laplacian_phi):
    for k in range(n_pop):
        # eq. 10 of Ref.2:
        Gamma = w[k] * (1 + (cx[k]*u1 + cy[k]*v1)/Cs2 + (cx[k]*u1 + cy[k]*v1)**2/(2*Cs4) - (u1*u1+v1*v1)/(2*Cs2))
        # g_eq -->  eq.(17) Ref. 2:
        g_eq[:,:,k] = P*w[k] + (Gamma-w[k])
    
    # force terms Eq.18 Ref.2:
    # Pressure Force Eq. 19 Ref.2:
    F_P_x = - P * Cs2 * ((Density_h-Density_l)/(PhaseField_h-PhaseField_l)) * phi_x_1st_derivative
    F_P_y = - P * Cs2 * ((Density_h-Density_l)/(PhaseField_h-PhaseField_l)) * phi_y_1st_derivative
    # chemical potential & free-energy pre-definition
    #free_energy = (phi**2) * (1-phi)**2     # eq.(4)
    chemical_potential = 4*betha*(phi-PhaseField_l)*(phi-PhaseField_h)*(phi-PhaseField_av) - kappa*laplacian_phi  + 0.  #eq.(5)
    # Surface tension (F_s) calculation --> in text after eq.(7)
    Fs_x = chemical_potential * phi_x_1st_derivative 
    Fs_y = chemical_potential * phi_y_1st_derivative 
    #Body Force:
    Fb_x= rho * g_b_x # -1.0*(rho-Density_l)*BuoyancyX + rho*g_b_x
    Fb_y= rho * g_b_y
    # viscous force
    
    ################################################################
    # relaxation time for hydrodynamic calculations --> eq.(25) Ref.2
    tau_rho = visc / (rho*Cs2)
    ################################################################
    Sv = 1/(tau_rho + 0.5)     # eq.(29) Ref.2

    # MRT collision operator for Hydrodynamic --> eq.(31) Ref.2
    for i in range(nx):
        for j in range(ny):
            if struct[j,i]:
                S_hat[7,7] = Sv[j,i]  #  eq.(28) Ref.2
                S_hat[8,8] = Sv[j,i]  #  eq.(28) Ref.2
                minus_M_S_M = (-1) * matmul( matmul(M_inverse, S_hat) , M )                
                collision_operator_F_mu[j,i,:] = matmul(minus_M_S_M , (g[j,i,:]-g_eq[j,i,:]))
    
    Sxx =  collision_operator_F_mu * cx * cx
    Sxy =  collision_operator_F_mu * cx * cy
    Syy =  collision_operator_F_mu * cy * cy
    sxx = sum(Sxx, axis=2)
    sxy = sum(Sxy, axis=2)
    syy = sum(Syy, axis=2)
    
    Fmu_x = (-tau_rho) * (sxx*phi_x_1st_derivative + sxy*phi_y_1st_derivative) * ((Density_h-Density_l)/(PhaseField_h-PhaseField_l)) # eq.31 Ref.2 based on Mitchel code TCLB
    Fmu_y = (-tau_rho) * (sxy*phi_x_1st_derivative + syy*phi_y_1st_derivative) * ((Density_h-Density_l)/(PhaseField_h-PhaseField_l)) # eq.31 Ref.2 based on Mitchel code TCLB 
    #Total Hydrodynamic Force
    Ft_x = Fmu_x + Fb_x + Fs_x + F_P_x
    Ft_y = Fmu_y + Fb_y + Fs_y + F_P_y
    
    return Ft_x, Ft_y


#3: Parameters (no need to change)
# general parameters
dx = 1                          # LB length scale
dt = 1                          # LB time scale
c = dx/dt                       # dx/dt
Cs2 = double(1/3)               # sound velocity square Cs^2
Cs4 = Cs2**2
tau_phi = mobility/(Cs2*dt)

# D2Q9 Parameters
w1 = double(1/9)
w2 = double(1/36)
w3 = double(4/9)
w = array([w3,w1,w1,w1,w1,w2,w2,w2,w2] , dtype=double) # velocity sets
cx = array([0,1,0,-1,0,1,-1,-1,1], dtype=int)          #weight coefficients 
cy = array([0,0,1,0,-1,1,1,-1,-1], dtype=int)          #weight coefficients 
n_pop= len(w);                                         # number of populations

# conversion from physical to LB: 
C_rho = Density_real_H/Density_h                       # density convertion factor: C_rho= rho/rho*          
C_L = dx_ph/dx                                         # Length convertion factor: C_L= L/L* eq. 7.2 Timm KrUger book
C_t = dt_ph/dt                                         # time convertion factor: C_t= t/t* eq. 7.12 Timm KrUger book

# other parameters:
nu_real_1 = visc_real_H/Density_real_H                 # kinematic viscosity m2/s
nu_real_0 = visc_real_L/Density_real_L                 # kinematic viscosity m2/s
contact_angle = double(contact_angle * pi / 180)       # converting to radian to use in simulation
Density_l = Density_real_L / C_rho                     # lighter component density
visc_H = visc_real_H * C_t / C_rho / C_L**2
visc_L = visc_real_L * C_t / C_rho / C_L**2
nu_1 = nu_real_1 * C_t / C_L**2
nu_0 = nu_real_0 * C_t / C_L**2
surface_tension = surface_tension_real * C_t**2 / (C_rho) / (C_L**3) 
#q_in = q_in_real * C_t / C_L**3     # lu_x**3/lu_t

# neccessary constants for calculation of chemical potential --> in text - above eq.(4)
betha = 12*surface_tension/interface_thickness
kappa = 3*surface_tension*interface_thickness/2
#dimentionless group



 
# MRT parameters --> from another paper from Fakhari, 2013 --> "Multiple-relaxation-time lattice Boltzmann method for immiscible fluids at high Reynolds numbers"
M = array([[ 1,  1,  1,  1,  1,  1,  1,  1,  1],
           [-4, -1, -1, -1, -1,  2,  2,  2,  2],
           [ 4, -2, -2, -2, -2,  1,  1,  1,  1],
           [ 0,  1,  0, -1,  0,  1, -1, -1,  1],
           [ 0, -2,  0,  2,  0,  1, -1, -1,  1],
           [ 0,  0,  1,  0, -1,  1,  1, -1, -1],
           [ 0,  0, -2,  0,  2,  1,  1, -1, -1],
           [ 0,  1, -1,  1, -1,  0,  0,  0,  0],
           [ 0,  0,  0,  0,  0,  1, -1,  1, -1]] , dtype=double) # LB book 10.30
M_inverse = array([[4,   -4,   4,   0,	 0,	 0,	 0,	 0,	 0],
                   [4,   -1,  -2,   6,  -6,	 0,	 0,	 9,	 0],
                   [4,   -1,  -2,   0,	 0,	 6,	-6,	-9,	 0],
                   [4,   -1,  -2,  -6,	 6,	 0,	 0,	 9,	 0],
                   [4,   -1,  -2,   0,	 0,	-6,	 6,	-9,	 0],
                   [4,    2,   1,   6,   3,	 6,	 3,	 0,	 9],
                   [4,    2,   1,  -6,  -3,	 6,	 3,	 0,	-9],
                   [4,    2,   1,  -6,  -3,	-6,	-3,	 0,	 9],
                   [4,    2,   1,   6,	 3,	-6,	-3,	 0,	-9]] , dtype=double) / 36  ## LB book 10.33

# parameter             # related to:                   # role of parameter in LB:
# from Fakhari and Bolster 
#omega0_rho = 0          # density               -->     physical
#omega1_e = 1            # energy                -->     tunning
#omega2_epsilon = 1      # enegy**2              -->     tunning
#omega3_jx = 0           # mass flux ()momentum  -->     physical
#omega4_qx = 1.7           # energy flux           -->     tunning
#omega5_jy = 0           # mass flux ()momentum  -->     physical
#omega6_qx = 1.7           # energy flux           -->     tunning
#omega7_pxx = 1          # digonal stress tensor -->     tunning
#omega8_pxy = 1          # off-digonal stress tensor --> tunning
# S_hat = diag((omega0_rho,omega1_e,omega2_epsilon,omega3_jx,omega4_qx,omega5_jy,omega6_qx,omega7_pxx,omega8_pxy))
#S_hat = identity(9, dtype=double) # it will be updated to eq. 28 Ref. 2. later
S_hat = diag((0,1,1,0,1.7,0,1.7,1,1))
#normal vectors: should work on them but its ok for now
# Identifing boundery nodes in form of 2D & 3D arrays + finding normal vector for all boundary solid nodes
    # coords_solid_boundaries[:,:,0] -->  "j"s of solid boundary nodes - if an element is equal to -100 it is not a solid boundary node
    # coords_solid_boundaries[:,:,1] -->  "i"s of solid boundary nodes - if an element is equal to -100 it is not a solid boundary node
    # neighbour_arrangement[j,i,:] --> contains all neighbours of the [j,i] nodes in direction of ":" - if all elements in "k" direction is zeros, the node is not a solid boundary node
    # solid_normal_vectors_2d[:,:,0] --> contains "j"s of all solid noundary nodes, if an element is equal to -100 it isn't a solid boundary node
    # solid_normal_vectors_2d[:,:,1] --> contains "i"s of all solid noundary nodes, if an element is equal to -100 it isn't a solid boundary node 
coords_solid_boundaries,neighbour_arrangement_3d = normal_vector().matrix_boundary_finder(struct)
solid_normal_vectors_2d = normal_vector().normal_bank_2d(neighbour_arrangement_3d)

# phase field for solid boundaries--> eq.(25) and eq.(26) 
# in matrix form
for j in range(ny):
    for i in range(nx):
        if coords_solid_boundaries[j,i,0] != -100 : # if the node is a solid boundary node
            newy = j - solid_normal_vectors_2d[j,i,0]
            newx = i + solid_normal_vectors_2d[j,i,1]
            normal_half_length = sqrt(solid_normal_vectors_2d[j,i,0]**2 + solid_normal_vectors_2d[j,i,1]**2) / 2
            if contact_angle == pi / 2:
                phi[j,i] = phi[newy,newx] #phi_s=phi_f
            else: 
                a1 = double(-1 * normal_half_length * sqrt(2*betha/kappa) * cos(contact_angle))
                phi[j,i] = 1/a1 * (1 + a1 - sqrt((1+a1)**2-4*a1*phi[newy,newx] )) -  phi[newy,newx]
                if (1+a1)**2-4*a1*phi[newy,newx] < 0:
                    print('j= {0} , i= {1} \n' .format(j,i))


# density pre-definition
rho = Density_l + (Density_h - Density_l) * (phi - PhaseField_l)/(PhaseField_h - PhaseField_l)
# viscosity pre-definition
visc = visc_L + (visc_H - visc_L) * (phi - PhaseField_l)/(PhaseField_h - PhaseField_l)

# calculating first derivative of x,y and laplacian ( isotrpic centered diffrence)
phi_x_1st_derivative,phi_y_1st_derivative,laplacian_phi = phi_1st_2nd_isotropic(phi,struct)
# normal vector --> eq.(2)
epsilon = double(10**(-32))     # in text - under eq.(2)
magnitude_of_grad_phi = sqrt(phi_x_1st_derivative**2 + phi_y_1st_derivative**2) + epsilon
normal_x = phi_x_1st_derivative / (magnitude_of_grad_phi)
normal_y = phi_y_1st_derivative / (magnitude_of_grad_phi)


# initialization of distribution function
g_eq = zeros((ny,nx,n_pop), dtype=double)
h_eq = zeros((ny,nx,n_pop), dtype=double)
F_phi= zeros((ny,nx,n_pop), dtype=double)

# equilibrium distribuion
for k in range(n_pop):
    # eq. 10 of Ref.2:
    Gamma = w[k] * (1 + (cx[k]*u1 + cy[k]*v1)/Cs2 + (cx[k]*u1 + cy[k]*v1)**2/(2*Cs4) - (u1*u1+v1*v1)/(2*Cs2))   
    # eq. 7 of Ref.2:
    F_phi[:,:,k]  = ( 1 - 4 *(phi-PhaseField_av)**2)/interface_thickness * w[k] *(cx[k]*normal_x + cy[k]*normal_y)
    # eq. 9 of Ref.2:
    h_eq[:,:,k] = Gamma*phi - 0.5 * F_phi[:,:,k]  
    # g_eq -->  eq.(17) Ref. 2:
    g_eq[:,:,k] = P*w[k] + (Gamma-w[k]) # - 0.5 * force[:,:,k]

g = deepcopy(g_eq)
h = deepcopy(h_eq)

collision_operator = zeros((ny,nx,n_pop), dtype= double)
collision_operator_F_mu = zeros((ny,nx,n_pop), dtype= double)
force = zeros((ny,nx,n_pop), dtype= double)

# update the force terms:
Ft_x, Ft_y = update_force(P,rho,phi,visc,u1,v1,g,phi_x_1st_derivative,phi_y_1st_derivative,laplacian_phi)
    
# main loop

for t in range(n_steps):
    #print("t= {0}" .format(t))
    
    ##collision: #######################################################
    
    # equilibrium distribuion
    for k in range(n_pop):
        # eq. 10 of Ref.2:
        Gamma = w[k] * (1 + (cx[k]*u1 + cy[k]*v1)/Cs2 + (cx[k]*u1 + cy[k]*v1)**2/(2*Cs4) - (u1*u1+v1*v1)/(2*Cs2))
        # eq. 7 of Ref.2:
        F_phi[:,:,k] = ( 1 - 4 *(phi-PhaseField_av)**2)/interface_thickness * w[k] *(cx[k]*normal_x + cy[k]*normal_y)
        # eq. 9 of Ref.2:
        h_eq[:,:,k] = Gamma*phi - 0.5 * F_phi[:,:,k]
        # force calculation --> eq.(15) Ref.2
        force[:,:,k] = dt * w[k] * (cx[k]*Ft_x + cy[k]*Ft_y)/ (rho*Cs2)
        # g_eq -->  eq.(16) Ref. 2:
        g_eq[:,:,k] = P*w[k] + (Gamma-w[k]) - 0.5 * force[:,:,k]
    
    # collision step for phase field --> eq.(6) Ref.2
    h_pc = h - (h-h_eq)/(tau_phi+0.5) + F_phi
    # collision step for hydrodynamics --> eq.(14) Ref.2
    tau_rho = visc / (rho*Cs2) # relaxation time --> eq.(25) Ref.2
    Sv = 1/(tau_rho + 0.5)     # # eq.(29) Ref.2
    
    # MRT collision operator for Hydrodynamic --> eq.(27)Ref.2
    for i in range(nx):
        for j in range(ny):
            if struct[j,i]:
                S_hat[7,7] = Sv[j,i]  #  eq.(28) Ref.2
                S_hat[8,8] = Sv[j,i]  #  eq.(28) Ref.2
                minus_M_S_M = (-1) * matmul( matmul(M_inverse, S_hat) , M )                
                collision_operator[j,i,:] = matmul(minus_M_S_M , (g[j,i,:]-g_eq[j,i,:]))
    
    # =============================================================================
    #     # Single Relaxation-time: BGK
    #     for j in range(ny):
    #         for i in range(nx):
    #             collision_operator[j,i,:] = - (g[j,i,:]-g_eq[j,i,:])/(tau_rho[j,i]+0.5)
    # =============================================================================
     
    
    g_pc = g + collision_operator + force # eq.(14)Ref.2
    ###########################################################################
    
    ## Streaming ##############################################################
    # B.C treatment
    # inlet and outlet: Periodic 
    g_pc[:,0,:] = g_pc[:,nx-2,:] # Left boundary periodic
    g_pc[:,nx-1,:] = g_pc[:,1,:] # right boundary periodic
    h_pc[:,0,:] = h_pc[:,nx-2,:]
    h_pc[:,nx-1,:] = h_pc[:,1,:]
    # solid boundaries: replacing the distribution of solid boundary nodes
    for j in range(ny):
        for i in range(nx): 
            if coords_solid_boundaries[j,i,0]!=-100 : # if this is a solid boundari node
                for k in range(1,n_pop):
                    newy = j - cy[k]
                    newx = i + cx[k]
                    if newx==-1:
                        newx=nx-1
                    if newx==nx:
                        newx=0
                    if newy==-1:
                        newy=ny-1
                    if newy==ny:
                        newy=0    
                    g_pc[j,i,k] = g_pc[newy,newx,BounceBackD2Q9(k)] + force[newy,newx,BounceBackD2Q9(k)] #eq C8 Y. Q. Zu and S. He, 2013
                    h_pc[j,i,k] = h_pc[newy,newx,BounceBackD2Q9(k)] #eq C7b Y. Q. Zu and S. He, 2013

    # streamimng the fluid and solid boundary nodes 
    for j in range(ny):
        for i in range(nx):
            if struct[j,i] or (coords_solid_boundaries[j,i,0]!=-100) :  # if this is a fluid node or a solid boundary node 
                for k in range(n_pop):
                    newy = j - cy[k]
                    newx = i + cx[k]
                    if newx==-1:
                        newx=nx-1
                    if newx==nx:
                        newx=0
                    if newy==-1:
                        newy=ny-1
                    if newy==ny:
                        newy=0 
                    g[newy,newx,k] = g_pc[j,i,k] 
                    h[newy,newx,k] = h_pc[j,i,k] 
    
    ## Updating Macroscopic variables#########################################
    ## phase field:
    phi = sum(h, axis=2)
    #B.C:
    phi[:,0] = phi[:,nx-2] # Left boundary periodic
    phi[:,nx-1] = phi[:,1] # Left boundary periodic
    
    #phi[logical_and (struct==0,(coords_solid_boundaries[:,:,0]==-100)) ]=0
      
    # phase field for solid boundaries--> eq.(25) and eq.(26) 
    for j in range(ny):
        for i in range(nx):
            if coords_solid_boundaries[j,i,0] != -100 :
                newy = j - solid_normal_vectors_2d[j,i,0]
                newx = i + solid_normal_vectors_2d[j,i,1]
                normal_half_length = sqrt(solid_normal_vectors_2d[j,i,0]**2 + solid_normal_vectors_2d[j,i,1]**2) / 2
                if contact_angle == pi / 2:
                    phi[j,i] = phi[newy,newx]
                else: 
                    a1 = double(-1 * normal_half_length * sqrt(2*betha/kappa) * cos(contact_angle))
                    phi[j,i] = 1/a1 * (1 + a1 - sqrt((1+a1)**2-4*a1*phi[newy,newx] )) -  phi[newy,newx] 
                    if (1+a1)**2-4*a1*phi[newy,newx] < 0:
                        print('j= {0} , i= {1} , t= {2}\n' .format(j,i,t))
    
    # =============================================================================
    #     # calculating first and second derivative: (simple central difference + ignoring solid boundary nodes)
    #     phi_x_1st_derivative,phi_y_1st_derivative,phi_x_2nd_derivative,phi_y_2nd_derivative = phi_1st_2nd_derivative(phi,struct)
    # =============================================================================
    # =============================================================================
    #     #calculating first and second derivatives (simple central difference + considering solid boundary nodes)
    #     phi_x_1st_derivative,phi_y_1st_derivative,phi_x_2nd_derivative,phi_y_2nd_derivative = phi_1st_2nd_SimplewithSolid(phi,struct)
    #     
    # =============================================================================
    # calculating first derivative of x,y and laplacian ( isotrpic centered diffrence)
    phi_x_1st_derivative,phi_y_1st_derivative,laplacian_phi = phi_1st_2nd_isotropic(phi,struct)
    # normal vector --> eq.(2)
    magnitude_of_grad_phi = sqrt(phi_x_1st_derivative**2 + phi_y_1st_derivative**2) + epsilon
    normal_x = phi_x_1st_derivative / (magnitude_of_grad_phi )
    normal_y = phi_y_1st_derivative / (magnitude_of_grad_phi )
    
    
    ## Density calculation --> eq.(13) which was incorrect in paper --> 
    rho = Density_l + (Density_h - Density_l) * (phi - PhaseField_l)/(PhaseField_h - PhaseField_l)
    ## viscosity calculation --> eq.(19)
    visc = visc_L + (visc_H - visc_L) * (phi - PhaseField_l)/(PhaseField_h - PhaseField_l)
    visc [phi<PhaseField_l]=visc_L
    visc [phi>PhaseField_h]=visc_H
    ## Pressure
    P = getPressure(g)
    ## Forces
    Ft_x, Ft_y = update_force(P,rho,phi,visc,u1,v1,g,phi_x_1st_derivative,phi_y_1st_derivative,laplacian_phi)
    
    ## velocity
    g_right = g[:,:,[1,5,8]]
    g_left = g[:,:,[3,6,7]]
    u1 = sum(g_right, axis=2) - sum(g_left, axis=2) + dt*Ft_x/(2*rho)     # eq.(32b) ref.2

    g_up = g[:,:,[2,5,6]]
    g_down = g[:,:,[4,7,8]]
    v1 = sum(g_up, axis=2) - sum(g_down, axis=2)  + dt*Ft_y/(2*rho)

    u1[struct==0]=0
    v1[struct==0]=0
    
    #B.C:
    u1[:,0] = u1[:,nx-2] # Left boundary periodic
    u1[:,nx-1] = u1[:,1] # right boundary periodic
    
    v1[:,0] = v1[:,nx-2] # Left boundary periodic
    v1[:,nx-1] = v1[:,1] # right boundary periodic

    # Setting velocity of inlet and outlet boundary condition
    #u1[1:ny-1,0] = U_inlet 
    #v1[:,0] = 0
    
    if (mod(t,50)==0)&(t!=0):
        xx = time()
        conv= (mean(sum(u1,axis=0))-mean(sum(u1_old,axis=0)))/mean(sum(u1,axis=0)) ;
        print('Conv= {0} \t\t runtime= {1}' .format(conv,xx-x))
        # Checking convergence
        if abs(conv)<1e-4:
            break;
        else:
            u1_old=u1
            
    
results={}
results['phase_field']=phi
results['velocity_x']=u1
results['velocity_y']=v1
results['density']=rho
results['Pressure']=P
fname_pkl='RESULTS/results.pkl'


with open(fname_pkl,'wb') as fwrite:
    pickle.dump(results,fwrite)
    


    

        













