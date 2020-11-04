from matplotlib.pyplot import imshow, subplots, subplot, scatter
from numpy import array , ones, zeros, sum, unique, array_equal, append, double
from tkinter import filedialog
from pandas import read_excel


#Adress = filedialog.askopenfilename()
#excel_file = read_excel(Adress, sheet_name='300_400' , header=None)
#Style = array(excel_file.values)
#Style = Style.astype(int)



# reading normal vector database
#adress_normal_vector = filedialog.askopenfilename()
#adress_normal_vector = 'C:/Users/Roozshenas/Desktop/LBMPython/Free-Energy/AllNormal.xlsx'
#normal_excel = read_excel(adress_normal_vector, sheet_name='All' , header=None)
#all_normal_coord = array(normal_excel.values)
#all_normal_coord = all_normal_coord.astype(int)
#normal_style_excel = read_excel(adress_normal_vector, sheet_name='AllStyles', header=None)
#all_normal_style = array(normal_style_excel.values)
#all_normal_style = all_normal_style.astype(int)
    
cx = array([0,1,0,-1,0,1,-1,-1,1], dtype=int)
cy = array([0,0,1,0,-1,1,1,-1,-1], dtype=int)
w = array([4/9,1/9,1/9,1/9,1/9,1/36,1/36,1/36,1/36], dtype=double)
n_pop = len(cx)
    
class normal_vector:
    def __init__(self):
        # D2Q9 Parameters
        self.cx = array([0,1,0,-1,0,1,-1,-1,1], dtype=int)
        self.cy = array([0,0,1,0,-1,1,1,-1,-1], dtype=int)
        self.w = array([4/9,1/9,1/9,1/9,1/9,1/36,1/36,1/36,1/36], dtype=double)
        self.n_pop = len(self.cx)
        
    # this mehod finds solid nodes with at least one fluid node neighbour
    def linear_boundary_finder(self,Style):
        self.ny = len(Style)
        self.nx = len(Style[0])
        x_boundary = ones(1, dtype=int) # to store the i component of a boundary node in Style
        y_boundary = ones(1, dtype=int) # to store the j component of a boundary node in Style
        neighbour_arrangement = zeros((1,self.n_pop), dtype=int ) # to stor the arrangement of the neighbours of boundary nodes
        
        #coordinates_of_solid_boundary = ones((self.ny,self.nx,2), dtype= int) * (-100)  # to store i and j of solid boundary in a 2d array form --> the -100 values denotes that this node is not a solid boundary node
                                                                  # coordinates_of_solid_boundary[:,:,0] == "j"s of sloid boundaries
                                                                  # coordinates_of_solid_boundary[:,:,1] == "i"s of solid boundary
        #neighbour_arrangement_3d = zeros((self.ny,self.nx,self.n_pop), dtype=int)    # to store the arrangement of neighbours of a solid boundary node in 3d form                                                   
        
        counter = int(0)    
        for i in range(1,self.nx-1):
            for j in range(1,self.ny-1):
                if Style[j,i]==0 : # if a node is solid, go through this block
                    neigh = array([0, # getting all the neighbour of a solid node
                             Style[j-self.cy[1],i+self.cx[1]],
                             Style[j-self.cy[2],i+self.cx[2]],
                             Style[j-self.cy[3],i+self.cx[3]],
                             Style[j-self.cy[4],i+self.cx[4]],
                             Style[j-self.cy[5],i+self.cx[5]],
                             Style[j-self.cy[6],i+self.cx[6]],
                             Style[j-self.cy[7],i+self.cx[7]],
                             Style[j-self.cy[8],i+self.cx[8]]], dtype = int)
                    if sum(neigh, axis=0) > 0: # checking if at least one one fluid node exists around the solid boundary node. if exists go through this:
                        if counter == 0 :
                            x_boundary[counter] = i # getting i component of a boundary node having at least one fluid neighbour 
                            y_boundary[counter] = j # getting j component of a boundary node having at least one fluid neighbour 
                            neighbour_arrangement[counter,:] = neigh # getting arrangement of all neighbour nodes for a solid node that at least has one fluid neighbour
                            #coordinates_of_solid_boundary[j,i,:] = [j,i]
                            #neighbour_arrangement_3d[j,i,:] = neigh
                            counter += 1 
                        else:
                            x_boundary = append(x_boundary, [i], axis=0)
                            y_boundary = append(y_boundary, [j], axis=0)
                            neighbour_arrangement = append(neighbour_arrangement, [neigh], axis=0)
                            #coordinates_of_solid_boundary[j,i,:] = [j,i]
                            #neighbour_arrangement_3d[j,i,:] = neigh
                            counter += 1
        # for upper most and lower most row --> which are fully solid
        walls_y = array([0, self.ny-1], dtype=int)                  
        for j in walls_y :
            for i in range(1,self.nx-1):
                if Style[j,i]==0 :
                    if j == 0 :
                        neigh = [0,0,0,0,Style[j-self.cy[4],i+self.cx[4]],0,0,Style[j-self.cy[7],i+self.cx[7]],Style[j-self.cy[8],i+self.cx[8]]]
                    elif j == self.ny-1 :
                        neigh = [0,0,Style[j-self.cy[2],i+self.cx[2]],0,0,Style[j-self.cy[5],i+self.cx[5]],Style[j-self.cy[6],i+self.cx[6]],0,0]        
                    
                    if counter == 0 :
                        x_boundary[counter] = i 
                        y_boundary[counter] = j
                        neighbour_arrangement[counter,:] = neigh
                        #coordinates_of_solid_boundary[j,i,:] = [j,i]
                        #neighbour_arrangement_3d[j,i,:] = neigh
                        counter += 1 
                    else:
                        x_boundary = append(x_boundary, [i], axis=0)
                        y_boundary = append(y_boundary, [j], axis=0)
                        neighbour_arrangement = append(neighbour_arrangement, [neigh], axis=0)
                        #coordinates_of_solid_boundary[j,i,:] = [j,i]
                        #neighbour_arrangement_3d[j,i,:] = neigh
                        counter += 1              
        
        # asumption: i=0,1,nx-1,nx-2 are all set to be fluid nodes ==> there is no normal vector in i=0,nx-1
        # for Conrners:
        corners_x = array([0, self.nx-1, 0, self.nx-1], dtype=int)
        corners_y = array([0, 0, self.ny-1, self.ny-1], dtype=int)
        for i,j in zip(corners_x,corners_y): 
            if Style[j,i]==0 :
                if i==0 and j==0 :
                    neigh = [0,Style[j-self.cy[1],i+self.cx[1]],0,0,Style[j-self.cy[4],i+self.cx[4]],0,0,0,Style[j-self.cy[8],i+self.cx[8]]]
                if i==self.nx-1 and j==0:
                    neigh = [0,0,0,Style[j-self.cy[3],i+self.cx[3]],Style[j-self.cy[4],i+self.cx[4]],0,0,Style[j-self.cy[7],i+self.cx[7]],0]
                if i==0 and j==self.ny-1:
                    neigh = [0,Style[j-self.cy[1],i+self.cx[1]],Style[j-self.cy[2],i+self.cx[2]],0,0,Style[j-self.cy[5],i+self.cx[5]],0,0,0]
                if i==self.nx-1 and j==self.ny-1:
                    neigh = [0,0,Style[j-self.cy[2],i+self.cx[2]],Style[j-self.cy[3],i+self.cx[3]],0,0,Style[j-self.cy[6],i+self.cx[6]],0,0]
                
                if counter == 0 :
                    x_boundary[counter] = i 
                    y_boundary[counter] = j
                    neighbour_arrangement[counter,:] = neigh
                    #coordinates_of_solid_boundary[j,i,:] = [j,i]
                    #neighbour_arrangement_3d[j,i,:] = neigh
                    counter += 1 
                else:
                    x_boundary = append(x_boundary, [i], axis=0)
                    y_boundary = append(y_boundary, [j], axis=0)
                    neighbour_arrangement = append(neighbour_arrangement, [neigh], axis=0)
                    #coordinates_of_solid_boundary[j,i,:] = [j,i]
                    #neighbour_arrangement_3d[j,i,:] = neigh
                    counter += 1       
                    
        return(x_boundary,y_boundary,neighbour_arrangement)
        
    def matrix_boundary_finder(self,Style):
        self.ny = len(Style)
        self.nx = len(Style[0])
        #x_boundary = ones(1, dtype=int) # to store the i component of a boundary node in Style
        #y_boundary = ones(1, dtype=int) # to store the j component of a boundary node in Style
        #neighbour_arrangement = zeros((1,self.n_pop), dtype=int ) # to stor the arrangement of the neighbours of boundary nodes
        coordinates_of_solid_boundary = ones((self.ny,self.nx,2), dtype= int) * (-100)  # to store i and j of solid boundary in a 2d array form --> the -100 values denotes that this node is not a solid boundary node
                                                                  # coordinates_of_solid_boundary[:,:,0] == "j"s of sloid boundaries
                                                                  # coordinates_of_solid_boundary[:,:,1] == "i"s of solid boundary
        neighbour_arrangement_3d = zeros((self.ny,self.nx,self.n_pop), dtype=int)    # to store the arrangement of neighbours of a solid boundary node in 3d form                                                   
        counter = int(0)    
        for i in range(1,self.nx-1):
            for j in range(1,self.ny-1):
                if Style[j,i]==0 : # if a node is solid, go through this block
                    neigh = array([0, # getting all the neighbour of a solid node
                             Style[j-self.cy[1],i+self.cx[1]],
                             Style[j-self.cy[2],i+self.cx[2]],
                             Style[j-self.cy[3],i+self.cx[3]],
                             Style[j-self.cy[4],i+self.cx[4]],
                             Style[j-self.cy[5],i+self.cx[5]],
                             Style[j-self.cy[6],i+self.cx[6]],
                             Style[j-self.cy[7],i+self.cx[7]],
                             Style[j-self.cy[8],i+self.cx[8]]], dtype = int)
                    if sum(neigh, axis=0) > 0: # checking if at least one one fluid node exists around the solid boundary node. if exists go through this:
                        if counter == 0 :
                            #x_boundary[counter] = i # getting i component of a boundary node having at least one fluid neighbour 
                            #y_boundary[counter] = j # getting j component of a boundary node having at least one fluid neighbour 
                            #neighbour_arrangement[counter,:] = neigh # getting arrangement of all neighbour nodes for a solid node that at least has one fluid neighbour
                            coordinates_of_solid_boundary[j,i,:] = [j,i]
                            neighbour_arrangement_3d[j,i,:] = neigh
                            counter += 1 
                        else:
                            #x_boundary = append(x_boundary, [i], axis=0)
                            #y_boundary = append(y_boundary, [j], axis=0)
                            #neighbour_arrangement = append(neighbour_arrangement, [neigh], axis=0)
                            coordinates_of_solid_boundary[j,i,:] = [j,i]
                            neighbour_arrangement_3d[j,i,:] = neigh
                            counter += 1
        # for upper most and lower most row --> which are fully solid
        walls_y = array([0, self.ny-1], dtype=int)                  
        for j in walls_y :
            for i in range(1,self.nx-1):
                if Style[j,i]==0 :
                    if j == 0 :
                        neigh = [0,0,0,0,Style[j-self.cy[4],i+self.cx[4]],0,0,Style[j-self.cy[7],i+self.cx[7]],Style[j-self.cy[8],i+self.cx[8]]]
                    elif j == self.ny-1 :
                        neigh = [0,0,Style[j-self.cy[2],i+self.cx[2]],0,0,Style[j-self.cy[5],i+self.cx[5]],Style[j-self.cy[6],i+self.cx[6]],0,0]        
                    
                    if sum(neigh, axis=0) > 0: # checking if at least one one fluid node exists around the solid boundary node. if exists go through this:
                        if counter == 0 :
                            #x_boundary[counter] = i 
                            #y_boundary[counter] = j
                            #neighbour_arrangement[counter,:] = neigh
                            coordinates_of_solid_boundary[j,i,:] = [j,i]
                            neighbour_arrangement_3d[j,i,:] = neigh
                            counter += 1 
                        else:
                            #x_boundary = append(x_boundary, [i], axis=0)
                            #y_boundary = append(y_boundary, [j], axis=0)
                            #neighbour_arrangement = append(neighbour_arrangement, [neigh], axis=0)
                            coordinates_of_solid_boundary[j,i,:] = [j,i]
                            neighbour_arrangement_3d[j,i,:] = neigh
                            counter += 1              
            
        # asumption: i=0,1,nx-1,nx-2 are all set to be fluid nodes ==> there is no normal vector in i=0,nx-1
        # for Conrners:
        corners_x = array([0, self.nx-1, 0, self.nx-1], dtype=int)
        corners_y = array([0, 0, self.ny-1, self.ny-1], dtype=int)
        for i,j in zip(corners_x,corners_y): 
            if Style[j,i]==0 :
                if i==0 and j==0 :
                    neigh = [0,Style[j-self.cy[1],i+self.cx[1]],0,0,Style[j-self.cy[4],i+self.cx[4]],0,0,0,Style[j-self.cy[8],i+self.cx[8]]]
                if i==self.nx-1 and j==0:
                    neigh = [0,0,0,Style[j-self.cy[3],i+self.cx[3]],Style[j-self.cy[4],i+self.cx[4]],0,0,Style[j-self.cy[7],i+self.cx[7]],0]
                if i==0 and j==self.ny-1:
                    neigh = [0,Style[j-self.cy[1],i+self.cx[1]],Style[j-self.cy[2],i+self.cx[2]],0,0,Style[j-self.cy[5],i+self.cx[5]],0,0,0]
                if i==self.nx-1 and j==self.ny-1:
                    neigh = [0,0,Style[j-self.cy[2],i+self.cx[2]],Style[j-self.cy[3],i+self.cx[3]],0,0,Style[j-self.cy[6],i+self.cx[6]],0,0]
                
                if sum(neigh, axis=0) > 0: # checking if at least one one fluid node exists around the solid boundary node. if exists go through this:
                    if counter == 0 :
                        #x_boundary[counter] = i 
                        #y_boundary[counter] = j
                        #neighbour_arrangement[counter,:] = neigh
                        coordinates_of_solid_boundary[j,i,:] = [j,i]
                        neighbour_arrangement_3d[j,i,:] = neigh
                        counter += 1 
                    else:
                        #x_boundary = append(x_boundary, [i], axis=0)
                        #y_boundary = append(y_boundary, [j], axis=0)
                        #neighbour_arrangement = append(neighbour_arrangement, [neigh], axis=0)
                        coordinates_of_solid_boundary[j,i,:] = [j,i]
                        neighbour_arrangement_3d[j,i,:] = neigh
                        counter += 1       
                    
        return(coordinates_of_solid_boundary, neighbour_arrangement_3d)

    
    # this method has the data bank of all normal vectors and all possible arrangements for neighbour nodes of a boundary node
    def normal_bank_1d(self,neighbour_arrangement): 
        all_normal_coord =array( [  # all possible normal vectors
                                	[1,1],	[1,-1],	[-1,-1],[-1,-1],[-1,1],	[-1,1], [-1,1],[-1,1],[1,1],[1,1],[1,1],[1,1],[1,1],[1,1],[1,1],[1,1],[0,-1],[0,-1],[0,-1],[0,-1],[0,-1],[0,-1],[0,-1],[0,-1],[0,-1],[0,-1],[0,-1],[0,-1],[0,-1],[0,-1],[0,-1],[0,-1],
                                	[-1,0],	[-1,0],	[-1,0],	[-1,0],	[-1,0],	[-1,0],[-1,0],[-1,0],[-1,0],[-1,0],[-1,0],[-1,0],[-1,0],[-1,0],[-1,0],[-1,0],[-1,0],[0,-1],[-1,-1],[-1,-1],[-1,0],[-1,0],[-1,-1],[-1,-1],[-1,0],[0,-1],[-1,-1],[-1,-1],[-1,0],[-1,0],[-1,-1],[-1,-1],
                                	[0,1],	[0,1],	[0,1],	[0,1],	[0,1],	[0,1],[0,1],[0,1],[0,1],[0,1],[0,1],[0,1],[0,1],[0,1],[0,1],[0,1],[0,1],[0,-1],[0,-1],[0,-1],[0,1],[0,1],[0,1],[0,-1],[0,1],[0,1],[0,1],[0,-1],[0,1],[0,1],[0,1],[0,1],
                                	[0,1],	[0,1],	[-1,0],	[-1,0],	[-1,1],	[-1,1],[-1,1],[-1,1],[0,1],[0,1],[0,1],[0,1],[-1,1],[-1,1],[-1,1],[-1,1],[0,1],[0,-1],[-1,-1],[-1,-1],[-1,1],[-1,1],[-1,0],[-1,0],[0,1],[0,1],[-1,-1],[-1,-1],[-1,1],[-1,1],[-1,0],[-1,0],
                                	[1,0],	[1,0],	[1,0],	[1,0],	[1,0],	[1,0],[1,0],[1,0],[1,0],[1,0],[1,0],[1,0],[1,0],[1,0],[1,0],[1,0],[1,0],[1,-1],[0,-1],[1,-1],[1,0],[1,-1],[0,-1],[1,-1],[1,0],[1,-1],[1,0],[1,-1],[1,0],[1,-1],[1,0],[1,-1],
                                	[1,0],	[1,0],	[-1,0],	[1,0],	[-1,0],	[1,0],[-1,0],[-1,0],[1,0],[1,0],[1,0],[1,0],[1,0],[1,0],[-1,0],[1,0],[1,0],[1,-1],[-1,-1],[0,-1],[-1,0],[1,-1],[-1,-1],[0,-1],[1,0],[1,-1],[-1,-1],[0,-1],[1,0],[1,-1],[-1,-1],[0,-1],
                                	[1,0],	[1,0],	[1,0],	[1,0],	[0,1],	[1,0],[0,1],[1,0],[1,1],[1,1],[1,1],[1,1],[1,1],[1,1],[1,1],[1,1],[1,0],[1,-1],[0,-1],[1,-1],[0,1],[1,-1],[0,1],[1,-1],[1,1],[1,0],[1,1],[1,0],[1,1],[1,0],[1,1],[1,0],
                                	[1,0],	[1,0],	[-1,0],	[1,0],	[-1,1],	[-1,1],[-1,1],[-1,1],[1,1],[1,1],[1,1],[1,1],[0,1],[0,1],[0,1],[0,1],[1,0],[1,-1],[-1,-1],[0,-1],[-1,1],[-1,1],[-1,0],[-1,-1],[1,1],[1,0],[1,1],[1,-1],[0,1],[1,1],[-1,1],[1,0]] , dtype=int )
                        
        all_normal_style = array([ # all possible arrangements of neighbour nodes
                                    [0,0,0,0,0,0,0,0],[0,0,0,0,0,0,0,1],[0,0,0,0,0,0,1,0],[0,0,0,0,0,0,1,1],[0,0,0,0,0,1,0,0],[0,0,0,0,0,1,0,1],[0,0,0,0,0,1,1,0],[0,0,0,0,0,1,1,1],
                                    [0,0,0,0,1,0,0,0],[0,0,0,0,1,0,0,1],[0,0,0,0,1,0,1,0],[0,0,0,0,1,0,1,1],[0,0,0,0,1,1,0,0],[0,0,0,0,1,1,0,1],[0,0,0,0,1,1,1,0],[0,0,0,0,1,1,1,1],
                                    [0,0,0,1,0,0,0,0],[0,0,0,1,0,0,0,1],[0,0,0,1,0,0,1,0],[0,0,0,1,0,0,1,1],[0,0,0,1,0,1,0,0],[0,0,0,1,0,1,0,1],[0,0,0,1,0,1,1,0],[0,0,0,1,0,1,1,1],
                                    [0,0,0,1,1,0,0,0],[0,0,0,1,1,0,0,1],[0,0,0,1,1,0,1,0],[0,0,0,1,1,0,1,1],[0,0,0,1,1,1,0,0],[0,0,0,1,1,1,0,1],[0,0,0,1,1,1,1,0],[0,0,0,1,1,1,1,1],
                                    [0,0,1,0,0,0,0,0],[0,0,1,0,0,0,0,1],[0,0,1,0,0,0,1,0],[0,0,1,0,0,0,1,1],[0,0,1,0,0,1,0,0],[0,0,1,0,0,1,0,1],[0,0,1,0,0,1,1,0],[0,0,1,0,0,1,1,1],
                                    [0,0,1,0,1,0,0,0],[0,0,1,0,1,0,0,1],[0,0,1,0,1,0,1,0],[0,0,1,0,1,0,1,1],[0,0,1,0,1,1,0,0],[0,0,1,0,1,1,0,1],[0,0,1,0,1,1,1,0],[0,0,1,0,1,1,1,1],
                                    [0,0,1,1,0,0,0,0],[0,0,1,1,0,0,0,1],[0,0,1,1,0,0,1,0],[0,0,1,1,0,0,1,1],[0,0,1,1,0,1,0,0],[0,0,1,1,0,1,0,1],[0,0,1,1,0,1,1,0],[0,0,1,1,0,1,1,1],
                                    [0,0,1,1,1,0,0,0],[0,0,1,1,1,0,0,1],[0,0,1,1,1,0,1,0],[0,0,1,1,1,0,1,1],[0,0,1,1,1,1,0,0],[0,0,1,1,1,1,0,1],[0,0,1,1,1,1,1,0],[0,0,1,1,1,1,1,1],
                                    [0,1,0,0,0,0,0,0],[0,1,0,0,0,0,0,1],[0,1,0,0,0,0,1,0],[0,1,0,0,0,0,1,1],[0,1,0,0,0,1,0,0],[0,1,0,0,0,1,0,1],[0,1,0,0,0,1,1,0],[0,1,0,0,0,1,1,1],
                                    [0,1,0,0,1,0,0,0],[0,1,0,0,1,0,0,1],[0,1,0,0,1,0,1,0],[0,1,0,0,1,0,1,1],[0,1,0,0,1,1,0,0],[0,1,0,0,1,1,0,1],[0,1,0,0,1,1,1,0],[0,1,0,0,1,1,1,1],
                                    [0,1,0,1,0,0,0,0],[0,1,0,1,0,0,0,1],[0,1,0,1,0,0,1,0],[0,1,0,1,0,0,1,1],[0,1,0,1,0,1,0,0],[0,1,0,1,0,1,0,1],[0,1,0,1,0,1,1,0],[0,1,0,1,0,1,1,1],
                                    [0,1,0,1,1,0,0,0],[0,1,0,1,1,0,0,1],[0,1,0,1,1,0,1,0],[0,1,0,1,1,0,1,1],[0,1,0,1,1,1,0,0],[0,1,0,1,1,1,0,1],[0,1,0,1,1,1,1,0],[0,1,0,1,1,1,1,1],
                                    [0,1,1,0,0,0,0,0],[0,1,1,0,0,0,0,1],[0,1,1,0,0,0,1,0],[0,1,1,0,0,0,1,1],[0,1,1,0,0,1,0,0],[0,1,1,0,0,1,0,1],[0,1,1,0,0,1,1,0],[0,1,1,0,0,1,1,1],
                                    [0,1,1,0,1,0,0,0],[0,1,1,0,1,0,0,1],[0,1,1,0,1,0,1,0],[0,1,1,0,1,0,1,1],[0,1,1,0,1,1,0,0],[0,1,1,0,1,1,0,1],[0,1,1,0,1,1,1,0],[0,1,1,0,1,1,1,1],
                                    [0,1,1,1,0,0,0,0],[0,1,1,1,0,0,0,1],[0,1,1,1,0,0,1,0],[0,1,1,1,0,0,1,1],[0,1,1,1,0,1,0,0],[0,1,1,1,0,1,0,1],[0,1,1,1,0,1,1,0],[0,1,1,1,0,1,1,1],
                                    [0,1,1,1,1,0,0,0],[0,1,1,1,1,0,0,1],[0,1,1,1,1,0,1,0],[0,1,1,1,1,0,1,1],[0,1,1,1,1,1,0,0],[0,1,1,1,1,1,0,1],[0,1,1,1,1,1,1,0],[0,1,1,1,1,1,1,1],
                                    [1,0,0,0,0,0,0,0],[1,0,0,0,0,0,0,1],[1,0,0,0,0,0,1,0],[1,0,0,0,0,0,1,1],[1,0,0,0,0,1,0,0],[1,0,0,0,0,1,0,1],[1,0,0,0,0,1,1,0],[1,0,0,0,0,1,1,1],
                                    [1,0,0,0,1,0,0,0],[1,0,0,0,1,0,0,1],[1,0,0,0,1,0,1,0],[1,0,0,0,1,0,1,1],[1,0,0,0,1,1,0,0],[1,0,0,0,1,1,0,1],[1,0,0,0,1,1,1,0],[1,0,0,0,1,1,1,1],
                                    [1,0,0,1,0,0,0,0],[1,0,0,1,0,0,0,1],[1,0,0,1,0,0,1,0],[1,0,0,1,0,0,1,1],[1,0,0,1,0,1,0,0],[1,0,0,1,0,1,0,1],[1,0,0,1,0,1,1,0],[1,0,0,1,0,1,1,1],
                                    [1,0,0,1,1,0,0,0],[1,0,0,1,1,0,0,1],[1,0,0,1,1,0,1,0],[1,0,0,1,1,0,1,1],[1,0,0,1,1,1,0,0],[1,0,0,1,1,1,0,1],[1,0,0,1,1,1,1,0],[1,0,0,1,1,1,1,1],
                                    [1,0,1,0,0,0,0,0],[1,0,1,0,0,0,0,1],[1,0,1,0,0,0,1,0],[1,0,1,0,0,0,1,1],[1,0,1,0,0,1,0,0],[1,0,1,0,0,1,0,1],[1,0,1,0,0,1,1,0],[1,0,1,0,0,1,1,1],
                                    [1,0,1,0,1,0,0,0],[1,0,1,0,1,0,0,1],[1,0,1,0,1,0,1,0],[1,0,1,0,1,0,1,1],[1,0,1,0,1,1,0,0],[1,0,1,0,1,1,0,1],[1,0,1,0,1,1,1,0],[1,0,1,0,1,1,1,1],
                                    [1,0,1,1,0,0,0,0],[1,0,1,1,0,0,0,1],[1,0,1,1,0,0,1,0],[1,0,1,1,0,0,1,1],[1,0,1,1,0,1,0,0],[1,0,1,1,0,1,0,1],[1,0,1,1,0,1,1,0],[1,0,1,1,0,1,1,1],
                                    [1,0,1,1,1,0,0,0],[1,0,1,1,1,0,0,1],[1,0,1,1,1,0,1,0],[1,0,1,1,1,0,1,1],[1,0,1,1,1,1,0,0],[1,0,1,1,1,1,0,1],[1,0,1,1,1,1,1,0],[1,0,1,1,1,1,1,1],
                                    [1,1,0,0,0,0,0,0],[1,1,0,0,0,0,0,1],[1,1,0,0,0,0,1,0],[1,1,0,0,0,0,1,1],[1,1,0,0,0,1,0,0],[1,1,0,0,0,1,0,1],[1,1,0,0,0,1,1,0],[1,1,0,0,0,1,1,1],
                                    [1,1,0,0,1,0,0,0],[1,1,0,0,1,0,0,1],[1,1,0,0,1,0,1,0],[1,1,0,0,1,0,1,1],[1,1,0,0,1,1,0,0],[1,1,0,0,1,1,0,1],[1,1,0,0,1,1,1,0],[1,1,0,0,1,1,1,1],
                                    [1,1,0,1,0,0,0,0],[1,1,0,1,0,0,0,1],[1,1,0,1,0,0,1,0],[1,1,0,1,0,0,1,1],[1,1,0,1,0,1,0,0],[1,1,0,1,0,1,0,1],[1,1,0,1,0,1,1,0],[1,1,0,1,0,1,1,1],
                                    [1,1,0,1,1,0,0,0],[1,1,0,1,1,0,0,1],[1,1,0,1,1,0,1,0],[1,1,0,1,1,0,1,1],[1,1,0,1,1,1,0,0],[1,1,0,1,1,1,0,1],[1,1,0,1,1,1,1,0],[1,1,0,1,1,1,1,1],
                                    [1,1,1,0,0,0,0,0],[1,1,1,0,0,0,0,1],[1,1,1,0,0,0,1,0],[1,1,1,0,0,0,1,1],[1,1,1,0,0,1,0,0],[1,1,1,0,0,1,0,1],[1,1,1,0,0,1,1,0],[1,1,1,0,0,1,1,1],
                                    [1,1,1,0,1,0,0,0],[1,1,1,0,1,0,0,1],[1,1,1,0,1,0,1,0],[1,1,1,0,1,0,1,1],[1,1,1,0,1,1,0,0],[1,1,1,0,1,1,0,1],[1,1,1,0,1,1,1,0],[1,1,1,0,1,1,1,1],
                                    [1,1,1,1,0,0,0,0],[1,1,1,1,0,0,0,1],[1,1,1,1,0,0,1,0],[1,1,1,1,0,0,1,1],[1,1,1,1,0,1,0,0],[1,1,1,1,0,1,0,1],[1,1,1,1,0,1,1,0],[1,1,1,1,0,1,1,1],
                                    [1,1,1,1,1,0,0,0],[1,1,1,1,1,0,0,1],[1,1,1,1,1,0,1,0],[1,1,1,1,1,0,1,1],[1,1,1,1,1,1,0,0],[1,1,1,1,1,1,0,1],[1,1,1,1,1,1,1,0],[1,1,1,1,1,1,1,1] ] , dtype=int)            
                
                
        normal_list = zeros((len(neighbour_arrangement), 2), dtype=int)
        counter = 0
        for i in range(len(neighbour_arrangement)):
            index = all_normal_style.tolist().index(neighbour_arrangement[i,1:self.n_pop].tolist())
            normal_list[counter,:] = all_normal_coord[index,:]
            counter+=1
        return(normal_list)
        
    def normal_bank_2d(self,neighbour_arrangement_3d): 
        all_normal_coord =array( [  # all possible normal vectors # in each binary: first is "i" second number is "j"
                                	[1,1],	[1,-1],	[-1,-1],[-1,-1],[-1,1],	[-1,1], [-1,1],[-1,1],[1,1],[1,1],[1,1],[1,1],[1,1],[1,1],[1,1],[1,1],[0,-1],[0,-1],[0,-1],[0,-1],[0,-1],[0,-1],[0,-1],[0,-1],[0,-1],[0,-1],[0,-1],[0,-1],[0,-1],[0,-1],[0,-1],[0,-1],
                                	[-1,0],	[-1,0],	[-1,0],	[-1,0],	[-1,0],	[-1,0],[-1,0],[-1,0],[-1,0],[-1,0],[-1,0],[-1,0],[-1,0],[-1,0],[-1,0],[-1,0],[-1,0],[0,-1],[-1,-1],[-1,-1],[-1,0],[-1,0],[-1,-1],[-1,-1],[-1,0],[0,-1],[-1,-1],[-1,-1],[-1,0],[-1,0],[-1,-1],[-1,-1],
                                	[0,1],	[0,1],	[0,1],	[0,1],	[0,1],	[0,1],[0,1],[0,1],[0,1],[0,1],[0,1],[0,1],[0,1],[0,1],[0,1],[0,1],[0,1],[0,-1],[0,-1],[0,-1],[0,1],[0,1],[0,1],[0,-1],[0,1],[0,1],[0,1],[0,-1],[0,1],[0,1],[0,1],[0,1],
                                	[0,1],	[0,1],	[-1,0],	[-1,0],	[-1,1],	[-1,1],[-1,1],[-1,1],[0,1],[0,1],[0,1],[0,1],[-1,1],[-1,1],[-1,1],[-1,1],[0,1],[0,-1],[-1,-1],[-1,-1],[-1,1],[-1,1],[-1,0],[-1,0],[0,1],[0,1],[-1,-1],[-1,-1],[-1,1],[-1,1],[-1,0],[-1,0],
                                	[1,0],	[1,0],	[1,0],	[1,0],	[1,0],	[1,0],[1,0],[1,0],[1,0],[1,0],[1,0],[1,0],[1,0],[1,0],[1,0],[1,0],[1,0],[1,-1],[0,-1],[1,-1],[1,0],[1,-1],[0,-1],[1,-1],[1,0],[1,-1],[1,0],[1,-1],[1,0],[1,-1],[1,0],[1,-1],
                                	[1,0],	[1,0],	[-1,0],	[1,0],	[-1,0],	[1,0],[-1,0],[-1,0],[1,0],[1,0],[1,0],[1,0],[1,0],[1,0],[-1,0],[1,0],[1,0],[1,-1],[-1,-1],[0,-1],[-1,0],[1,-1],[-1,-1],[0,-1],[1,0],[1,-1],[-1,-1],[0,-1],[1,0],[1,-1],[-1,-1],[0,-1],
                                	[1,0],	[1,0],	[1,0],	[1,0],	[0,1],	[1,0],[0,1],[1,0],[1,1],[1,1],[1,1],[1,1],[1,1],[1,1],[1,1],[1,1],[1,0],[1,-1],[0,-1],[1,-1],[0,1],[1,-1],[0,1],[1,-1],[1,1],[1,0],[1,1],[1,0],[1,1],[1,0],[1,1],[1,0],
                                	[1,0],	[1,0],	[-1,0],	[1,0],	[-1,1],	[-1,1],[-1,1],[-1,1],[1,1],[1,1],[1,1],[1,1],[0,1],[0,1],[0,1],[0,1],[1,0],[1,-1],[-1,-1],[0,-1],[-1,1],[-1,1],[-1,0],[-1,-1],[1,1],[1,0],[1,1],[1,-1],[0,1],[1,1],[-1,1],[1,0]] , dtype=int )
            
            # each array contains 8 numbers which denote this order of D2Q9 directions: [1,2,3,4,5,6,7,8] , zero direction is igored here
        all_normal_style = array([ # all possible arrangements of neighbour nodes
                                    [0,0,0,0,0,0,0,0],[0,0,0,0,0,0,0,1],[0,0,0,0,0,0,1,0],[0,0,0,0,0,0,1,1],[0,0,0,0,0,1,0,0],[0,0,0,0,0,1,0,1],[0,0,0,0,0,1,1,0],[0,0,0,0,0,1,1,1],
                                    [0,0,0,0,1,0,0,0],[0,0,0,0,1,0,0,1],[0,0,0,0,1,0,1,0],[0,0,0,0,1,0,1,1],[0,0,0,0,1,1,0,0],[0,0,0,0,1,1,0,1],[0,0,0,0,1,1,1,0],[0,0,0,0,1,1,1,1],
                                    [0,0,0,1,0,0,0,0],[0,0,0,1,0,0,0,1],[0,0,0,1,0,0,1,0],[0,0,0,1,0,0,1,1],[0,0,0,1,0,1,0,0],[0,0,0,1,0,1,0,1],[0,0,0,1,0,1,1,0],[0,0,0,1,0,1,1,1],
                                    [0,0,0,1,1,0,0,0],[0,0,0,1,1,0,0,1],[0,0,0,1,1,0,1,0],[0,0,0,1,1,0,1,1],[0,0,0,1,1,1,0,0],[0,0,0,1,1,1,0,1],[0,0,0,1,1,1,1,0],[0,0,0,1,1,1,1,1],
                                    [0,0,1,0,0,0,0,0],[0,0,1,0,0,0,0,1],[0,0,1,0,0,0,1,0],[0,0,1,0,0,0,1,1],[0,0,1,0,0,1,0,0],[0,0,1,0,0,1,0,1],[0,0,1,0,0,1,1,0],[0,0,1,0,0,1,1,1],
                                    [0,0,1,0,1,0,0,0],[0,0,1,0,1,0,0,1],[0,0,1,0,1,0,1,0],[0,0,1,0,1,0,1,1],[0,0,1,0,1,1,0,0],[0,0,1,0,1,1,0,1],[0,0,1,0,1,1,1,0],[0,0,1,0,1,1,1,1],
                                    [0,0,1,1,0,0,0,0],[0,0,1,1,0,0,0,1],[0,0,1,1,0,0,1,0],[0,0,1,1,0,0,1,1],[0,0,1,1,0,1,0,0],[0,0,1,1,0,1,0,1],[0,0,1,1,0,1,1,0],[0,0,1,1,0,1,1,1],
                                    [0,0,1,1,1,0,0,0],[0,0,1,1,1,0,0,1],[0,0,1,1,1,0,1,0],[0,0,1,1,1,0,1,1],[0,0,1,1,1,1,0,0],[0,0,1,1,1,1,0,1],[0,0,1,1,1,1,1,0],[0,0,1,1,1,1,1,1],
                                    [0,1,0,0,0,0,0,0],[0,1,0,0,0,0,0,1],[0,1,0,0,0,0,1,0],[0,1,0,0,0,0,1,1],[0,1,0,0,0,1,0,0],[0,1,0,0,0,1,0,1],[0,1,0,0,0,1,1,0],[0,1,0,0,0,1,1,1],
                                    [0,1,0,0,1,0,0,0],[0,1,0,0,1,0,0,1],[0,1,0,0,1,0,1,0],[0,1,0,0,1,0,1,1],[0,1,0,0,1,1,0,0],[0,1,0,0,1,1,0,1],[0,1,0,0,1,1,1,0],[0,1,0,0,1,1,1,1],
                                    [0,1,0,1,0,0,0,0],[0,1,0,1,0,0,0,1],[0,1,0,1,0,0,1,0],[0,1,0,1,0,0,1,1],[0,1,0,1,0,1,0,0],[0,1,0,1,0,1,0,1],[0,1,0,1,0,1,1,0],[0,1,0,1,0,1,1,1],
                                    [0,1,0,1,1,0,0,0],[0,1,0,1,1,0,0,1],[0,1,0,1,1,0,1,0],[0,1,0,1,1,0,1,1],[0,1,0,1,1,1,0,0],[0,1,0,1,1,1,0,1],[0,1,0,1,1,1,1,0],[0,1,0,1,1,1,1,1],
                                    [0,1,1,0,0,0,0,0],[0,1,1,0,0,0,0,1],[0,1,1,0,0,0,1,0],[0,1,1,0,0,0,1,1],[0,1,1,0,0,1,0,0],[0,1,1,0,0,1,0,1],[0,1,1,0,0,1,1,0],[0,1,1,0,0,1,1,1],
                                    [0,1,1,0,1,0,0,0],[0,1,1,0,1,0,0,1],[0,1,1,0,1,0,1,0],[0,1,1,0,1,0,1,1],[0,1,1,0,1,1,0,0],[0,1,1,0,1,1,0,1],[0,1,1,0,1,1,1,0],[0,1,1,0,1,1,1,1],
                                    [0,1,1,1,0,0,0,0],[0,1,1,1,0,0,0,1],[0,1,1,1,0,0,1,0],[0,1,1,1,0,0,1,1],[0,1,1,1,0,1,0,0],[0,1,1,1,0,1,0,1],[0,1,1,1,0,1,1,0],[0,1,1,1,0,1,1,1],
                                    [0,1,1,1,1,0,0,0],[0,1,1,1,1,0,0,1],[0,1,1,1,1,0,1,0],[0,1,1,1,1,0,1,1],[0,1,1,1,1,1,0,0],[0,1,1,1,1,1,0,1],[0,1,1,1,1,1,1,0],[0,1,1,1,1,1,1,1],
                                    [1,0,0,0,0,0,0,0],[1,0,0,0,0,0,0,1],[1,0,0,0,0,0,1,0],[1,0,0,0,0,0,1,1],[1,0,0,0,0,1,0,0],[1,0,0,0,0,1,0,1],[1,0,0,0,0,1,1,0],[1,0,0,0,0,1,1,1],
                                    [1,0,0,0,1,0,0,0],[1,0,0,0,1,0,0,1],[1,0,0,0,1,0,1,0],[1,0,0,0,1,0,1,1],[1,0,0,0,1,1,0,0],[1,0,0,0,1,1,0,1],[1,0,0,0,1,1,1,0],[1,0,0,0,1,1,1,1],
                                    [1,0,0,1,0,0,0,0],[1,0,0,1,0,0,0,1],[1,0,0,1,0,0,1,0],[1,0,0,1,0,0,1,1],[1,0,0,1,0,1,0,0],[1,0,0,1,0,1,0,1],[1,0,0,1,0,1,1,0],[1,0,0,1,0,1,1,1],
                                    [1,0,0,1,1,0,0,0],[1,0,0,1,1,0,0,1],[1,0,0,1,1,0,1,0],[1,0,0,1,1,0,1,1],[1,0,0,1,1,1,0,0],[1,0,0,1,1,1,0,1],[1,0,0,1,1,1,1,0],[1,0,0,1,1,1,1,1],
                                    [1,0,1,0,0,0,0,0],[1,0,1,0,0,0,0,1],[1,0,1,0,0,0,1,0],[1,0,1,0,0,0,1,1],[1,0,1,0,0,1,0,0],[1,0,1,0,0,1,0,1],[1,0,1,0,0,1,1,0],[1,0,1,0,0,1,1,1],
                                    [1,0,1,0,1,0,0,0],[1,0,1,0,1,0,0,1],[1,0,1,0,1,0,1,0],[1,0,1,0,1,0,1,1],[1,0,1,0,1,1,0,0],[1,0,1,0,1,1,0,1],[1,0,1,0,1,1,1,0],[1,0,1,0,1,1,1,1],
                                    [1,0,1,1,0,0,0,0],[1,0,1,1,0,0,0,1],[1,0,1,1,0,0,1,0],[1,0,1,1,0,0,1,1],[1,0,1,1,0,1,0,0],[1,0,1,1,0,1,0,1],[1,0,1,1,0,1,1,0],[1,0,1,1,0,1,1,1],
                                    [1,0,1,1,1,0,0,0],[1,0,1,1,1,0,0,1],[1,0,1,1,1,0,1,0],[1,0,1,1,1,0,1,1],[1,0,1,1,1,1,0,0],[1,0,1,1,1,1,0,1],[1,0,1,1,1,1,1,0],[1,0,1,1,1,1,1,1],
                                    [1,1,0,0,0,0,0,0],[1,1,0,0,0,0,0,1],[1,1,0,0,0,0,1,0],[1,1,0,0,0,0,1,1],[1,1,0,0,0,1,0,0],[1,1,0,0,0,1,0,1],[1,1,0,0,0,1,1,0],[1,1,0,0,0,1,1,1],
                                    [1,1,0,0,1,0,0,0],[1,1,0,0,1,0,0,1],[1,1,0,0,1,0,1,0],[1,1,0,0,1,0,1,1],[1,1,0,0,1,1,0,0],[1,1,0,0,1,1,0,1],[1,1,0,0,1,1,1,0],[1,1,0,0,1,1,1,1],
                                    [1,1,0,1,0,0,0,0],[1,1,0,1,0,0,0,1],[1,1,0,1,0,0,1,0],[1,1,0,1,0,0,1,1],[1,1,0,1,0,1,0,0],[1,1,0,1,0,1,0,1],[1,1,0,1,0,1,1,0],[1,1,0,1,0,1,1,1],
                                    [1,1,0,1,1,0,0,0],[1,1,0,1,1,0,0,1],[1,1,0,1,1,0,1,0],[1,1,0,1,1,0,1,1],[1,1,0,1,1,1,0,0],[1,1,0,1,1,1,0,1],[1,1,0,1,1,1,1,0],[1,1,0,1,1,1,1,1],
                                    [1,1,1,0,0,0,0,0],[1,1,1,0,0,0,0,1],[1,1,1,0,0,0,1,0],[1,1,1,0,0,0,1,1],[1,1,1,0,0,1,0,0],[1,1,1,0,0,1,0,1],[1,1,1,0,0,1,1,0],[1,1,1,0,0,1,1,1],
                                    [1,1,1,0,1,0,0,0],[1,1,1,0,1,0,0,1],[1,1,1,0,1,0,1,0],[1,1,1,0,1,0,1,1],[1,1,1,0,1,1,0,0],[1,1,1,0,1,1,0,1],[1,1,1,0,1,1,1,0],[1,1,1,0,1,1,1,1],
                                    [1,1,1,1,0,0,0,0],[1,1,1,1,0,0,0,1],[1,1,1,1,0,0,1,0],[1,1,1,1,0,0,1,1],[1,1,1,1,0,1,0,0],[1,1,1,1,0,1,0,1],[1,1,1,1,0,1,1,0],[1,1,1,1,0,1,1,1],
                                    [1,1,1,1,1,0,0,0],[1,1,1,1,1,0,0,1],[1,1,1,1,1,0,1,0],[1,1,1,1,1,0,1,1],[1,1,1,1,1,1,0,0],[1,1,1,1,1,1,0,1],[1,1,1,1,1,1,1,0],[1,1,1,1,1,1,1,1] ] , dtype=int)            
                
        self.ny = len(neighbour_arrangement_3d)
        self.nx = len(neighbour_arrangement_3d[0])        
        normal_list = ones((self.ny,self.nx, 2), dtype=int) * (-100)  # k=0 plane is for "j"s  and k=1 plane is for "i"s of normal vectors
        for j in range(self.ny):
            for i in range(self.nx):
                neigh = neighbour_arrangement_3d[j,i,1:self.n_pop].tolist()
                index = all_normal_style.tolist().index(neigh)
                normal_list[j,i,0] = all_normal_coord[index,1] # j of normal vector
                normal_list[j,i,1] = all_normal_coord[index,0] # i of normal vector
        return(normal_list)
        
def BounceBackD2Q9(b):
    if b==1:
        a=3;
    elif b==2:
        a=4;
    elif b==3:
        a=1;
    elif b==4:
        a=2;
    elif b==5:
        a=7;
    elif b==6:
        a=8;
    elif b==7:
        a=5;
    elif b==8:
        a=6;
    elif b==0:
        a=0;     
    return(a)

def phi_1st_2nd_derivative(phi,Style):
    ny = len(Style)
    nx = len(Style[0])
    
    x_1st_derivative = zeros((ny,nx), dtype=double) * 100
    y_1st_derivative = zeros((ny,nx), dtype=double) * 100
    x_2nd_derivative = zeros((ny,nx), dtype=double) * 100 
    y_2nd_derivative = zeros((ny,nx), dtype=double) * 100
    
    
    for j in range(1,ny-1):
        for i in range(0,nx): # excluding inlet and outlet nodes 
            if Style[j,i]: # only if the node is a fluid nide, we calculate derivative
                
                if i!=0 and i!=nx-1 :    # excluding inlet and outlet from calculating x derivatives
                    # first derivative in x dir
                    if Style[j,i+1] and Style[j,i-1]: 
                        # central difference ==> 2nd order error
                        x_1st_derivative[j,i] = (phi[j,i+1] - phi[j,i-1])/2 
                    elif Style[j,i+1] and (Style[j,i-1]==0) :
                        # forward difference ==> 1st order error
                        x_1st_derivative[j,i] = (phi[j,i+1] - phi[j,i])/1
                    elif (Style[j,i+1]==0) and Style[j,i-1]: 
                        # backward difference ==> 1st order error
                        x_1st_derivative[j,i] = (phi[j,i] - phi[j,i-1])/1
                    else: # Style[j,i+1]==0 and Style[j,i-1]==0
                        x_1st_derivative[j,i] = 0
                        
                    # second derivative in x dir
                    if Style[j,i+1] and Style[j,i-1]:
                        # central ==> 2nd order error
                        x_2nd_derivative[j,i] = (phi[j,i+1] - 2*phi[j,i] + phi[j,i-1])/(1**2)
                    elif (i!=nx-2) and (Style[j,i+1] and Style[j,i+2]) and (Style[j,i-1]==0):
                        # forward ==> 1st order 
                        x_2nd_derivative[j,i] = (phi[j,i+2] - 2*phi[j,i+1] + phi[j,i])/(1**2)
                    elif (i!=1) and (Style[j,i+1]==0) and (Style[j,i-1] and Style[j,i-2]): 
                        # backward ==> 1st order
                        x_2nd_derivative[j,i] = (phi[j,i] - 2*phi[j,i-1] + phi[j,i-2])/1**2
                    else:   # for the nodes between two solid nodes, and for those node that dont have enough data to calculate forward and backward second derivatives
                        x_2nd_derivative[j,i] = 0 
                        
                elif i==0:
                    # first derivative in x dir
                    if Style[j,i+1]  :
                        # forward difference ==> 1st order error
                        x_1st_derivative[j,i] = (phi[j,i+1] - phi[j,i])/1
                    else: # Style[j,i+1]==0 
                        x_1st_derivative[j,i] = 0
                    # second derivative in x dir
                    if (Style[j,i+1] and Style[j,i+2]):
                        # forward ==> 1st order error
                        x_2nd_derivative[j,i] = (phi[j,i+2] - 2*phi[j,i+1] + phi[j,i])/(1**2)
                    else:   # for the nodes between two solid nodes, and for those node that dont have enough data to calculate forward and backward second derivatives
                        x_2nd_derivative[j,i] = 0
                
                elif i==nx-1:
                    # first derivative in x dir
                    if Style[j,i-1]: 
                        # backward difference ==> 1st order error
                        x_1st_derivative[j,i] = (phi[j,i] - phi[j,i-1])/1
                    else: # Style[j,i-1]==0
                        x_1st_derivative[j,i] = 0
                    # second derivative in x dir
                    if (Style[j,i-1] and Style[j,i-2]):
                        # backward ==> 1st order
                        x_2nd_derivative[j,i] = (phi[j,i] - 2*phi[j,i-1] + phi[j,i-2])/1**2
                    else:   # for the nodes between two solid nodes, and for those node that dont have enough data to calculate forward and backward second derivatives
                        x_2nd_derivative[j,i] = 0
                    
                                    
                # first derivative in y dir
                if Style[j-1,i] and Style[j+1,i]:
                    # central difference
                    y_1st_derivative[j,i] = (phi[j-1,i] - phi[j+1,i])/2
                elif Style[j-1,i] and Style[j+1,i]==0:
                    # forward difference
                    y_1st_derivative[j,i] = (phi[j-1,i] - phi[j,i])/1
                elif Style[j-1,i]==0 and Style[j+1,i]:
                    # backward difference
                    y_1st_derivative[j,i] = (phi[j,i]-phi[j+1,i])/1
                else: # Style[j-1,i]==0 and Style[j+1,i]
                    y_1st_derivative[j,i] = 0
                    
                    
                # second derivative in y dir
                if Style[j-1,i] and Style[j+1,i]:
                    # central
                    y_2nd_derivative[j,i] = (phi[j-1,i] - 2*phi[j,i] + phi[j+1,i]) / (1**2)
                elif (j!=1) and (Style[j-1,i] and Style[j-2,i]) and Style[j+1,i]==0:
                    # forward
                    y_2nd_derivative[j,i] = (phi[j-2,i] - 2*phi[j-1,i] + phi[j,i])/(1**2)
                elif (j!=ny-2) and Style[j-1,i]==0 and (Style[j+1,i] and Style[j+2,i]):
                    # backward
                    y_2nd_derivative[j,i] = (phi[j,i] - 2*phi[j+1,i] + phi[j+2,i])/(1**2)
                else:   # for the nodes between two solid nodes, and for those node that dont have enough data to calculate forward and backward second derivatives
                    y_2nd_derivative[j,i] = 0
    
                    
    return(x_1st_derivative,y_1st_derivative,x_2nd_derivative,y_2nd_derivative)

def phi_1st_2nd_SimplewithSolid(phi,Style):
    """
    this method takes into acount the phase field of solid nodes for dericvatives calculation.
    the derivative calculation method is simple central difference for interior nodes, forward difference for inlet and backward difference for outlet nodes
    """
    ny = len(Style)
    nx = len(Style[0])
    c=dt=1
    
    x_1st_derivative = zeros((ny,nx), dtype=double) * 100
    y_1st_derivative = zeros((ny,nx), dtype=double) * 100
    x_2nd_derivative = zeros((ny,nx), dtype=double) * 100 
    y_2nd_derivative = zeros((ny,nx), dtype=double) * 100
    
    # first derivative phase-field' on every lattice node
    # x direction derivatives
    x_1st_derivative[:,1:nx-1] = (phi[:,2:nx] - phi[:,0:nx-2])/(2*c*dt) # Central Difference in x-dir
    x_1st_derivative[:,0] = (phi[:,1] - phi[:,0])/(1*c*dt)              # left boundary --> Forward difference
    x_1st_derivative[:,nx-1] = (phi[:,nx-1] - phi[:,nx-2])/(1*c*dt)     # right boundary --> Backward 
    # y direction derivatives
    y_1st_derivative[1:ny-1,:] = (phi[0:ny-2,:] - phi[2:ny,:])/(2*c*dt)   # Central in y-dir
    y_1st_derivative[0,:] = (phi[0,:] - phi[1,:])/(1*c*dt)              # upper boundary --> Backward
    y_1st_derivative[ny-1,:] = (phi[ny-2,:] - phi[ny-1,:])/(1*c*dt)       # lower boundary --> Forward    
    
    # second derivative of 'phase-field' on every lattice node
    # x direction drivatives
    x_2nd_derivative[:,1:nx-1] = (phi[:,2:nx] - 2*phi[:,1:nx-1] + phi[:,0:nx-2])/(c*dt)**2   # Central --> (phi1-2phi0+phi3)/dx**2
    x_2nd_derivative[:,0] = (phi[:,2] - 2*phi[:,1] + phi[:,0]) / (c*dt)**2                   # left boundary --> forward --> (phi2-2*ph1+phi0)/dx**2
    x_2nd_derivative[:,nx-1] = (phi[:,nx-1] -  2*phi[:,nx-2] + phi[:,nx-3])/(c*dt)**2        # right boundary --> backward   
    # y direction derivatives
    y_2nd_derivative [1:ny-1,:] = (phi[0:ny-2,:] - 2*phi[1:ny-1,:] + phi[2:ny,:])/(c*dt)**2   # central
    y_2nd_derivative [0,:] = (phi[0,:] - 2*phi[1,:] + phi[2,:]) / (c*dt)**2                   # upper boundary --> backward
    y_2nd_derivative [ny-1,:] = (phi[ny-3,:] - 2*phi[ny-2,:] - phi[ny-1,:]) / (c*dt)**2       # lower boundary --> forward
    
    return x_1st_derivative,y_1st_derivative,x_2nd_derivative,y_2nd_derivative

def phi_1st_2nd_isotropic(phi,Style):
    """
    
    """
    
    ny = len(Style)
    nx = len(Style[0])
    c=dt=1
    Cs2=1./3
    
    x_1st_derivative = zeros((ny,nx), dtype=double)
    y_1st_derivative = zeros((ny,nx), dtype=double)
#    x_2nd_derivative = zeros((ny,nx)) * 100 
#    y_2nd_derivative = zeros((ny,nx)) * 100
    laplacian_phi = zeros((ny,nx), dtype=double)
    
    #eq 34 and 35 Ref. 2 
    for j in range(1,ny-1): 
        for i in range(1,nx-1):     # excluding inlet and outlet - we implement forward and backward difference on these nodes
            if Style[j,i]:  # if the current node is a fluid node
                x_1st_derivative[j,i] = (c/Cs2) * ( cx[1]*w[1]*phi[j-cy[1],i+cx[1]] + cx[2]*w[2]*phi[j-cy[2],i+cx[2]] + \
                                                    cx[3]*w[3]*phi[j-cy[3],i+cx[3]] + cx[4]*w[4]*phi[j-cy[4],i+cx[4]] + cx[5]*w[5]*phi[j-cy[5],i+cx[5]] + \
                                                    cx[6]*w[6]*phi[j-cy[6],i+cx[6]] + cx[7]*w[7]*phi[j-cy[7],i+cx[7]] + cx[8]*w[8]*phi[j-cy[8],i+cx[8]] )
                
                y_1st_derivative[j,i] = (c/Cs2) * ( cy[1]*w[1]*phi[j-cy[1],i+cx[1]] + cy[2]*w[2]*phi[j-cy[2],i+cx[2]] + \
                                                    cy[3]*w[3]*phi[j-cy[3],i+cx[3]] + cy[4]*w[4]*phi[j-cy[4],i+cx[4]] + cy[5]*w[5]*phi[j-cy[5],i+cx[5]] + \
                                                    cy[6]*w[6]*phi[j-cy[6],i+cx[6]] + cy[7]*w[7]*phi[j-cy[7],i+cx[7]] + cy[8]*w[8]*phi[j-cy[8],i+cx[8]])
                   
                laplacian_phi[j,i] = (2/Cs2) * ( w[1]*( phi[j-cy[1],i+cx[1]] - phi[j,i] ) + w[2]*( phi[j-cy[2],i+cx[2]] - phi[j,i] ) + \
                                            w[3]*( phi[j-cy[3],i+cx[3]] - phi[j,i] ) + w[4]*( phi[j-cy[4],i+cx[4]] - phi[j,i] ) + \
                                            w[5]*( phi[j-cy[5],i+cx[5]] - phi[j,i] ) + w[6]*( phi[j-cy[6],i+cx[6]] - phi[j,i] ) + \
                                            w[7]*( phi[j-cy[7],i+cx[7]] - phi[j,i] ) + w[8]*( phi[j-cy[8],i+cx[8]] - phi[j,i] ) )
    
    # top and bottom boundaries: (Nonslip B.C.) it is not important at all in solid nodes
    x_1st_derivative[0,1:nx-2]=(phi[0,2:nx-1] - phi[0,0:nx-3])/3.0
    y_1st_derivative[0,:]=0
    
    x_1st_derivative[ny-1,1:nx-2]=(phi[ny-1,2:nx-1] - phi[ny-1,0:nx-3])/3.0
    y_1st_derivative[ny-1,:]=0
    
    
    # inlet and outlet for Periodic B.C: # I do not think it really matters
    #inlet:
    #x_1st_derivative[1:ny-1,0] = (phi[1:ny-1,1] - phi[1:ny-1,-1])/3.0 + (phi[0:ny-2,1] - phi[2:ny,-1] + phi[2:ny,1] - phi[0:ny-2,-1])/12.0
    #y_1st_derivative[1:ny-1,0] = (phi[0:ny-2,0] - phi[2:ny,0])/3.0 + (phi[0:ny-2,1] - phi[2:ny,-1] + phi[0:ny-2,-1] - phi[2:ny,1] )/12.0
    #x_1st_derivative[0,0] = (phi[0,1] - phi[0,-1])/3.0 + ( phi[1,1] - phi[1,-1])/12.0
    #y_1st_derivative[0,0] = 0 #???
    
    x_1st_derivative[:,0] = x_1st_derivative[:,nx-2]
    y_1st_derivative[:,0] = y_1st_derivative[:,nx-2]
    
    x_1st_derivative[:,nx-1] = x_1st_derivative[:,1]
    y_1st_derivative[:,nx-1] = y_1st_derivative[:,1]
    
    laplacian_phi[:,0] = laplacian_phi[:,nx-2]
    laplacian_phi[:,nx-1] = laplacian_phi[:,1]
    # calculating inlet and aoutlet values (using forward and backward for x-derivative , and , central difference for y-derivative)
    #x_1st_derivative[1:ny-1,0] = (phi[1:ny-1,1] - phi[1:ny-1,0])/(1*c*dt)              # left boundary --> Forward difference
    #x_1st_derivative[1:ny-1,nx-1] = (phi[1:ny-1,nx-1] - phi[1:ny-1,nx-2])/(1*c*dt)     # right boundary --> Backward 
    # y direction derivatives
    #y_1st_derivative[1:ny-1,0] = (phi[0:ny-2,0] - phi[2:ny,0])/(2*c*dt)   # central # excluding upper-most and lower-most solid nodes
    #y_1st_derivative[1:ny-1,nx-1] = (phi[0:ny-2,nx-1] - phi[2:ny,nx-1])/(2*c*dt)   # central # excluding upper-most and lower-most solid nodes
    
    # second derivative of 'phase-field' on every lattice node
    # x direction drivatives
#    phi_x_2nd_inlet = (phi[1:ny-1,2] - 2*phi[1:ny-1,1] + phi[1:ny-1,0]) / (c*dt)**2                   # left boundary --> forward --> (phi2-2*ph1+phi0)/dx**2
#    phi_x_2nd_outlet = (phi[1:ny-1,nx-1] -  2*phi[1:ny-1,nx-2] + phi[1:ny-1,nx-3])/(c*dt)**2        # right boundary --> backward   
    # y direction derivatives
#    phi_y_2nd_inlet = (phi[0:ny-2,0] - 2*phi[1:ny-1,0] + phi[2:ny,0])/(c*dt)**2   # left boundary --> central
#    phi_y_2nd_outlet = (phi[0:ny-2,nx-1] - 2*phi[1:ny-1,nx-1] + phi[2:ny,nx-1])/(c*dt)**2   # right boundary --> central
    # calculating laplacian
#    laplacian_phi[1:ny-1,0] = phi_x_2nd_inlet + phi_y_2nd_inlet
#    laplacian_phi[1:ny-1,nx-1] = phi_x_2nd_outlet + phi_y_2nd_outlet
    return x_1st_derivative,y_1st_derivative,laplacian_phi

#def phi_1st_2nd_derivative_isotropic(phi,Style) 
    

#x_solid_boundary,y_solid_boundary,neighbour_arrangement = normal_vector().boundary_finder(Style)
#solid_normal_vectors = normal_vector().normal_bank(neighbour_arrangement)

#subplots(9,10)

#for i in range(450,540):
#    subplot(9,10,i+1 - 450)
#    IM = [[neighbour_arrangement[i,6],neighbour_arrangement[i,2],neighbour_arrangement[i,5]],
#          [neighbour_arrangement[i,3],             0            ,neighbour_arrangement[i,1]],
#          [neighbour_arrangement[i,7],neighbour_arrangement[i,4],neighbour_arrangement[i,8]]]
#    imshow(IM, cmap='gray'); scatter([solid_normal_vectors[i,0]+1],[-solid_normal_vectors[i,1]+1])
    
    