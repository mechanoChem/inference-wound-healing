from ufl import *
from dolfin import *
import numpy as np
import os
import pickle
print('~~~~~~~~~~~~~~~~~~~~~~~~~~~ STARTING GENERATE BASIS ~~~~~~~~~~~~~~~~~~~~~~~~~~~')

with open('vsi_settings.pkl' , 'rb') as settingsFile:
    vps = pickle.load( settingsFile)    
# Create mesh and define function space
#0.4375 x 0.5 x 0.5
micron_per_pix = 3.45
mag = 10
pix_stride = 100
zoom_um_per_pix = micron_per_pix/mag
zoom_um_per_pix = 0.5           #Overriding this value based on Patrick's email June 29th 
len_stride = pix_stride * zoom_um_per_pix

len_stride = vps.len_stride
numXPoints = vps.numXPoints
lastPixel = vps.lastPixel
oneD_flag = vps.oneD_flag
halfDataFlag = vps.halfDataFlag
wells = vps.wells
nWells = vps.nWells
time_gap = vps.time_gap
total_time = vps.total_time
sigma = vps.sigma
rolling_win = vps.rolling_win
total_num_step = total_time
NumBasis = vps.total_num_basis

for w in range(nWells):
    wellCur = wells[w]

    point0=Point(0,0)
    point1=Point(numXPoints*len_stride,lastPixel*len_stride)

    mesh=Mesh()
    hdf5path = vps.dataSavePathPreProc.format(wellCur)
    hdf5name = vps.dataSaveNameRWh5.format(sigma, sigma, rolling_win)
    hdf5 = HDF5File(MPI.comm_world, hdf5path+hdf5name, 'r')
    hdf5.read(mesh, '/mesh', False)


    V=FunctionSpace(mesh, "Lagrange", 1)

    C = Function(V)
    #erk = Function(V) 
    #akt = Function(V) 

    print('total_DOF=',C.vector()[:].size)

    w= TestFunction(V)

    y_index=1
    #old
    #v=Constant((1.0,0.0)) # original: v=Constant((0.0, 1.0)) But I think I need to do the opposite bc my data is oriented the other way
    #new
    #C_n = C at previous time step - using C_n  = C because it's easier
    C_n = C
    grad_C_n = grad(C_n)
    v = -grad_C_n/(1e-6 + pow(dot(grad_C_n, grad_C_n), 0.5))

    BC1 =  CompiledSubDomain("near(x[1], side) && on_boundary", side = 0)
    BC2 =  CompiledSubDomain("near(x[1], side) && on_boundary", side = lastPixel*len_stride)
    if oneD_flag:
      y_index=0
      #v=Constant((1.0,))
     
      BC1 =  CompiledSubDomain("near(x[0], side) && on_boundary", side = 0)
      BC2 =  CompiledSubDomain("near(x[0], side) && on_boundary", side = lastPixel*len_stride)

    
    bcl_1 = DirichletBC(V, 0, BC1)	#C=C_data on bc1
    if halfDataFlag:
      bcs = [bcl_1]	#No flux on bc2
    else: 
      bcl_2 = DirichletBC(V, 0, BC2)   #C=C_data on bc2
      bcs=[bcl_1, bcl_2]

    grad_C = grad(C)
    grad_w = grad(w)
    #grad_erk = grad(erk)
    #grad_akt = grad(akt)

    x=SpatialCoordinate(mesh)
    y_m=lastPixel*len_stride/2
    y_l=lastPixel*len_stride


    if oneD_flag:
      v=Constant((1.0,))

    print(len(grad_w))
    print(len(grad_C))
    y_coord=x[0]
    y_tilda = (y_coord-y_m)/y_l
    def assemble_R(basis_id):
    #Diffusion   
      if basis_id==0:
        R = -1*inner(grad_w,grad_C)*dx
      elif basis_id==1:
        R = -C*inner(grad_w,grad_C)*dx
      elif basis_id==2:	
        R = -C*C*inner(grad_w,grad_C)*dx

        
    #Advection   
      elif basis_id==3:
        R = C*inner(grad_w,v)*1*dx
      elif basis_id==4:
        R = C*inner(grad_w,v)*C*dx
      elif basis_id==5:
        R = C*inner(grad_w,v)*C*C*dx
      

    #Reaction
      elif basis_id==6:  
        R = w*C*dx
      elif basis_id==7:  
        R = w*C*C*dx


    #For keeping track of density	
      elif basis_id==8:
        R = C*w*dx	
        
      R_=assemble(R)
      for bc in bcs:
        bc.apply(R_)
      R_value=R_.get_local()
      
      return R_value


    saveDirectory = vps.basisSavePath.format(wellCur, sigma, sigma, rolling_win)
    print(f'basisSavePath is {vps.basisSavePath}')
    print(f'save directory is : {saveDirectory}')
    if not os.path.exists(saveDirectory):
      os.makedirs(saveDirectory)
      print(f'New directory created: {saveDirectory}')
      
    for step in range(total_num_step):
      print(f'well = {wellCur}, step = {step}')
      hdf5.read(C,'density_'+str(step))
      #hdf5.read(erk,'erk_'+str(step))
      #hdf5.read(akt,'akt_'+str(step))
      
      
      basis=np.column_stack([assemble_R(basis_id) for basis_id in range(NumBasis+1)]) #An extra one to keep track of C
      # 
      basisFilename = vps.basisSaveName.format(step)
      

      np.savetxt(saveDirectory+basisFilename,basis)






print('~~~~~~~~~~~~~~~~~~~~~~~~~~~ ENDING GENERATE BASIS ~~~~~~~~~~~~~~~~~~~~~~~~~~~')




