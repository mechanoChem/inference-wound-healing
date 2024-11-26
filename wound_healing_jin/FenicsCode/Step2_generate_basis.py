from ufl import *
from dolfin import *
import numpy as np
import os

#import h5py as h5


#0.4375 x 0.5 x 0.5
micron_per_pix = 3.45
mag = 10
pix_stride = 100
zoom_um_per_pix = micron_per_pix/mag
zoom_um_per_pix = 0.5           #Overriding this value based on Patrick's email June 29th 
len_stride = pix_stride * zoom_um_per_pix
#each bin is 50 um in other data
len_stride = 0.05
len_stride = 50

oneD_flag=True
halfDataFlag = False
if oneD_flag:
    if halfDataFlag: 
    	lastPixel = 10	#For half length in 1D
    	numXPoints = 38
    else:
	    lastPixel = 19  #For full length in 1D
	    numXPoints = 37
else: #Need to write mesh refinement to do high numXPoints
    if halfDataFlag: 
    	lastPixel = 10	#For half length in 1D
    	numXPoints = 10
    else:
	    lastPixel = 19  #For full length in 1D
	    numXPoints = 19


point0=Point(0,0)
point1=Point(37*len_stride,lastPixel*len_stride)

#mesh = RectangleMesh(MPI.comm_world,point0,point1,60, 19)
mesh=Mesh()
sigma=3
rolling_win=1
total_densities = 6
dense_vec = [20000, 18000, 16000, 14000, 12000, 10000]

for ii in range(total_densities):

    preprocDataLoc = '../results/PreProcess/density' + str(dense_vec[ii]) + '/'
    
    if oneD_flag:
      hdf5=HDF5File(MPI.comm_world, preprocDataLoc + 'cell_density_1D_'+str(sigma)+'_'+str(sigma)+'_rolling_win'+str(rolling_win)+'_refine4'+'.h5','r')
      hdf5.read(mesh,'/mesh',False)
    else:
      hdf5=HDF5File(MPI.comm_world, '../results/PreProcess/cell_density_'+str(sigma)+'_'+str(sigma)+'_rolling_win'+str(rolling_win)+'.h5','r')
      hdf5.read(mesh,'/mesh',False)

    #if halfDataFlag ==True:
    #  boundaries = MeshFunction("size_t", mesh, mesh.topology().dim()-1, 0)
    #  BC2.mark(boundaries, 1)
    #  ds=Measure('ds',subdomain_data = boundaries, metadata={'quadrature_degree': 2})

    V=FunctionSpace(mesh, "Lagrange", 1)
    W=FunctionSpace(mesh, "DG", 0)

    C = Function(V)
    erk = Function(V) 
    akt = Function(V) 

    print('total_DOF=',C.vector()[:].size)

    w= TestFunction(V)

    y_index=1
    v=Constant((0.0,1.0))
    BC1 =  CompiledSubDomain("near(x[1], side) && on_boundary", side = 0)
    BC2 =  CompiledSubDomain("near(x[1], side) && on_boundary", side = lastPixel*len_stride)
    if oneD_flag:
      y_index=0
      v=Constant((1.0,))
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


    x=SpatialCoordinate(mesh)
    y_m=lastPixel*len_stride/2
    y_l=lastPixel*len_stride





    print(len(grad_w))
    print(len(grad_C))
    # y_coord=x[0]
    # y_tilda = (y_coord-y_m)/y_l

    # if oneD_flag:
      # y_coord=x
      # y_tilda = (y_coord-y_m)/y_l      
      # v_=project(conditional(gt(x[0], y_m), -1, 1), W)
      # v=Expression('conditional(gt(x[0], ym), -1, 1)',y_m)
    v =  -grad_C/(1e-10 + pow(dot(grad_C, grad_C), 0.5))

    def assemble_R(basis_id):
    #TODO where is dx defined?
    #Diffusion   
      if basis_id==0:
        R = -1*inner(grad_w,grad_C)*dx
      elif basis_id==1:
        R = -C*inner(grad_w,grad_C)*dx
      elif basis_id==2:	
        R = -C*C*inner(grad_w,grad_C)*dx

            
    #Advection   
      elif basis_id==3:
        R = C*inner(grad_w,v)*dx
      elif basis_id==4:
        R = C*inner(grad_w,v)*C*dx
      elif basis_id==5:
        R = C*inner(grad_w,v)*C*C*dx


    #Reaction
##      elif basis_id==6:  	
##        R = w*dx
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


    total_num_step=5
    basisDirDensity = 'density' + str(dense_vec[ii])
    if oneD_flag:
      directory = '../results/basis/Physics_Based/'+basisDirDensity + '/basis_1D_'+str(sigma)+'_'+str(sigma)+'_rolling_win'+str(rolling_win)
    else: 
      directory = '../results/basis/Physics_Based/basis_2D_'+str(sigma)+'_'+str(sigma)+'_rolling_win'+str(rolling_win)


    if not os.path.exists(directory):
      os.makedirs(directory)
      print('New directory created: %s'%directory)
      
    NumBasis=8
    for step in range(total_num_step):
      print('step=',step)
      hdf5.read(C,'density_'+str(step))

      
      
      basis=np.column_stack([assemble_R(basis_id) for basis_id in range(NumBasis+1)]) #An extra one to keep track of C
      # 
      if oneD_flag:
        filename = directory+'/basis_step_'+str(step)+'_refine4'+'.dat'
      else: 
        filename = directory+'/basis_step_'+str(step)+'.dat'
      
      
      np.savetxt(filename,basis)








