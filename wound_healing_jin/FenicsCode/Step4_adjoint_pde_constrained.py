from ufl import *
from dolfin import *
import numpy as np
#import h5py as h5
from dolfin_adjoint import *
set_log_active(False)
import os
import sys

def forward_model(data_filename, theta, learn_time_steps, oneD_flag=True, halfDataFlag=False, plotFlag = False):
  
  ################################# Define mesh #################################
 
  mesh = Mesh()
  hdf5=HDF5File(MPI.comm_world, data_filename,'r')
  hdf5.read(mesh,'/mesh',False)  

  ########################### Define Functional spaces ##########################

  V=FunctionSpace(mesh, "Lagrange", 1)
  C = Function(V)
  C_n = Function(V)
  C_BC = Function(V)

  w= TestFunction(V)

  grad_w = grad(w)
  grad_C=grad(C)


  #dxx = dx(metadata={'quadrature_degree': 1})
  dxx = dx

  V_data=FunctionSpace(mesh, "Lagrange", 1)
  C_data_n = Function(V_data)
  C_data_np1 = Function(V_data)
  Error_data_n = Function(V_data)

  grad_C_data_np1=grad(C_data_np1)

  x=SpatialCoordinate(mesh)
  y_m=lastPixel*len_stride/2
  y_l=lastPixel*len_stride

  if oneD_flag:
    y_coord=x[0]
    v=Constant((1.0,))
  else:
    y_coord=x[1]
    v=Constant((0.0,1.0))
  y_tilda = (y_coord-y_m)/y_l

  ############################## Define Boundaries ##############################

  if oneD_flag==True:
    BC1 =  CompiledSubDomain("near(x[0], side) && on_boundary", side = 0)
    BC2 =  CompiledSubDomain("near(x[0], side) && on_boundary", side = lastPixel*len_stride)

  else:
    BC1 =  CompiledSubDomain("near(x[1], side) && on_boundary", side = 0)
    BC2 =  CompiledSubDomain("near(x[1], side) && on_boundary", side = lastPixel*len_stride)



  ################# Define residue (without Boundary condition) #################


  num_basis=len(theta)
  #theta=[]
  #for i in range(num_basis):
  #  theta.append(Constant(0.0))  
  
  basis_pool=[0]*num_basis
  basis_pool[0]=theta[0]*1*inner(grad_w,grad_C)*dxx
  basis_pool[1]=theta[1]*C*inner(grad_w,grad_C)*dxx
  basis_pool[2]=theta[2]*C*C*inner(grad_w,grad_C)*dxx


  diff_supg = theta[0] + theta[1]*C + theta[4]*C*C 
    
  basis_pool[3]=-theta[3]*C*inner(grad_w,v)*1*dxx
  basis_pool[4]=-theta[4]*C*inner(grad_w,v)*C*dxx
  basis_pool[5]=-theta[5]*C*inner(grad_w,v)*C*C*dxx
    
  adv_supg=v*(theta[3] + theta[4]*C + theta[5]*C*C) 

     
  basis_pool[6]=-theta[6]*w*C*dxx 
  basis_pool[7]=-theta[7]*w*C*C*dxx 


  reac_supg= theta[6]*C + theta[7]*C*C
  
  dt=Constant(1)
  R=(C-C_n)*w/dt*dxx
  for i in range(num_basis):
    R+=basis_pool[i]
  
  #SUPG stablilization
  tau_ini=100.
  tau=Constant(tau_ini)
  
  W1 = VectorFunctionSpace(mesh, "CG",1)
  W2 = FunctionSpace(mesh, "DG",0)		
  lap_c = project(div(project(grad(C),W1)),W2)
  R+=tau*len_stride*inner(grad_w,adv_supg)*((C-C_n)/dt -diff_supg*lap_c +div(adv_supg*C)-reac_supg )*dxx  
  
  #Jacobian for nonlinear solver
  J=derivative(R, C)

  # Time-stepping
  t = 0
  results = []
  loss = 0.0

  hdf5.read(C_n,'density_'+str(0))
  C_data_np1.assign(C_n)

  
  
  for step in range(total_num_step):
    print('======== step=%f ======== '%step)
    C_data_n.assign(C_data_np1)


    #Experimental and boundary values
    hdf5.read(C_data_np1,'density_'+str(step+1))

    if reinitializeFlag and step%reinitializeInterval==0: 
      C_n.assign(C_data_n)	
  
    #for i in range(num_basis):
    #  theta[i].assign(gamma[i])  
  
  
    #For stability issues uniformly divide the time steps into smaller pieces
    for mini_step in range(numSmallStep):
      print('mini step =  ', mini_step+1)
    
      C_BC.assign( ((mini_step+1)*C_data_np1 + (numSmallStep - mini_step - 1)*C_data_n)/numSmallStep )
      dt.assign((time[step+1]-time[step])/numSmallStep)
	
      bcl_1 = DirichletBC(V, C_BC, BC1)	#C=C_data on bc1
      if halfDataFlag:
        bcs = [bcl_1]	#No flux on bc2
        #vel=theta[4]*grad_C + theta[5]*grad_erk + theta[6]*grad_akt 
      else: 
        bcl_2 = DirichletBC(V, C_BC, BC2)   #C=C_data on bc2
        bcs=[bcl_1, bcl_2]
      
      problem = NonlinearVariationalProblem(R, C,bcs,J)
      solver = NonlinearVariationalSolver(problem)
      prm = solver.parameters
      prm["newton_solver"]["absolute_tolerance"] = 1E-8
      prm["newton_solver"]["relative_tolerance"] = 1E-9
      prm["newton_solver"]["maximum_iterations"] = 300
      prm['newton_solver']["error_on_nonconvergence"] = False

      solver.solve()
      C_n.assign(C)
    if step in learn_time_steps: 
      loss += loss_func(C_data_np1, C, grad_C_data_np1, grad_C)
    results.append(C.copy())
  if plotFlag:
    plot_C = Function(V)
    file_plotC = XDMFFile(MPI.comm_world,'../results/Adjoint_gamma_matrix/Physics_Based/final_density.xdmf') 
    for step in range(total_num_step):
      plot_C.assign(results[step])
      file_plotC.write(plot_C,step+1)
    file_plotC.close()
  return loss, results
  
def loss_func(C_exp, C_sim, grad_C_exp, grad_C_sim):
  factor = 0
  loss = assemble( ((C_exp - C_sim)**2 + factor*inner(grad_C_exp - grad_C_sim, grad_C_exp - grad_C_sim))  *dx)
  #print(assemble( (C_exp - C_sim)**2*dx))
  #print(assemble( inner(grad_C_exp - grad_C_sim, grad_C_exp - grad_C_sim)*dx))
  return loss

def iter_cb(m):
  print ("m = ", m)

def eval_cb(j, m):
  print ("j = %f, m = %f." % (j, m))

def derivative_cb(j, dj, m):
  print("j = %f, dj = %s, m = %s." %(j, [float(k) for k in dj], [float(k) for k in m]))   
###############################################################################
################################ Main Code ####################################
###############################################################################


################################ Problem Setup ################################

sigma=3
rolling_win=1
F=200000
best_F = F

total_num_step = 4
numSmallStep = 5
learnStepInterval = 1
maxiter=500

learn_time_steps = np.arange(1,total_num_step,learnStepInterval)

oneD_flag = True
halfDataFlag = False
#Set newMeshFlag = True if do not want to read erk and akt data
reinitializeFlag = False
reinitializeInterval = 500
if oneD_flag:
    if halfDataFlag: 
    	lastPixel = 10	#For half length in 1D
    	numXPoints = 37
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

micron_per_pix = 3.45
mag = 10
pix_stride = 100
zoom_um_per_pix = micron_per_pix/mag
len_stride = pix_stride * zoom_um_per_pix
len_stride = 50

all_densities = [20000, 18000, 16000, 14000, 12000, 10000]
# all_densities = [18000]
for ii, density in enumerate(all_densities):
  if oneD_flag:
    data_filename = '../results/PreProcess/density' + str(density) + '/cell_density_1D_'+str(sigma)+'_'+str(sigma)+'_rolling_win'+str(rolling_win)+'_refine4'+'.h5'
  else:
    data_filename = '../results/PreProcess/cell_density_'+str(sigma)+'_'+str(sigma)+'_rolling_win'+str(rolling_win)+'.h5'
    
  time=np.loadtxt('../data/time.dat',delimiter=",")
  #name_VSI='1D_Group_'+str(sigma)+'_'+str(sigma)+'_rolling_win'+str(rolling_win)+'_F'+str(F)+'_refine4'
  #gamma_matrix_VSI=np.loadtxt('../results/VSI_gamma_matrix/Advection_spatial/gamma_'+name_VSI+'.dat')

  if oneD_flag:
      folder_ID = 'Physics_Based_Time_Independent_1D/density' + str(density) + '/'
      file_ID = 'gamma_history_Group_'+str(sigma)+'_'+str(sigma)+'_rolling_win'+str(rolling_win)+'_F'+str(F)+'_refine4'+'.dat'
  else:
      folder_ID = 'Physics_Based_Time_Independent_2D/'
      file_ID = 'gamma_Group_'+str(sigma)+'_'+str(sigma)+'_rolling_win'+str(rolling_win)+'_F'+str(F)+'.dat'

  filename_VSI = '../results/VSI_gamma_matrix/'+folder_ID+file_ID
  gamma_matrix_VSI=np.loadtxt(filename_VSI)
  
  gamma = np.array(gamma_matrix_VSI[:,5].transpose())
  print('####')
  print(gamma)
  # gamma[0] = 100
  #set control
  num_basis = 8

  control_index=[]
  theta=[]
  for i in range(num_basis):
    theta.append(Constant(gamma[i]))  
    if(np.abs(gamma[i])>1.0e-12):
      control_index.extend([i])			#To make sure only non-zero gamma are fine tuned


  #Set forward


  loss, learned_data = forward_model(data_filename, theta, learn_time_steps, oneD_flag, halfDataFlag)

     
  print(f'ii = {ii}, Density = {density}')
  print(control_index,gamma[control_index])

  #do adjoint
  control_parameter=[Control(theta[i]) for i in control_index]  


  #Setting up bounds
  diffu_index=[i for i in range(3)]
  adv_index=[i for i in range(3,6)]
  reac_index=[i for i in range(6,8)]
          

  bounds=np.zeros((2,len(control_index)))
  bounds[0,:]=-np.inf
  bounds[1,:]=np.inf

  # for i in range(len(control_index)):
  #   index = control_index[i]
  #   if index in diffu_index:
  #     bounds[0,i] = 1e-3*gamma[index]
  #     bounds[1,i] = 1e+3*gamma[index]
  #   else:
  #     if gamma[index]>0:   
  #         bounds[0,i] = -1e-3*(gamma[index])
  #         bounds[1,i] =  1e+3*(gamma[index])
  #     else:
  #         bounds[0,i] = 1e+3*(gamma[index])
  #         bounds[1,i] = -1e-3*(gamma[index])
  for i in range(len(control_index)):
    index = control_index[i]
    if index in diffu_index:
      bounds[0,i] = 1e-6
      bounds[1,i] = 1e+10
  #control_index_diff = []
  #control_index_adv = []
  #for i in range(len(control_index)):
  #  index = control_index[i]
  #  if index in diffu_index: 
  #    bounds[0,i] = 1e-6
  #    #control_index_diff.append(i)
  #  if index in adv_index:
  #    bounds[1,i] = 1e-6
  #    #control_index_adv.append(i)
          
    
  method = "L-BFGS-B"
  reduced_functional = ReducedFunctional(loss, control_parameter, derivative_cb_post=derivative_cb)
  results_opt = minimize(reduced_functional, method = method, bounds=bounds,tol=1.0e-9, options = {'ftol':1.0e-11,'maxiter':maxiter, 'disp': True},callback = iter_cb)
  converge_flag=True 

  tem=gamma
  if len(control_index)==1:
    tem[control_index[i]]=results_opt
  else:
    for i in range(len(control_index)):
      tem[control_index[i]]=results_opt[i]
  print(tem)

  gamma_matrix_all = [tem for i in range(total_num_step)]




  if oneD_flag:
      directory = '../results/Adjoint_gamma_matrix/Physics_Based_Time_Independent_1D/density'+str(density)
  else:
      directory = '../results/Adjoint_gamma_matrix/Physics_Based_Time_Independent_2D'

  if not os.path.exists(directory):
    os.makedirs(directory)
    print('New directory created: %s'%directory)
  #np.savetxt(directory + '/gamma_1D_Group_'+str(sigma)+'_'+str(sigma)+'_rolling_win'+str(rolling_win)+'_refine4'+'.dat',gamma_matrix_all)

  if oneD_flag:
      np.savetxt(directory + '/gamma_Group_'+str(sigma)+'_'+str(sigma)+'_rolling_win'+str(rolling_win)+'_F'+str(best_F)+'_refine4'+'.dat',gamma_matrix_all)
      #np.savetxt(directory + '/loss_Group_'+str(sigma)+'_'+str(sigma)+'_rolling_win'+str(rolling_win)+'_F'+str(best_F)+'_refine4'+'.dat',loss_all)
  else:
      np.savetxt(directory + '/gamma_Group_'+str(sigma)+'_'+str(sigma)+'_rolling_win'+str(rolling_win)+'_F'+str(best_F)+'.dat',gamma_matrix_all)
      #np.savetxt(directory + '/loss_Group_'+str(sigma)+'_'+str(sigma)+'_rolling_win'+str(rolling_win)+'_F'+str(best_F)+'.dat',loss_all)	

  for i in range(num_basis): 
    theta[i].assign(tem[i])

  ##plotFlag = True
  ##loss, learned_data = forward_model(data_filename, theta, learn_time_steps, oneD_flag, halfDataFlag, plotFlag)

