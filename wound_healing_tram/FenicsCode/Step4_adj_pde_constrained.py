from ufl import *
from dolfin import *
import numpy as np
#import h5py as h5
from dolfin_adjoint import *
set_log_active(False)
import os
import sys
import pickle
import time as timelib

testcase_num = 0 #Choose a number 0-6 representing the well group number


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

  diff_supg = theta[0] + theta[1]*C + theta[2]*C*C 
      
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
    # print('======== step=%f ======== '%step)
    C_data_n.assign(C_data_np1)
   
    #Experimental and boundary values
    hdf5.read(C_data_np1,'density_'+str(step+1))
    if reinitializeFlag and step%reinitializeInterval==0: 
      C_n.assign(C_data_n)	
  
    #for i in range(num_basis):
    #  theta[i].assign(gamma[i])  
  
  
    #For stability issues uniformly divide the time steps into smaller pieces
    for mini_step in range(numSmallStep):
      
      # print('mini step =  ', mini_step+1)
    
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
  factor = 1
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
print('~~~~~~~~~~~~~~~~~~~~~~~~~~~ STARTING ADJ OPT ~~~~~~~~~~~~~~~~~~~~~~~~~~~')
adjoint_restart_flag = True

if adjoint_restart_flag:
    print(f'Running testcase {testcase_num} with adjoint restart parameters')
else: 
    print(f'Running testcase {testcase_num} with VSI parameters')
   
# with open('vsi_settings.pkl' , 'rb') as settingsFile:
#     vps = pickle.load( settingsFile)   

# well = 'A03' #SS
sigma = 3
rolling_win = 5
F = 200000
best_F = F
total_num_step = 108
numSmallStep = 2
learnStepInterval = 9
reinitializeInterval=10
maxiter=100

learn_time_steps = np.array([learnStepInterval+k*reinitializeInterval for k in range(total_num_step) if learnStepInterval+k*reinitializeInterval<=total_num_step] )

num_basis = 8
method = 'L-BFGS-B'
converge_flag = True


oneD_flag = False
halfDataFlag = False
#Set newMeshFlag = True if do not want to read erk and akt data
reinitializeFlag = True
reinitializeInterval = 10
lastPixel = 41
numXPoints = 79
len_stride = 96*0.5
# data_path = vps.dataSavePathPreProc.format(well)
# data_file = vps.dataSaveNameRWh5.format(sigma, sigma, rolling_win)
# data_filename = data_path + data_file
  
time=np.loadtxt('../data/time.dat',delimiter=",")
# wells = vps.wells
# nWells = vps.nWells
# wellGroups = vps.wellGroups
# nGroups = vps.nGroups
# print('Wells are:')
# print(wells)
# print(f'There are {nGroups} well groups')
#Iterate through groups
# for gg in range(nGroups):
wellGroups = [['A01', 'B01', 'C01', 'D01'], ['A02', 'B02', 'C02', 'D02'], ['A03', 'B03', 'C03', 'D03'], \
        ['A04', 'B04', 'C04', 'D04'], ['A05', 'B05', 'C05', 'D05'], ['A06', 'B06', 'C06', 'D06']]
wellindex = [0, 1, 2, 3, 4, 5]

gg = wellindex[testcase_num]
groupCur = wellGroups[testcase_num]
loss = [0]*len(groupCur)
print(f'Group: {gg+1}/{6} with wells: {groupCur}')
#Get gamma matrix for this group
filepath_VSI = '../results/VSI_gamma_matrix/{}/Physics_Based_Time_Independent/'.format(f'group{gg+1}_{groupCur[0]}_to_{groupCur[-1]}')
filename_VSI = filepath_VSI + 'gamma_Group_{}_{}_rolling_win{}_F{}.dat'.format(sigma, sigma, rolling_win, F)
filename_VSI_history = filepath_VSI + 'gamma_history_Group_{}_{}_rolling_win{}_F{}.dat'.format(sigma,sigma, rolling_win, F)
gamma_matrix_VSI = np.loadtxt(filename_VSI_history)

#Iterate through bases to optimize
for kk in range(5, 8):

    #Filename to be saved
    adjSaveDirectory = '../results/Adjoint_gamma_matrix/{}/Physics_Based_Time_Independent/'.format(f'group{gg+1}_{groupCur[0]}_to_{groupCur[-1]}_step{kk}')
    if not os.path.exists(adjSaveDirectory):
        os.makedirs(adjSaveDirectory)
        print(f'New dictionary created: {adjSaveDirectory}')
        
    adjSaveName = 'gamma_Group_{}_{}_rolling_win{}_F{}.dat'.format(sigma, sigma, rolling_win, best_F) 
    adjSaveName_temp = 'temp_gamma_Group_{}_{}_rolling_win{}_F{}.dat'.format(sigma, sigma, rolling_win, best_F) 


    gamma = np.array(gamma_matrix_VSI[:,kk].transpose())

    control_index = []
    theta = []
    if adjoint_restart_flag:
      gamma_matrix_VSI_temp = np.loadtxt(adjSaveDirectory+adjSaveName_temp)
      iii=0
      for ii in range(num_basis):
        if(np.abs(gamma[ii])>1.0e-12):
          gamma[ii] = gamma_matrix_VSI_temp[iii]
          iii = iii+1

    for ii in range(num_basis):
        theta.append(Constant(gamma[ii]))
        if(np.abs(gamma[ii])>1.0e-12):
            control_index.extend([ii])

    print(f'for group {gg+1}/{6}, basis {kk}, gamma = {gamma}')

    #Import all wells in this group, calculate loss, and compile cumulative loss
    loss_sum = 0
    for ww, wellCur in enumerate(groupCur):
        data_filepath = '../results/PreProcess/{}/'.format(wellCur)
        data_filename = 'cell_density_{}_{}_rolling_win{}.h5'.format(sigma, sigma, rolling_win)
        loss[ww], learned_data = forward_model(data_filepath+data_filename, theta, learn_time_steps, oneD_flag, halfDataFlag)
        print(f'ww = {ww}, well = {wellCur}, loss = {loss[ww]}')
        loss_sum = loss_sum + loss[ww]
    
    
    #Set up bounds
    control_parameter = [Control(theta[i]) for i in control_index]
    diffu_index = [i for i in range(3)]
    adv_index = [i for i in range(3,6)]
    reac_index = [i for i in range(6, num_basis)]
    bounds=np.zeros((2,len(control_index)))
    bounds[0,:]=-np.inf
    bounds[1,:]=np.inf

    for i in range(len(control_index)):
        index = control_index[i]
        if index in diffu_index:
            bounds[0,i] = 1e-3*gamma[index]
            bounds[0,i] = 0;
            bounds[1,i] = 1e+3*gamma[index]
        else:
        #SS: changed the bounds to \pm 30% correction
            bounds[0,i] = gamma[index]-0.3*np.abs(gamma[index])
            bounds[1,i] = gamma[index]+0.3*np.abs(gamma[index])


    def eval_cb_save(j, m):
        with open(adjSaveDirectory+adjSaveName_temp, "wb") as f:
            f.write(b'\n')
            np.savetxt(f, m)

    #reduced_functional = ReducedFunctional(loss_sum, control_parameter, derivative_cb_post=derivative_cb)
    reduced_functional = ReducedFunctional(loss_sum, control_parameter, eval_cb_post=eval_cb_save)
    
    main_start_optim = timelib.perf_counter()	
    results_opt = minimize(reduced_functional, method = method, bounds=bounds,tol=1.0e-9, options = {'ftol':1.0e-11,'maxiter':maxiter, 'disp': True},callback = iter_cb)
    main_stop_optim = timelib.perf_counter()	
    if(MPI.comm_world.Get_rank() == 0):
        print('===============================================')
        print('Reinitialization after every: %i steps'%(reinitializeInterval))
        print('Main optim time: %f seconds'%(float(main_stop_optim) -float(main_start_optim)))
        print('===============================================')


    tem = gamma
    """
    #Old - why does next 2 lines work? tem[control_index[i]] references previous loop
        if len(control_index)==1:
        tem[control_index[i]]=results_opt
    else:
        for i in range(len(control_index)):
        tem[control_index[i]]=results_opt[i]
    """
    if len(control_index)==1:
        tem[control_index] = results_opt
    else:
        for ii in range(len(control_index)):
            tem[control_index[ii]] = results_opt[ii]
    print(f'After minimize, tem = {tem}')
    #This replicates tem (gamma matrix) for all timesteps, could be better way to do it
    gamma_matrix_all = [tem for ii in range(total_num_step)]
    
    #Save gamma matrix all
    np.savetxt(adjSaveDirectory+adjSaveName, gamma_matrix_all)
    print(f'Done with step {kk+1}/{num_basis} with wells: {groupCur}')
print(f'Done with group: {gg+1}/{6} with wells: {groupCur}')
print('~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~')
print('~~~~~~~~~~~~~~~~~~~~~~~~~~~ ENDING ADJ OPT ~~~~~~~~~~~~~~~~~~~~~~~~~~~')
