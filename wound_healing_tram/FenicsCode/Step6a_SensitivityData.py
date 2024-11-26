print('starting sensitivity code before imports')
from ufl import *
from dolfin import *
import numpy as np
from dolfin_adjoint import *
set_log_active(False)
import os
import sys
import pickle
print('after imports')
"""
pseudocode for fwd model
1. define parameters
2. Open mesh from preprocess location
3. Define functional space
4. Define boundary conditions
5. Define residue (bases are here)
6. Load bases (dependent on run type)
7. Run fwd model for each basis loaded
8. Save fwd model and other files
"""

"""
New code:
define parameters

for gg in groups:
  for ss in steps:
    import basis for group gg, step ss
    for ww in length(groups[gg]):
      Define functional spaces
      Define boundaries
      Define residue
      Import mesh for well ww
      run Fwd model with basis[gg,ss] on init cond ww
      save
"""
print('~~~~~~~~~~~~~~~~~~~~~~~~~~~ STARTING SENSITIVITY CODE ~~~~~~~~~~~~~~~~~~~~~~~~~~~')

with open('vsi_settings.pkl' , 'rb') as settingsFile:
    vps = pickle.load( settingsFile)    
print('after open vps')
    
#Note - this code only runs for 1 well at a time - need to update for one well then wrap a for loop 

################################ Problem Setup ################################

sigma = vps.sigma
rolling_win = vps.rolling_win
F = vps.best_F

total_num_step = vps.total_time
numSmallStep = vps.numSmallSteps

oneD_flag = vps.oneD_flag
temporalParam_flag = vps.temporalParamFlag
halfDataFlag = vps.halfDataFlag
newMeshFlag = vps.newMeshFlag
reinitializeFlag = vps.reinitializeFlag
reinitializeInterval = vps.reinitializeIntervalFwd

lastPixel = vps.lastPixel
numXPoints = vps.numXpoints                
len_stride = vps.len_stride
nWells = vps.nWells
wells = vps.wells
wellGroups = vps.wellGroups
nGroups = vps.nGroups
nSteps = vps.endBasis-vps.startBasis

#For loss estimation
learnStepInterval = vps.learnStepInterval
learn_time_steps = np.arange(1,total_num_step,learnStepInterval)

counter = 0
# Iterate through each group (collection of wells with identical conditions)
print(f'about to start iterating')
print(f'nGroups: {nGroups}')
for gg in range(nGroups): 
  print(f'inside first loop: gg = {gg}')
  groupCur = wellGroups[gg]
  loss_list = []
  #Select the N-1 timestep (when there are 2 bases remaining)
  print(f'about to start second loop: [nSteps-2] = {nSteps-2}')
  for ss in [nSteps-3]:#range(nSteps): #SS: Fixing stepwise step for sensititity plot to last-1
   
    print(f'Group: {gg+1}/{nGroups} ---- step: {ss+1}/{nSteps} ---- Opening gamma matrix')
    stepCur = ss+vps.startBasis
    time = np.loadtxt('../data/time.dat', delimiter = ",")
    plot_tag = vps.adjPlotTag
    #Build names for adjoint gamma location and open
    adjSaveDirectory = vps.adjSavePath.format(f'group{gg+1}_{groupCur[0]}_to_{groupCur[-1]}_step{stepCur}')
    adjSaveName = vps.adjSaveNameGamma.format(sigma, sigma, rolling_win, F)
    gamma_matrix_original = np.loadtxt(adjSaveDirectory+adjSaveName)
    temp = f'sensdata_Group_{sigma}_{sigma}_rolling_win{rolling_win}_F{F}.dat'

    sensitivity_data_filename = adjSaveDirectory+temp
    param1_filename = adjSaveDirectory+f'senseparam1_Group_{sigma}_{sigma}_rolling_win{rolling_win}_F{F}.dat'
    param2_filename = adjSaveDirectory+f'senseparam2_Group_{sigma}_{sigma}_rolling_win{rolling_win}_F{F}.dat'
    print(f'adjSaveName: {adjSaveName}')
    print(f'temp: {temp}')
    print(f'sense data filename: {sensitivity_data_filename}')
    
    print('~~~~~~~~~~~')
    print(f'gamma size = {gamma_matrix_original.shape}')
    print(f'gamma row 0 = {gamma_matrix_original[0,:]}')
    print('~~~~~~~~~~~')
     
    #SS: mean param values
    #D0_mean = 6.8
    #C1_mean = 0.05
    #C2_mean = -34.0
    #TODO: update to vary each around the Adjoint selected values
    #TODO: update to vary the number of points in each (use linspace or similar)
    
    #Need to be equally placed for contour function
    D0_range = [1,6,11,16,21]
    C1_range = [0.01,0.03,0.05,0.07,0.09]
    C2_range = [10.,20.,30.,40.,50.]

    D0_range = [6, 11, 16]
    C1_range = [0.03, 0.05, 0.07]
    C2_range = [20., 30., 40.]
    D0_range = [6, 6]
    C1_range = [0.03, 0.03]
    #SS: Run combination: 
    #SS: Currently we will only vary D0 and C1 but change run_combo for more options
    
    # 0: D0, C1
    # 1: D0, C2
    # 2: C1, C2
    run_combo = 0
    if run_combo==0:
      param_pos = [0,6]
      param_vals1 = D0_range
      param_vals2 = C1_range
    elif run_combo==1:
      param_pos = [0,7]
      param_vals1 = D0_range
      param_vals2 = C2_range
    elif run_combo==2:
      param_pos = [6,7]
      param_vals1 = C1_range
      param_vals2 = C2_range
    else:  
      exit()
    
    with open(param1_filename, 'wb') as f:
      np.save(f, param_vals1)
    with open(param2_filename, 'wb') as f:
      np.save(f, param_vals2)
    print(f'starting param iteration for run_combo = {run_combo}')
    
    for param1 in param_vals1: 
      for param2 in param_vals2:
        
        print(f'param 1 = {param1}, param2 = {param2}')
        loss = 0.0 #Initialize zero loss
        
        #Define parameters for this run. 
        gamma_matrix = gamma_matrix_original
        gamma_matrix[:,param_pos[0]] = param1
        gamma_matrix[:,param_pos[1]] = param2

        #Iterate through all wells in group
        for ww in range(len(groupCur)):
            print(f'Applying basis to well {groupCur[ww]}')
            wellCur = groupCur[ww]

            
            #~~~~~~~~~~~~~~~~ Import mesh ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
            meshLoc = vps.dataSavePathPreProc.format(wellCur)
            meshName = vps.dataSaveNameRWh5.format(sigma, sigma, rolling_win)
            if newMeshFlag: 
                point0=Point(0,0) #What does 60 mean on next line? TODO
                point1=Point(60*len_stride,lastPixel*len_stride)
                mesh = RectangleMesh(MPI.comm_world,point0,point1,60, lastPixel)

            if oneD_flag:
                point0=Point(0,0)
                point1=Point(len_stride,lastPixel*len_stride)
                mesh = IntervalMesh(MPI.comm_world,numXPoints, 0,lastPixel*len_stride)
            else:
                mesh = Mesh()
                hdf5 = HDF5File(MPI.comm_world, meshLoc+meshName, 'r')
                hdf5.read(mesh, '/mesh', False)
            print('done importing mesh')    



            ########################### Define Functional spaces ##########################
            V=FunctionSpace(mesh, "Lagrange", 1)
            C = Function(V)
            C_n = Function(V)
            C_BC = Function(V)
            w= TestFunction(V)

            velocity_plot1 = Function(V)
            velocity_plot2 = Function(V)
            diffusion_plot = Function(V)
            reaction_plot = Function(V)

            grad_w = grad(w)
            grad_C=grad(C)

            #dxx = dx(metadata={'quadrature_degree': 1})
            dxx = dx

            V_data=FunctionSpace(mesh, "Lagrange", 1)
            C_data_n = Function(V_data)
            C_data_np1 = Function(V_data)
            grad_C_data_np1=grad(C_data_np1)
            Error_data_n = Function(V_data)


            x=SpatialCoordinate(mesh)
            y_m=lastPixel*len_stride/2
            y_l=lastPixel*len_stride

            if oneD_flag:
                y_coord=x[0]
                v1=Constant((1.0,))
            else:
                y_coord=x[1]
                v1=Constant((0.0,1.0))
                v2=Constant((1.0,0.0))

            y_tilda = (y_coord-y_m)/y_l
            v=v1
            print('done defining functional space')
            ############################## Define Boundaries ##############################
            #Shouldn't change
            if oneD_flag:
                BC1 =  CompiledSubDomain("near(x[0], side) && on_boundary", side = 0)
                BC2 =  CompiledSubDomain("near(x[0], side) && on_boundary", side = lastPixel*len_stride)
            else:
                BC1 =  CompiledSubDomain("near(x[1], side) && on_boundary", side = 0)
                BC2 =  CompiledSubDomain("near(x[1], side) && on_boundary", side = lastPixel*len_stride)
            ################# Define residue (without Boundary condition) #################
            print('done defining boundaries')

            num_basis=vps.total_num_basis
            theta=[]
            for i in range(num_basis):
                theta.append(Constant(0.0))  
                
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
            tau_ini=vps.tau_ini
            tau=Constant(tau_ini)

            W1 = VectorFunctionSpace(mesh, "CG",1)
            W2 = FunctionSpace(mesh, "DG",0)		
            lap_c = project(div(project(grad(C),W1)),W2)
            R+=tau*len_stride/4*inner(grad_w,adv_supg)*((C-C_n)/dt -diff_supg*lap_c +div(adv_supg*C)-reac_supg )*dxx

            #Jacobian for nonlinear solver
            J=derivative(R, C)
            
            #Loss function
            def loss_func(C_exp, C_sim, grad_C_exp, grad_C_sim):
                factor = 1
                loss = assemble( ((C_exp - C_sim)**2 + factor*inner(grad_C_exp - grad_C_sim, grad_C_exp - grad_C_sim))  *dx)
                return loss
            
            ####################### Run FWD solution ####################### 
            print(f'running forward solution for {plot_tag}, VSI step {stepCur}')
            
            directory = vps.fwdSavePath.format(plot_tag, wellCur, stepCur)
            if not os.path.exists(directory):
                os.makedirs(directory)
                print(f'New directory created: {directory}')
                

            
            if newMeshFlag:
                C_n.assign(Constant(1.0))
                C_data_np1.assign(C_n)
                
            else:
                hdf5.read(C_n, 'density_'+str(0))
                C_data_np1.assign(C_n)
            
            for step in range(total_num_step-1): 
                print(f'~~~~~ step = {step} ~~~~~')
                C_data_n.assign(C_data_np1)

                if newMeshFlag: 
                    #Write expressions of C_data, erk_data and akt_data in terms of coordinates x and step (time)
                    C_data_np1.assign(Constant(1.0))

                else:
                    #Experimental and boundary values
                    hdf5.read(C_data_np1,'density_'+str(step+1))

                if reinitializeFlag and step%reinitializeInterval==0: 
                    C_n.assign(C_data_n)	
                
                for i in range(num_basis):
                    theta[i].assign(np.array(gamma_matrix[step])[i])
                
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
                    prm['newton_solver']["error_on_nonconvergence"] = True

                    solver.solve()
                    C_n.assign(C)
                    if step in learn_time_steps: 
                        loss += loss_func(C_data_np1, C, grad_C_data_np1, grad_C)
            del mesh, V, C, C_n, C_BC, w, velocity_plot1, velocity_plot2, diffusion_plot, reaction_plot, grad_w, grad_C, V_data, C_data_n, C_data_np1, grad_C_data_np1, Error_data_n, x, y_m, y_l, v1 , v2, BC1, BC2, theta, basis_pool, diff_supg, adv_supg, reac_supg, dt, R, tau, W1, W2, lap_c, J 
        print('appending to loss list')
        loss_list.append([param1, param2, loss])
        print(f'loss list is {loss_list}')
        


        counter+=1
        temp = f'sensdatacurr{counter}_Group_{sigma}_{sigma}_rolling_win{rolling_win}_F{F}.dat'
        curr_sensitivity_data_filename = adjSaveDirectory + temp 
        loss_list_current = np.array(loss_list)

        with open(curr_sensitivity_data_filename, 'wb') as f:
            np.save(f, loss_list_current)
        print(f' would try to save to location: {sensitivity_data_filename}')
        
             
      print('end of param2 loop')
      
      print('just applied np.array to loss_list')
      print(f'saving loss list to locatoin: {sensitivity_data_filename}')
      print('done saving loss list')
      
    print('end of param1 loop')
    loss_list = np.array(loss_list)
    with open(sensitivity_data_filename, 'wb') as f:
      np.save(f, loss_list)
    #end of the param loops
            
print('~~~~~~~~~~~~~~~~~~~~~~~~~~~ ENDING Loss data generation for sensititvity analysis ~~~~~~~~~~~~~~~~~~~~~~~~~~~')
        
