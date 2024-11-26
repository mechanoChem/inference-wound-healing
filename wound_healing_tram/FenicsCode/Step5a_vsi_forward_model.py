from ufl import *
from dolfin import *
import numpy as np
from dolfin_adjoint import *
set_log_active(False)
import os
import sys
import pickle

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

with open('vsi_settings.pkl' , 'rb') as settingsFile:
    vps = pickle.load( settingsFile)    
    
#Note - this code only runs for 1 well at a time - need to update for one well then wrap a for loop 

################################ Problem Setup ################################

sigma = vps.sigma
rolling_win = vps.rolling_win
F = vps.best_F

total_num_step = vps.total_time
numSmallStep = vps.numSmallSteps

runtype = 2 #update to remove VSI code

oneD_flag = vps.oneD_flag
temporalParam_flag = vps.temporalParamFlag
halfDataFlag = vps.halfDataFlag
newMeshFlag = vps.newMeshFlag
reinitializeFlag = vps.reinitializeFlag
reinitializeInterval = vps.reinitializeIntervalFwd
basisString = 'A' #Need to figure out how basisString is used: ONLY used for loading basis - 
                    #Can probably replace with vps.basisSavePath and vps.basisSaveName

lastPixel = vps.lastPixel
numXPoints = vps.numXpoints                    
len_stride = vps.len_stride
nWells = vps.nWells
wells = vps.wells
wellGroups = vps.wellGroups
nGroups = vps.nGroups
nSteps = vps.total_num_basis
total_num_basis = vps.total_num_basis
nSteps = vps.endBasis - vps.startBasis
#All of the below replaced with vps
"""
sigma=3
rolling_win=3
F=200000


total_num_step=4
numSmallStep = 5

runtype = 2 #1 for VSI, 2 for Adjoint, anything else for user supplied (see near line 175 and change variable: gamma_list_userSupplied)
oneD_flag = True
temporalParam_flag=False
halfDataFlag = False
#Set newMeshFlag = True if do not want to read erk and akt data
newMeshFlag = False
reinitializeFlag = False
reinitializeInterval = 1
basisString = 'Simul'
initCondDensity = 20000

print(f'Running Adjoint fwd models')

if oneD_flag:
    if halfDataFlag: 
    	lastPixel = 10	#For half length in 1D
    	numXPoints = 37
    else:
	    lastPixel = 37  #For full length in 1D
	    numXPoints = 37
else: #Need to write mesh refinement to do high numXPoints
    if halfDataFlag: 
    	lastPixel = 10	#For half length in 1D
    	numXPoints = 10
    else:
	    lastPixel = 19  #For full length in 1D
	    numXPoints = 19


##micron_per_pix = 3.45
##mag = 10
##pix_stride = 100
##zoom_um_per_pix = micron_per_pix/mag
##zoom_um_per_pix = 0.5           #Overriding this value based on Patrick's email June 29th 
##len_stride = pix_stride * zoom_um_per_pix
len_stride = 50
"""
print('~~~~~~~~~~~~~~~~~~~~~~~~~~~ Starting VSI FWD ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~')

for gg in range(nGroups):
  groupCur = wellGroups[gg]
  for ss in range(nSteps):
    stepCur = ss+vps.startBasis
    print(f'Group: {gg+1}/{nGroups} ---- step: {ss+1}/{nSteps} ---- Opening gamma matrix')
   

    time = np.loadtxt('../data/time.dat', delimiter = ",")
    plot_tag = vps.vsiPlotTag
    #Build names for adjoint gamma location and open
#    adjSaveDirectory = vps.adjSavePath.format(f'group{gg+1}_{groupCur[0]}_to_{groupCur[-1]}_step{stepCur}')
#    adjSaveName = vps.adjSaveNameGamma.format(sigma, sigma, rolling_win, F)
    vsiSaveDirectory = vps.vsiSavePath.format(f'group{gg+1}_{groupCur[0]}_to_{groupCur[-1]}')
    vsiSaveName = vps.vsiSaveNamegammahist.format(sigma, sigma, rolling_win, F)
    gamma_hist = np.loadtxt(vsiSaveDirectory+vsiSaveName)
    gamma_values = gamma_hist[:,stepCur]
    gamma_matrix = np.tile(gamma_values, (total_num_step,1))
    print('~~~~~~~~~~~')
    print(f'gamma size = {gamma_matrix.shape}')
    print(f'gamma row 0 = {gamma_matrix[0,:]}')
    print('~~~~~~~~~~~')
    
    #Iterate through all wells in group
    for ww in range(len(groupCur)):
      print(f'Applying basis to well {groupCur[ww]} out of group {groupCur}, {ww+1}/{len(groupCur)}')
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
      
  ############################## Define Boundaries ##############################
  #Shouldn't change
      if oneD_flag:
        BC1 =  CompiledSubDomain("near(x[0], side) && on_boundary", side = 0)
        BC2 =  CompiledSubDomain("near(x[0], side) && on_boundary", side = lastPixel*len_stride)
      else:
        BC1 =  CompiledSubDomain("near(x[1], side) && on_boundary", side = 0)
        BC2 =  CompiledSubDomain("near(x[1], side) && on_boundary", side = lastPixel*len_stride)
    ################# Define residue (without Boundary condition) #################


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
      
      ####################### Run FWD solution ####################### 
      print(f'running forward solution for {plot_tag}, VSI step {stepCur}')
      directory = vps.fwdSavePath.format(plot_tag, wellCur, stepCur)
      if not os.path.exists(directory):
        os.makedirs(directory)
        print(f'New directory created: {directory}')
        
      file_plotC = XDMFFile(MPI.comm_world,directory+'density.xdmf') 
      file_plotVel = XDMFFile(MPI.comm_world,directory+'velocity.xdmf') 
      file_plotDiff = XDMFFile(MPI.comm_world,directory+'diffusion.xdmf') 
      file_plotReact = XDMFFile(MPI.comm_world,directory+'reaction.xdmf') 
      file_plotError = XDMFFile(MPI.comm_world,directory+'Error.xdmf')       
      
      
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
          prm["newton_solver"]["maximum_iterations"] = 1000
          prm['newton_solver']["error_on_nonconvergence"] = True

          solver.solve()
          C_n.assign(C)
        
        #Plotting data
        
        #Advection velocity
        if oneD_flag:
          adv_vel1_= inner(adv_supg,v1)
          adv_vel1 = project(adv_vel1_, V)
          velocity_plot1.vector()[:] = adv_vel1.vector()[:]
          file_plotVel.write(velocity_plot1, step+1)
        else: 
          adv_vel1_ = inner(adv_supg,v1) 
          adv_vel2_ = inner(adv_supg,v2)
          adv_vel1 = project(adv_vel1_, V)
          adv_vel2 = project(adv_vel2_, V)
          velocity_plot1.vector()[:] = adv_vel1.vector()[:]
          velocity_plot2.vector()[:] = adv_vel2.vector()[:]
          file_plotVel.write(velocity_plot1, step+1)
          file_plotVel.write(velocity_plot2, step+1)	
        
        #Density
        file_plotC.write(C,step+1)
        #Diffusion
        diff_ = theta[0] + theta[1]*C + theta[2]*C*C
        diff = project(diff_,V)
        diffusion_plot.vector()[:] = diff.vector()[:]
        file_plotDiff.write(diffusion_plot,step+1)
        #Reaction
        reac_ = theta[6]*C + theta[7]*C*C
        reac = project(reac_,V)
        reaction_plot.vector()[:] = reac.vector()[:]  
        file_plotReact.write(reaction_plot,step+1)
        #Error
        print(C_data_n.vector()[:])
        Error_data_n.vector()[:] = 100*np.abs( (C.vector()[:]-C_data_n.vector()[:]) /(0.0000001+C_data_n.vector()[:]))
        file_plotError.write(Error_data_n,step+1)
      file_plotVel.close() 
      file_plotC.close()
      file_plotDiff.close()
      file_plotReact.close() 
      file_plotError.close()
      print('Done with iteration, closed files')
          
        