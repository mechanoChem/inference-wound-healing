from ufl import *
from dolfin import *
import numpy as np
from dolfin_adjoint import *
set_log_active(False)
import os
import sys
import pickle
import matplotlib.pyplot as plt


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
print('~~~~~~~~~~~~~~~~~~~~~~~~~~~ STARTING ADJ FWD MODEL ~~~~~~~~~~~~~~~~~~~~~~~~~~~')

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
nSteps = vps.endBasis-vps.startBasis

#For loss estimation
learnStepInterval = vps.learnStepInterval
learn_time_steps = np.arange(1,total_num_step,learnStepInterval)
#All of the below replaced with vps

for gg in range(nGroups): 
  groupCur = wellGroups[gg]
  for ss in [nSteps-1]:

    print(f'Group: {gg+1}/{nGroups} ---- step: {ss+1}/{nSteps} ---- Opening gamma matrix')
    stepCur = ss+vps.startBasis

    time = np.loadtxt('../data/time.dat', delimiter = ",")
    plot_tag = vps.adjPlotTag
    #Build names for adjoint gamma location and open
    adjSaveDirectory = vps.adjSavePath.format(f'group{gg+1}_{groupCur[0]}_to_{groupCur[-1]}_step{stepCur}')
    adjSaveName = vps.adjSaveNameGamma.format(sigma, sigma, rolling_win, F)
    gamma_matrix_original = np.loadtxt(adjSaveDirectory+adjSaveName)
    temp = 'sensdata_Group_{sigma}_{sigma}_rolling_win{rolling_win}_F{F}.dat'
    sensitivity_data_filename = adjSaveDirectory+temp
    print('~~~~~~~~~~~')
    print(f'gamma size = {gamma_matrix_original.shape}')
    print(f'gamma row 0 = {gamma_matrix_original[0,:]}')
    print('~~~~~~~~~~~')

    #SS: mean param values
    #D0_mean = 6.8
    #C1_mean = 0.05
    #C2_mean = -34.0
    
    #Need to be equally placed for contour function
    D0_range = [1,6,11,16,21]
    C1_range = [0.01,0.03,0.05,0.07,0.09]
    C2_range = [10.,20.,30.,40.,50.]

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
      xlabel = 'D0'
      ylabel = 'C1'
    elif run_combo==1:
      param_pos = [0,7]
      param_vals1 = D0_range
      param_vals2 = C2_range      
      xlabel = 'D0'
      ylabel = 'C2'
    elif run_combo==2:
      param_pos = [6,7]
      param_vals1 = C1_range
      param_vals2 = C2_range
      xlabel = 'C1'
      ylabel = 'C2'
    else:  
      exit()
    NumP1 = len(param_vals1)
    NumP2 = len(param_vals2)

    loss_list = []
    #Temp loss calculation - replace with reading data
    #for param1 in param_vals1: 
    #  for param2 in param_vals2:
    #    loss_list.append([param1, param2, (param1/11-1)**2+400*(param2/0.05-1)**2])
    with open(sensitivity_data_filename, 'rb') as f:
      loss_list = np.load(f)   

    #Plotting code goes here
    X = loss_list[:,0].reshape((NumP1,NumP2))
    Y = loss_list[:,1].reshape((NumP1,NumP2))
    Z = loss_list[:,2].reshape((NumP1,NumP2))

    fig,ax=plt.subplots(1,1)
    cp = ax.contourf(X, Y, Z)
    fig.colorbar(cp) # Add a colorbar to a plot
    ax.set_title('Filled Contours Plot')
    ax.set_xlabel(xlabel)
    ax.set_ylabel(ylabel)
    plt.savefig('test.png')


