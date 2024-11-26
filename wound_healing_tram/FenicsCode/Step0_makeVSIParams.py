import numpy as np
from vsiParams import vsiParamStruct as vps
import pickle

#Instantiate vsi params
v = vps()
#general params
v.len_stride = 96*0.5 #microns
v.numXPoints = 79
v.lastPixel = 41
#v.total_densities = 6
#v.dense_vec = ['20000', '18000', '16000', '14000', '12000', '10000']
v.time_gap = 1
v.total_time = 109
v.sigma = 3
v.rolling_win = 5
#Paths and names for files
"""
v.dataLoadPathPreProc = '../data/cell_density/density/{}/'
v.dataSavePathPreProc = '../results/PreProcess/well{}/'
v.dataSaveNameRWh5 = 'cell_density_1D_{}_{}_rolling_win{}_refine4.h5'
v.dataSaveNameRWxdmf = 'density_1D_{}_{}_rolling_win{}_refine4.xdmf'
"""
#TODO: Determine what should be set here and what should be set in default struct



#VSI settings
v.total_num_basis = 8
v.total_steps = v.total_time
v.temporalParamFlag = False
v.oneD_flag = False
v.bestF = 200000
v.vsi_Lsq_bounds=np.zeros((2,v.total_num_basis))
v.vsi_Lsq_bounds[0,:]=-np.inf
v.vsi_Lsq_bounds[1,:]=np.inf
v.vsi_Lsq_bounds[0,0:3]=1e-8
v.vsi_group = [[0]]
v.vsi_Drop_strategy = 'most_insignificant'
v.vsi_Linear_regressor = 'lsq_linear'
v.vsi_ridge_cv = np.linspace(-5,1,1)
v.vsi_sigma_n = 1.0e-14

#adj settings
v.numSmallSteps = 2
v.learnStepInterval = 9 #Only the starting point
v.startBasis = 5
v.endBasis = v.total_num_basis
v.reinitializeInterval = 10


#SS: Reinit test
v.maxiter = 10
v.reinitializeFlag = True
#v.reinitializeInterval = 20


with open('vsi_settings.pkl' , 'wb') as settingsFile:
    pickle.dump(v, settingsFile)    
