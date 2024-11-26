# -*- coding: utf-8 -*-
"""
Created on Wed Apr 20 20:03:01 2022

@author: pkinn
"""
import numpy as np

class vsiParamStruct():
    def __init__(self):
        self.len_stride = 96*0.5 #Length repesented by one pixel
        self.oneD_flag = False #Is data 1D or 2D
        self.halfDataFlag = False #use half the data or all of it (assume symmetry)
        self.numXpoints = []
        self.lastPixel = []
        self.wells = []
        self.wells = ['%s%02d' % (r, c) for r in 'ABCD' for c in range(1, 7)]
        self.nWells = len(self.wells)
        self.time_gap = 1
        self.total_time = [] #Same as total steps
        
        self.sigma = 3 #spatial smoothing
        self.rolling_win = 5 #Temporal smoothing window
        self.wellGroups = [['A02', 'A03'], ['A04', 'A05', 'A06'], ['B01', 'B02', 'B03'], \
        ['B04', 'B05', 'B06'], ['C01', 'C02', 'C03'], ['C04', 'C05', 'C06']]
        self.wellGroups = [['A02', 'A03'], ['A04', 'A05', 'A06'], ['B01', 'B02', 'B03'], \
        ['B04', 'B05', 'B06']]
        self.wellGroups = [['A01', 'B01', 'C01', 'D01']]
        self.wellGroups = [['A01', 'B01', 'C01', 'D01'], ['A02', 'B02', 'C02', 'D02'], ['A03', 'B03', 'C03', 'D03'], \
        ['A04', 'B04', 'C04', 'D04'], ['A05', 'B05', 'C05', 'D05'], ['A06', 'B06', 'C06', 'D06']]
        self.nGroups = len(self.wellGroups)
        
        # All Paths
        self.dataLoadPathPreProc = '../data/cell_density/density/{}/'
        self.dataSavePathPreProc = '../results/PreProcess/{}/'
        self.CdataName = 'density_filter_{}_t{}.dat'
        self.dataSaveNameRWh5 = 'cell_density_{}_{}_rolling_win{}.h5'
        self.dataSaveNameRWxdmf = 'density_{}_{}_rolling_win{}.xdmf'
        self.basisSavePath = '../results/basis/Physics_Based/{}/basis_{}_{}_rolling_win{}/'
        self.basisSaveName = 'basis_step_{}.dat'
        self.vsiSavePath = '../results/VSI_gamma_matrix/{}/Physics_Based_Time_Independent/'
        self.vsiSaveNamegammahist = 'gamma_history_Group_{}_{}_rolling_win{}_F{}.dat'
        self.vsiSaveNamegamma = 'gamma_Group_{}_{}_rolling_win{}_F{}.dat'
        self.vsiSaveNameloss = 'loss_Group_{}_{}_rolling_win{}_F{}.dat'
        self.adjSavePath = '../results/Adjoint_gamma_matrix/{}/Physics_Based_Time_Independent/'
        self.adjSaveNameGamma = 'gamma_Group_{}_{}_rolling_win{}_F{}.dat'
        #Note - adj and vsi fwd solutions are separated by the plot tag (first format brackets)
        self.fwdSavePath = '../results/forward_solution/{}/well{}/step{}/' 
        
        
        #From VSI
        self.total_num_basis = 8
        self.total_steps = [] # This is just number of times
        self.temporalParamFlag = False #If false, parameters are time Independent. If true they can vary w/ time.
        self.oneD_flag = False #If true, data is 1-Dimensional. If False, it's 2-D.
        self.best_F = 200000 #F score cutoff defines when VSI stops running. If set to a very high number, VSI will proceed until it has 1 parameter
        self.vsi_Lsq_bounds = []
        self.vsi_group = [[0]] #List of lists. One parameter from each sublist will be maintained during VSI
        self.vsi_Drop_strategy = 'most_insignificant' #Or aggressive - see stepwiseRegression.py
        self.vsi_Linear_regressor = 'lsq_linear' #'lsq_linear' or 'ridge_cv' used, other options in stepwiseRegression.py
        self.vsi_ridge_cv = np.linspace(-5,1,1)
        self.vsi_sigma_n = 1.0e-14
        #From Adj
        self.numSmallSteps = []
        self.learnStepInterval = []
        self.reinitializeFlag = True
        self.startBasis = 0 #The most complex model you want to run adjoint on
        self.endBasis = self.total_num_basis #The least complex model you want to run adjoint on
        self.adjMethod = 'L-BFGS-B'
        self.adjConvFlag = True
        self.maxiter = 50        
        self.reinitializeInterval = 50
        
        #From forward solution
        self.newMeshFlag = False
        self.reinitializeIntervalFwd = 1 #Pretty sure this is different than the previous reinitializeInterval
        self.tau_ini = 100. #for supg stabiliziation in fwd model. 
        self.adjPlotTag = 'Adjoint_2D_Time_Independent' 
        self.vsiPlotTag = 'VSI_2D_Time_Independent'
        
        
