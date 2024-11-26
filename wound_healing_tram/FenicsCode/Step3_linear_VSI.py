import numpy as np
import sys
import stepwiseRegression as ST
import LeastR as LR
import os
import pickle
print('~~~~~~~~~~~~~~~~~~~~~~~~~~~ STARTING LINEAR VSI ~~~~~~~~~~~~~~~~~~~~~~~~~~~')

with open('vsi_settings.pkl' , 'rb') as settingsFile:
    vps = pickle.load( settingsFile)   
    
    
    


total_num_basis=vps.total_num_basis
total_steps=vps.total_time

temporalParam_flag=False
oneD_flag=vps.oneD_flag
best_F=vps.best_F


time=np.loadtxt('../data/time.dat',delimiter=",")
time.shape
sigma=vps.sigma
rolling_win=vps.rolling_win
wells = vps.wells
nWells = vps.nWells
wellGroups = vps.wellGroups
nGroups = vps.nGroups
print('Wells are:')
print(wells)
print(f'There are {nGroups} well groups')
print(f'Groups are {wellGroups}')
#Iterate through groups

for gg in range(nGroups):
  groupCur = wellGroups[gg]
  data=[]
  dataList = []
  #Iterate through each well in current group
  for ii in range(len(groupCur)):
      
      #dataC is an array of all data for a single well 
      dataC = []
      wellCur = groupCur[ii]
      print(f'reading data for well {wellCur} in group {gg+1}/{nGroups}')
      basisFolder = vps.basisSavePath.format(wellCur, sigma, sigma, rolling_win)
      for step in range(total_steps):
          basisFile = vps.basisSaveName.format(step)
          if step == 0:
              dataC = np.array(np.loadtxt(basisFolder + basisFile))   
          else:
              dataC = np.array(np.append(dataC, np.loadtxt(basisFolder + basisFile), 0))
          #print(f'on step {step+1}/{total_steps}, dataC shape: {dataC.shape}')
      dataList.append(dataC)
  dataAll = np.vstack(dataList)
  data = dataList[0]


  print(f'data has length {data.shape}')
  print(f'dataC has size {dataC.shape}')
  print(f'dataList[0] has size{dataList[0].shape}')
  print(f'dataAll has size{dataAll.shape} (should be 2 x dataC)')
  total_dof=int(data.shape[0]/total_steps)
  print(f'total_dof is {total_dof}')      

  total_dof=int(data.shape[0]/total_steps)

  C_dof=total_num_basis #13th basis term is only to evaluate C
  dt_list=time[1:]-time[:-1]
  used_steps=len(dt_list)
  #C_ave=np.zeros(used_steps)
  #Positive diffusivity
  """
  lsq_bounds=np.zeros((2,total_num_basis))
  lsq_bounds[0,:]=-np.inf
  lsq_bounds[1,:]=np.inf
  lsq_bounds[0,0:3]=1e-8
  #lsq_bounds[0,0:10]=1e-14
  #One element of each of the group survives
  #group=[[0,1,2,3],[4,5,6],[7,8,9,10,11,12]]
  group=[[0]] #Want some diffusion for stability, other terms are optional 
  drop_strategy='most_insignificant'
  linear_regressor='lsq_linear'
  #linear_regressor='ridge_cv'
  ridge_cv=np.logspace(-5,1, 1)      
  sigma_n=1.0e-12
  """
  lsq_bounds = vps.vsi_Lsq_bounds
  group = vps.vsi_group
  drop_strategy = vps.vsi_Drop_strategy
  linear_regressor = vps.vsi_Linear_regressor
  ridge_cv = vps.vsi_ridge_cv
  sigma_n = vps.vsi_sigma_n  
  
  vsiSaveDirectory = vps.vsiSavePath.format(f'group{gg+1}_{groupCur[0]}_to_{groupCur[-1]}')
  print(f'VsiSaveDirectory is: {vsiSaveDirectory}')
  
  if temporalParam_flag:
      for i in range(used_steps):
          y=(data[total_dof*(i+1):total_dof*(i+2),C_dof]-data[total_dof*i:total_dof*(i+1),C_dof])/dt_list[i]
          theta_matrix=data[total_dof*(i+1):total_dof*(i+2),0:total_num_basis]
          #active points - removes end points for sure and maybe more
          active_index=np.reshape(np.where(np.abs(y)>1.0e-15),(-1))
          y=y[active_index]
          
          theta_matrix=theta_matrix[active_index,:]
          
      #C_ave[i]=np.sum(theta_matrix[:,C_dof])
          #continue
      
    
          #y.flatten()
          print('Step num = %i / %i'%(i, used_steps))
          #exit()
          if linear_regressor=='lsq_linear':
              model = ST.stepwiseR(F_criteria=[best_F],linear_regressor=linear_regressor,group=group,lsq_bounds=lsq_bounds,basis_drop_strategy=drop_strategy,sigma_n=sigma_n)
          elif linear_regressor=='ridge_cv':
              model = ST.stepwiseR(F_criteria=[best_F],linear_regressor=linear_regressor,group=group,basis_drop_strategy=drop_strategy, ridge_cv= ridge_cv,sigma_n=sigma_n)

      
          model.stepwiseR_fit(theta_matrix,y)
      
      
          if i==0:
              gamma_matrix_all=np.reshape(model.gamma_matrix[:,-1],(1,-1))
              loss_all=np.array([model.loss[-1]])
          else:
              gamma_matrix_all=np.append(gamma_matrix_all,np.reshape(model.gamma_matrix[:,-1],(1,-1)),0)
              loss_all=np.append(loss_all,model.loss[-1])
      
  else:
      print('Running time independent VSI')
      yAll = []
      thetaMatrixAll = []
      for ii in range(len(groupCur)):
          print(f'ii = {ii}')
          dataC = dataList[ii]
          yc = np.zeros((total_dof*used_steps, 1))
          gamma_time_extend_map = np.zeros((total_num_basis*used_steps,total_num_basis))
          print(f'yc size is {yc.shape}')
          for i in range(used_steps):
              print(f'for i={i} in range(used_steps): {used_steps}')
              #yc[total_dof*(i):total_dof*(i+1)] =(data[total_dof*(i+1):total_dof*(i+2),C_dof]-data[total_dof*i:total_dof*(i+1),C_dof])/dt_list[i]
              #gamma_time_extend_map[total_num_basis*(i):total_num_basis*(i+1),:] = np.eye(total_num_basis)
              yc[total_dof*(i):total_dof*(i+1),0]=(dataC[total_dof*(i+1):total_dof*(i+2),C_dof]-dataC[total_dof*i:total_dof*(i+1),C_dof])/dt_list[i]
          theta_matrix=dataC[0:total_dof*used_steps,0:total_num_basis]
          yAll.append(yc)
          thetaMatrixAll.append(theta_matrix)          
      y = np.vstack(yAll)
      y = np.squeeze(y)
      theta_matrix = np.vstack(thetaMatrixAll)
      
      print('~~~~~~~~~~~~~~~~')
      """
      print(f'dataAll should be bigger than theta_matrix by {len(groupCur)} * nDataPerTime')
      print(f'dataAll size = {dataAll.shape}')
      print(f'theta_matrix size = {theta_matrix.shape}')
      print(f'yc size = {yc.shape}')
      print(f'y size = {y.shape}')
      print(f'starting fit')
      print(f'values of y are {y}')
      """
    
      if linear_regressor=='lsq_linear':
          model = ST.stepwiseR(F_criteria=[best_F],linear_regressor=linear_regressor,group=group,lsq_bounds=lsq_bounds,basis_drop_strategy=drop_strategy,sigma_n=sigma_n)
      elif linear_regressor=='ridge_cv':
          model = ST.stepwiseR(F_criteria=[best_F],linear_regressor=linear_regressor,group=group,basis_drop_strategy=drop_strategy, ridge_cv= ridge_cv,sigma_n=sigma_n)

      
      model.stepwiseR_fit(theta_matrix,y)
      
      print(model.gamma_matrix)
      gamma_matrix_all = np.ones((used_steps,1))*(model.gamma_matrix[:,-1].flatten())
      loss_all=np.array([model.loss])
      """
      print(loss_all)
      print('###########################################')
      print(model.gamma_matrix_history)
      print('###########################################')
      print('###########################################')
      print(model.loss_history)
      print('###########################################')
      print(len(model.gamma_matrix_history))
      for ii in range(len(model.gamma_matrix_history)):
          print(model.gamma_matrix_history[ii].shape)
      print('#####################')
      """
      print(f'For group {gg+1}, gamma_hist is:')
      print(model.gamma_matrix_history[ii])

  if not os.path.exists(vsiSaveDirectory):
    os.makedirs(vsiSaveDirectory)
    print(f'New directory created: {vsiSaveDirectory}')
  else:
    print(f'Saving to: {vsiSaveDirectory}')
  
  np.savetxt(vsiSaveDirectory + vps.vsiSaveNamegammahist.format(sigma, sigma, rolling_win, best_F), model.gamma_matrix_history[0])
  np.savetxt(vsiSaveDirectory + vps.vsiSaveNamegamma.format(sigma, sigma, rolling_win, best_F), gamma_matrix_all)
  np.savetxt(vsiSaveDirectory + vps.vsiSaveNameloss.format(sigma, sigma, rolling_win, best_F), loss_all)
  print('~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~')
print('~~~~~~~~~~~~~~~~~~~~~~~~~~~ ENDING LINEAR VSI ~~~~~~~~~~~~~~~~~~~~~~~~~~~')

"""
  if oneD_flag:
      np.savetxt(vsiSaveDirectory + '/gamma_Group_'+str(sigma)+'_'+str(sigma)+'_rolling_win'+str(rolling_win)+'_F'+str(best_F)+'_refine4'+'.dat',gamma_matrix_all)
      np.savetxt(vsiSaveDirectory + '/loss_Group_'+str(sigma)+'_'+str(sigma)+'_rolling_win'+str(rolling_win)+'_F'+str(best_F)+'_refine4'+'.dat',loss_all)
      np.savetxt(vsiSaveDirectory + '/gamma_history_Group_'+str(sigma)+'_'+str(sigma)+'_rolling_win'+str(rolling_win)+'_F'+str(best_F)+'_refine4'+'.dat',model.gamma_matrix_history[0])

  else:
      np.savetxt(vsiSaveDirectory + '/gamma_Group_'+str(sigma)+'_'+str(sigma)+'_rolling_win'+str(rolling_win)+'_F'+str(best_F)+'.dat',gamma_matrix_all)
      np.savetxt(vsiSaveDirectory + '/loss_Group_'+str(sigma)+'_'+str(sigma)+'_rolling_win'+str(rolling_win)+'_F'+str(best_F)+'.dat',loss_all)	
      np.savetxt(vsiSaveDirectory + '/gamma_history_Group_'+str(sigma)+'_'+str(sigma)+'_rolling_win'+str(rolling_win)+'_F'+str(best_F)+'.dat',model.gamma_matrix_history[0])
"""
