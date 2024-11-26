"""
Zhenlin Wang 2019
"""

import numpy as np
import LeastR as LR 
class stepwiseR(object):
  def __init__(self,F_criteria=[1.0e10],F_switch=[],basis_drop_strategy='aggressive',linear_regressor='LinearRegression',sigma_n=1.0e-16,anchor_index=[-1],group=[-1],alpha_lasso=0,alpha_ridge=0, ridge_cv=[0],lsq_bounds=(0),threshold_d=1.0e-14,n_jobs=1):
    self.F_criteria=F_criteria
    self.F_switch=F_switch
    self.sigma_n=sigma_n
    self.anchor_index=anchor_index
    self.alpha_lasso=alpha_lasso
    self.alpha_ridge=alpha_ridge
    self.ridge_cv=ridge_cv
    self.lsq_bounds=lsq_bounds
    self.n_jobs=n_jobs
    self.threshold_d=threshold_d
    self.basis_drop_strategy=basis_drop_strategy
    self.linear_regressor=linear_regressor
    self.last_F=0
    self.group=group

    
  def test(self):
    print('stepwiseR test_pass')
    
  def stepwiseR_fit_aggressive(self, theta_matrix, X_matrix):
    _,n_base_orign=theta_matrix.shape
    self.anchor=np.zeros(n_base_orign)
    if self.anchor_index[0]!=-1:
      for key in self.anchor_index:
        self.anchor[key]=1

    self.loss=np.zeros(n_base_orign)
    self.score=np.zeros(n_base_orign)
    self.F_index=np.zeros(n_base_orign)
    self.gamma_matrix=np.zeros((n_base_orign,n_base_orign))
    alpha_sum=self.alpha_lasso+self.alpha_ridge+self.ridge_cv[0]
    threshold_d=self.threshold_d
    self.best_alpha=np.zeros(n_base_orign)
    
    # local_to_global_index
    local_to_global_index=np.arange(n_base_orign)
    F_threshold=self.F_criteria[0]
    
    #########
    #first LS_regression
    #########
    num_column=0
    if self.linear_regressor=='LinearRegression':
        [gamma_vector,self.loss[0]]=LR.fit(theta_matrix,X_matrix)
    if self.linear_regressor=='lasso':
        [gamma_vector,self.loss[0]]=LR.fit_lasso(theta_matrix,X_matrix, alpha=self.alpha_lasso)
    if self.linear_regressor=='ridge':
        [gamma_vector,self.loss[0],self.score[0]]=LR.fit_ridge(theta_matrix,X_matrix, alpha=self.alpha_ridge)
    if self.linear_regressor=='ridge_cv':
        [gamma_vector,self.loss[0],self.score[0],self.best_alpha[0]]=LR.fit_ridge_cv(theta_matrix,X_matrix, alpha=self.ridge_cv)
    if self.linear_regressor=='lsq_linear':
        [gamma_vector,self.loss[0]]=LR.fit_lsq_linear(theta_matrix,X_matrix, bounds=self.lsq_bounds)
        
    self.gamma_matrix[local_to_global_index,num_column]=gamma_vector
    
    #########
    #stepwise
    #########
    num_column=num_column+1;
    num_canditate_basis=n_base_orign;
    frozen_index=[]
    stepwiseStep = 0
    while num_canditate_basis>1:
      stepwiseStep+=1
      print('stepwiseStep = ')
      print(stepwiseStep)
 
      #get current F_criteria
      for i in range(len(self.F_switch)):
        if num_column>self.F_switch[i]:
          F_threshold=self.F_criteria[i+1]
        else:
          break
          
      find_flag=False 
      # put anchor index into frozen_index
      for i in range(local_to_global_index.size) :
        if self.anchor[local_to_global_index[i]]==1 :
          frozen_index.append(i)
          
      # begin to do basis reduction
      for j in range(gamma_vector.size):
        # continue if j is in the frozen_index
        if j in frozen_index:
          continue
        # calculate the min of gamma_vector except the frozen_index
        gamma_vector_min=gamma_vector;
        gamma_vector_min=np.delete(gamma_vector_min, frozen_index)
        gamma_criteria=min(abs(gamma_vector_min) )+threshold_d;
        theta_matrix_try=theta_matrix;
        
        #tentative delete the basis
        if abs(gamma_vector[j])<gamma_criteria :
          frozen_index.append(j)
          find_flag=True
          # delete the corresponding column
          theta_matrix_try=np.delete(theta_matrix_try,j,1)
            
            
          if self.linear_regressor=='LinearRegression':
            [gamma_vector_try,loss_try]=LR.fit(theta_matrix_try,X_matrix)
          if self.linear_regressor=='lasso':
            [gamma_vector_try,loss_try]=LR.fit_lasso(theta_matrix_try,X_matrix, alpha=self.alpha_lasso)
          if self.linear_regressor=='ridge':
            [gamma_vector_try,loss_try,score_tem]=LR.fit_ridge(theta_matrix_try,X_matrix, alpha=self.alpha_ridge)
          if self.linear_regressor=='ridge_cv':
              [gamma_vector_try,loss_try,score_tem,best_alpha_tem]=LR.fit_ridge_cv(theta_matrix_try,X_matrix, alpha=self.ridge_cv)
          if self.linear_regressor=='lsq_linear':
              [gamma_vector_try,loss_try]=LR.fit_lsq_linear(theta_matrix_try,X_matrix, bounds=np.delete(self.lsq_bounds,j,1) )
                
          F=(loss_try-self.loss[num_column-1])/self.loss[num_column-1]*(n_base_orign-local_to_global_index.size+1)   
          if(F>self.last_F):
            self.last_F=F
                   
          # do F_test
          if F<F_threshold or loss_try<self.sigma_n:
            theta_matrix=np.delete(theta_matrix,j,1)
            local_to_global_index=np.delete(local_to_global_index,j)
            if self.linear_regressor=='lsq_linear':
              self.lsq_bounds=np.delete(self.lsq_bounds,drop_index,1) 
            self.F_index[num_column]=F
            if(len(self.ridge_cv)>2):
                self.best_alpha[num_column]=best_alpha_tem
            self.loss[num_column]=loss_try
            gamma_vector=gamma_vector_try
            self.gamma_matrix[local_to_global_index,num_column]=gamma_vector
            num_column=num_column+1
            num_canditate_basis=num_canditate_basis-1
            frozen_index=[]
            
        # break tentative deleting basis
        if find_flag==True:
          break
      #stop the algorithm   
      if find_flag==0 or gamma_vector_min.size<1:
        break
        
    self.gamma_matrix=np.delete(self.gamma_matrix,np.arange(num_column,n_base_orign), axis=1)
    self.loss=np.delete(self.loss,np.arange(num_column,n_base_orign))   
    self.F_index=np.delete(self.F_index,np.arange(num_column,n_base_orign))   
    self.best_alpha=np.delete(self.best_alpha,np.arange(num_column,n_base_orign))  
                    
  
  def stepwiseR_fit_most_insignificant(self, theta_matrix, X_matrix):
    _,n_base_orign=theta_matrix.shape
    self.anchor=np.zeros(n_base_orign)
    if self.anchor_index[0]!=-1:
      for key in self.anchor_index:
        self.anchor[key]=1

    self.loss=np.zeros(n_base_orign)
    self.score=np.zeros(n_base_orign)
    self.F_index=np.zeros(n_base_orign)
    self.gamma_matrix=np.zeros((n_base_orign,n_base_orign))
    alpha_sum=self.alpha_lasso+self.alpha_ridge+self.ridge_cv[0]
    threshold_d=self.threshold_d
    self.best_alpha=np.zeros(n_base_orign)
    
    # local_to_global_index
    local_to_global_index=np.arange(n_base_orign)
    F_threshold=self.F_criteria[0]
    
    #########
    #first LS_regression
    #########
    num_column=0
    if self.linear_regressor=='LinearRegression':
        [gamma_vector,self.loss[0]]=LR.fit(theta_matrix,X_matrix)
    if self.linear_regressor=='lasso':
        [gamma_vector,self.loss[0]]=LR.fit_lasso(theta_matrix,X_matrix, alpha=self.alpha_lasso)
    if self.linear_regressor=='ridge':
        [gamma_vector,self.loss[0],self.score[0]]=LR.fit_ridge(theta_matrix,X_matrix, alpha=self.alpha_ridge)
    if self.linear_regressor=='ridge_cv':
        [gamma_vector,self.loss[0],self.score[0],self.best_alpha[0]]=LR.fit_ridge_cv(theta_matrix,X_matrix, alpha=self.ridge_cv)
    if self.linear_regressor=='lsq_linear':
        [gamma_vector,self.loss[0]]=LR.fit_lsq_linear(theta_matrix,X_matrix, bounds=self.lsq_bounds)
        
    self.gamma_matrix[local_to_global_index,num_column]=gamma_vector
    self.gamma_matrix_history = []
    self.loss_history = []
    self.gamma_matrix_history.append(self.gamma_matrix)
    self.loss_history.append(self.loss)
    #########
    #stepwise
    #########
    num_column=num_column+1;
    num_canditate_basis=n_base_orign;
    frozen_index=[]
    stepwiseStep = 0
    while num_canditate_basis>1:
      stepwiseStep+=1
      print('stepwiseStep = ')
      print(stepwiseStep)
      #get current F_criteria
      for i in range(len(self.F_switch)):
        if num_column>self.F_switch[i]:
          F_threshold=self.F_criteria[i+1]
        else:
          break
          
      find_flag=False 
      
      if self.group[0]!=-1:
        for group_index in range(len(self.group)):
          if len(self.group[group_index])==1:
            self.anchor[self.group[group_index][0]]=1
            
      # put anchor index into frozen_index
      for i in range(local_to_global_index.size) :
        if self.anchor[local_to_global_index[i]]==1 :
          frozen_index.append(i)
          
      # begin to do basis reduction
      loss_tem=np.ones(gamma_vector.size)*1.0e10
      best_alpha_tem=np.zeros(gamma_vector.size)
      score_tem=np.zeros(gamma_vector.size)
      gamma_matrix_try=np.zeros((gamma_vector.size-1,gamma_vector.size))
      for j in range(gamma_vector.size):
        # continue if j is in the frozen_index
        if j in frozen_index:
          continue
        
        theta_matrix_try=np.delete(theta_matrix,j,1)
        if self.linear_regressor=='LinearRegression':
          [gamma_matrix_try[:,j],loss_tem[j] ]=LR.fit(theta_matrix_try,X_matrix)
        if self.linear_regressor=='lasso':
          [gamma_vector_try[:,j],loss_tem[j] ]=LR.fit_lasso(theta_matrix_try,X_matrix, alpha=self.alpha_lasso)
        if self.linear_regressor=='ridge':
          [gamma_matrix_try[:,j],loss_tem[j] ,score_tem[j] ]=LR.fit_ridge(theta_matrix_try,X_matrix, alpha=self.alpha_ridge)
        if self.linear_regressor=='ridge_cv':
          [gamma_matrix_try[:,j],loss_tem[j] ,score_tem[j],best_alpha_tem[j] ]=LR.fit_ridge_cv(theta_matrix_try,X_matrix, alpha=self.ridge_cv)
        if self.linear_regressor=='lsq_linear':
          [gamma_matrix_try[:,j],loss_tem[j]]=LR.fit_lsq_linear(theta_matrix_try,X_matrix, bounds=np.delete(self.lsq_bounds,j,1) )
            
      drop_index=np.argmin(loss_tem)  
      loss_try=loss_tem[drop_index]  
      F=(loss_try-self.loss[num_column-1])/self.loss[num_column-1]*(n_base_orign-local_to_global_index.size+1) 
      ###
      print(f'New F = {F} last F = {self.last_F} Fcutoff = {F_threshold} Loss_try = {loss_try}')
      ###
      if(F>self.last_F):
        self.last_F=F
      # do F_test
      if F<F_threshold or loss_try<self.sigma_n:
        find_flag=True
        theta_matrix=np.delete(theta_matrix,drop_index,1)
        if self.group[0]!=-1:
          for group_index in range(len(self.group)):
            if local_to_global_index[drop_index] in self.group[group_index]:
              self.group[group_index].remove(local_to_global_index[drop_index])
              break
              
        local_to_global_index=np.delete(local_to_global_index,drop_index)
        if self.linear_regressor=='lsq_linear':
          self.lsq_bounds=np.delete(self.lsq_bounds,drop_index,1) 
        self.F_index[num_column]=F
        if(len(self.ridge_cv)>2):
          self.best_alpha[num_column]=best_alpha_tem[drop_index]
        self.loss[num_column]=loss_try
        gamma_vector=gamma_matrix_try[:,drop_index]
        self.gamma_matrix[local_to_global_index,num_column]=gamma_vector
        num_column=num_column+1
        num_canditate_basis=num_canditate_basis-1
        frozen_index=[]
      
      self.gamma_matrix_history.append(self.gamma_matrix)
      self.loss_history.append(self.loss)
      #stop the algorithm of no operator can be eliminated or only one operator left.
      if find_flag==False:
        break
        
    self.gamma_matrix=np.delete(self.gamma_matrix,np.arange(num_column,n_base_orign), axis=1)
    self.loss=np.delete(self.loss,np.arange(num_column,n_base_orign))   
    self.F_index=np.delete(self.F_index,np.arange(num_column,n_base_orign))   
    self.best_alpha=np.delete(self.best_alpha,np.arange(num_column,n_base_orign))  
      
    
  ###################################################            
  def stepwiseR_fit(self, theta_matrix, X_matrix):
    if self.basis_drop_strategy=='aggressive':
      self.stepwiseR_fit_aggressive(theta_matrix, X_matrix)
    elif self.basis_drop_strategy=='most_insignificant':
      self.stepwiseR_fit_most_insignificant(theta_matrix, X_matrix)


    

    
    
      
    
    