import numpy as np
import sys
import stepwiseRegression as ST
import LeastR as LR
import os

total_num_basis=8
total_steps=5

temporalParam_flag=False
oneD_flag=True
best_F=200000


time=np.loadtxt('../data/time.dat',delimiter=",")
time.shape
sigma=3
rolling_win=1
total_densities = 6
dense_vec = [20000, 18000, 16000, 14000, 12000, 10000]
for ii in range(total_densities):
    basisFolder = '../results/basis/Physics_Based/density' + str(dense_vec[ii]) + '/'
    data=[]
    for step in range(total_steps):
        if step==0:
            if oneD_flag:
                data=np.array(np.loadtxt(basisFolder + 'basis_1D_'+str(sigma)+'_'+str(sigma)+'_rolling_win'+str(rolling_win)+'/basis_step_'+str(step)+'_refine4'+'.dat'))
            else: 
                data=np.array(np.loadtxt(basisFolder + 'basis_2D_'+str(sigma)+'_'+str(sigma)+'_rolling_win'+str(rolling_win)+'/basis_step_'+str(step)+'.dat'))
        else:
            if oneD_flag:
                data=np.array(np.append(data,np.loadtxt(basisFolder + 'basis_1D_'+str(sigma)+'_'+str(sigma)+'_rolling_win'+str(rolling_win)+'/basis_step_'+str(step)+'_refine4'+'.dat'),0))
            else:
                data=np.array(np.append(data,np.loadtxt(basisFolder + 'basis_2D_'+str(sigma)+'_'+str(sigma)+'_rolling_win'+str(rolling_win)+'/basis_step_'+str(step)+'.dat'),0))		


    total_dof=int(data.shape[0]/total_steps)

    C_dof=8 #13th basis term is only to evaluate C
    dt_list=time[1:]-time[:-1]
    used_steps=len(dt_list)
    #C_ave=np.zeros(used_steps)
    #Positive diffusivity
    lsq_bounds=np.zeros((2,total_num_basis))
    lsq_bounds[0,:]=-np.inf
    lsq_bounds[1,:]=np.inf
    lsq_bounds[0,0:10]=1e-8
    #One element of each of the group survives
    #group=[[0,1,2,3],[4,5,6],[7,8,9,10,11,12]]
    group=[[0],[1,2]]; #Want some diffusion for stability, other terms are optional 
    drop_strategy='most_insignificant'
    linear_regressor='lsq_linear'
    #linear_regressor='ridge_cv'
    ridge_cv=np.logspace(-5,1, 1)
        
    sigma_n=1.0e-14

    if temporalParam_flag:
        if oneD_flag:
             directory = '../results/VSI_gamma_matrix/Physics_Based_Time_Dependent_1D/density' + str(dense_vec[ii])
        else:
             directory = '../results/VSI_gamma_matrix/Physics_Based_Time_Dependent_2D/density' + str(dense_vec[ii])
                     
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
        if oneD_flag:
             directory = '../results/VSI_gamma_matrix/Physics_Based_Time_Independent_1D/density' + str(dense_vec[ii])
        else:
             directory = '../results/VSI_gamma_matrix/Physics_Based_Time_Independent_2D/density' + str(dense_vec[ii])	
                     
        y = np.zeros((total_dof*used_steps))
        gamma_time_extend_map = np.zeros((total_num_basis*used_steps,total_num_basis))
        for i in range(used_steps):
            y[total_dof*(i):total_dof*(i+1)] =(data[total_dof*(i+1):total_dof*(i+2),C_dof]-data[total_dof*i:total_dof*(i+1),C_dof])/dt_list[i]
            gamma_time_extend_map[total_num_basis*(i):total_num_basis*(i+1),:] = np.eye(total_num_basis)
        theta_matrix=data[0:total_dof*used_steps,0:total_num_basis]
            

            
        if linear_regressor=='lsq_linear':
            model = ST.stepwiseR(F_criteria=[best_F],linear_regressor=linear_regressor,group=group,lsq_bounds=lsq_bounds,basis_drop_strategy=drop_strategy,sigma_n=sigma_n)
        elif linear_regressor=='ridge_cv':
            model = ST.stepwiseR(F_criteria=[best_F],linear_regressor=linear_regressor,group=group,basis_drop_strategy=drop_strategy, ridge_cv= ridge_cv,sigma_n=sigma_n)

        
        model.stepwiseR_fit(theta_matrix,y)
            
        
        #gamma_matrix_all=np.reshape(model.gamma_matrix[:,-1],(1,-1))
        #loss_all=np.array([model.loss[-1]])

        print(model.gamma_matrix)
        gamma_matrix_all = np.ones((used_steps,1))*(model.gamma_matrix[:,-1].flatten())
        loss_all=np.array([model.loss])
        #print(loss_all)
        print('###########################################')
        #print(model.gamma_matrix_history)
        print('###########################################')
        print('###########################################')
        #print(model.loss_history)
        print('###########################################')
        #print(len(model.gamma_matrix_history))
        for ii in range(len(model.gamma_matrix_history)):
            #print(model.gamma_matrix_history[ii].shape)
            print('#####################')
            #print(model.gamma_matrix_history[ii])

    if not os.path.exists(directory):
      os.makedirs(directory)
      print('New directory created: %s'%directory)

    if oneD_flag:
        np.savetxt(directory + '/gamma_Group_'+str(sigma)+'_'+str(sigma)+'_rolling_win'+str(rolling_win)+'_F'+str(best_F)+'_refine4'+'.dat',gamma_matrix_all)
        np.savetxt(directory + '/loss_Group_'+str(sigma)+'_'+str(sigma)+'_rolling_win'+str(rolling_win)+'_F'+str(best_F)+'_refine4'+'.dat',loss_all)
        np.savetxt(directory + '/gamma_history_Group_'+str(sigma)+'_'+str(sigma)+'_rolling_win'+str(rolling_win)+'_F'+str(best_F)+'_refine4'+'.dat',model.gamma_matrix_history[0])
    else:
        np.savetxt(directory + '/gamma_Group_'+str(sigma)+'_'+str(sigma)+'_rolling_win'+str(rolling_win)+'_F'+str(best_F)+'.dat',gamma_matrix_all)
        np.savetxt(directory + '/loss_Group_'+str(sigma)+'_'+str(sigma)+'_rolling_win'+str(rolling_win)+'_F'+str(best_F)+'.dat',loss_all)	
