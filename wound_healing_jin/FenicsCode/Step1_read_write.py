from ufl import *
from dolfin import *
import numpy as np
import csv
import os

def central_rolling_mean(x,window):
    assert (window%2 == 1), "Rolling average window must be odd"
    N = int((window - 1)/2)
    total_time=len(x)
    shape=x[0].shape
    xx=[0]*total_time
    for i in range(total_time):
        left=i-N
        if left<0:
            left=0
        right=i+N
        if right>total_time-1:
            right=total_time-1
        tem=np.zeros(shape)
        for j in np.arange(left,right+1):
            tem += x[j]/(right-left+1)
        xx[i] = tem
    
    return xx

# Create mesh and define function space
#0.4375 x 0.5 x 0.5
micron_per_pix = 3.45
mag = 10
pix_stride = 100
zoom_um_per_pix = micron_per_pix/mag
zoom_um_per_pix = 0.5           #Overriding this value based on Patrick's email June 29th 
len_stride = pix_stride * zoom_um_per_pix
#each bin is 50 um in other data
len_stride = 0.05 #in millimeters PCK 12/8
len_stride = 50 # in microns

oneD_flag=True
halfDataFlag = False
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


total_densities = 6
dense_vec = [20000, 18000, 16000, 14000, 12000, 10000]

for ii in range(total_densities):
    dataLoadPath = '../data/cell_density/initCells' + str(dense_vec[ii]) + '/'
    dataSavePath =  '../results/PreProcess/density'+str(dense_vec[ii])
    if not os.path.isdir(dataSavePath):
        print('making directory: ' + dataSavePath)
        os.makedirs(dataSavePath)
    else:
        print('directory: ' + dataSavePath + ' already exists.')
    point0=Point(0,0)
    #TODO why is next line 60?
    numMeshPts = 37


    if oneD_flag:
      point0=Point(0,0)
      point1=Point(len_stride,lastPixel*len_stride)
      mesh = IntervalMesh(MPI.comm_world,numXPoints, 0,lastPixel*len_stride)
      #mesh = RectangleMesh(MPI.comm_world,point0,point1,1, 19)
      

    V = FunctionSpace(mesh, "Lagrange", 1)
    TS= TensorFunctionSpace(mesh, "CG", 1)
    C = Function(V) 
    erk = Function(V) 
    akt = Function(V) 

    d2v = dof_to_vertex_map(V)

    time_gap=1

    total_time=5

    C_data_all=[0]*total_time
    sigma=3
    if oneD_flag:
        for t in range(total_time):
            C_data_all[t]=np.loadtxt(dataLoadPath + 'density_1D_filter_'+str(sigma)+'_t'+str(t+1)+'_refine4'+'.dat',delimiter=",")[:numXPoints+1]
           
    else:
        for t in range(total_time):
            C_data_all[t]=np.loadtxt('../data/cell_density/density/density_filter_'+str(sigma)+'_t'+str(t+1)+'.dat',delimiter=",")
            #To read only half the data


    rolling_win=1
    C_data_all=central_rolling_mean(C_data_all,rolling_win) 


    if oneD_flag:

        file = HDF5File(MPI.comm_world, dataSavePath+ '/cell_density_1D_'+str(sigma)+'_'+str(sigma)+'_rolling_win'+str(rolling_win)+'_refine4'+'.h5', 'w')
        file.write(mesh,'/mesh')
        file_C = XDMFFile(MPI.comm_world,dataSavePath+ '/density_1D_'+str(sigma)+'_'+str(sigma)+'_rolling_win'+str(rolling_win)+'_refine4'+'.xdmf');
        
    else:
        file = HDF5File(MPI.comm_world, '../results/PreProcess/cell_density_'+str(sigma)+'_'+str(sigma)+'_rolling_win'+str(rolling_win)+'.h5', 'w')
        file.write(mesh,'/mesh')
        file_C = XDMFFile(MPI.comm_world,'../results/PreProcess/density_'+str(sigma)+'_'+str(sigma)+'_rolling_win'+str(rolling_win)+'.xdmf');


    for t in range(total_time):
      if t%1==0:
        C_data=C_data_all[t]
        C_array=np.reshape(C_data,(-1))
        C.vector().set_local(C_array[d2v])
        C.vector().apply("insert")
        
       
      
        file_C.write(C,t);

        file.write(C,'/density_'+str(t))

