from ufl import *
from dolfin import *
import numpy as np
import csv
import pickle

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
    
print('~~~~~~~~~~~~~~~~~~~~~~~~~~~ STARTING READ WRITE ~~~~~~~~~~~~~~~~~~~~~~~~~~~')
    
with open('vsi_settings.pkl' , 'rb') as settingsFile:
    vps = pickle.load( settingsFile)    

"""
# Create mesh and define function space
#0.4375 x 0.5 x 0.5
micron_per_pix = 3.45
mag = 10
pix_stride = 100
zoom_um_per_pix = micron_per_pix/mag
zoom_um_per_pix = 0.5           #Overriding this value based on Patrick's email June 29th 
len_stride = pix_stride * zoom_um_per_pix

wells = ['%s%02d' % (r, c) for r in 'ABC' for c in range(1, 7)]
wells = wells[1:]
print('Wells are:')
print(wells)
nWells = len(wells)
"""
len_stride = vps.len_stride
numXPoints = vps.numXPoints
lastPixel = vps.lastPixel
oneD_flag = vps.oneD_flag
halfDataFlag = vps.halfDataFlag
wells = vps.wells
nWells = vps.nWells
time_gap = vps.time_gap
total_time = vps.total_time
sigma = vps.sigma
rolling_win = vps.rolling_win


#to remove
"""
oneD_flag=False
halfDataFlag = False
if oneD_flag:
    if halfDataFlag: 
    	lastPixel = 10	#For half length in 1D
    	numXPoints = 43
    else:
	    lastPixel = 19  #For full length in 1D
	    numXPoints = 76
else: #Need to write mesh refinement to do high numXPoints
    if halfDataFlag: 
    	lastPixel = 10	#For half length in 1D
    	numXPoints = 10
    else:
	    lastPixel =  42-1 #For full length in 1D
	    numXPoints = 80-1
"""

for w in range(nWells):
    wellCur = wells[w]
    point0=Point(0,0)
    point1=Point(numXPoints*len_stride,lastPixel*len_stride)
    mesh = RectangleMesh(MPI.comm_world,point0,point1,numXPoints, lastPixel)


    V = FunctionSpace(mesh, "Lagrange", 1)
    TS= TensorFunctionSpace(mesh, "CG", 1)
    C = Function(V) 
    erk = Function(V) 
    akt = Function(V) 

    d2v = dof_to_vertex_map(V)

    C_data_all=[0]*total_time
    #erk_all=[0]*total_time
    #akt_all=[0]*total_time
    dataLoadPath = vps.dataLoadPathPreProc.format(wellCur)
    for t in range(total_time):
        C_data_all[t]=np.loadtxt(dataLoadPath + vps.CdataName.format(sigma, t+1),delimiter=",")
        #erk_all[t]=np.loadtxt('../data/cell_density/erk/erk_filter_'+str(sigma)+'_t'+str(t+1)+'.dat',delimiter=",")
        #akt_all[t]=np.loadtxt('../data/cell_density/akt/akt_filter_'+str(sigma)+'_t'+str(t+1)+'.dat',delimiter=",")		

    C_data_all=central_rolling_mean(C_data_all,rolling_win) 
    #erk_all=central_rolling_mean (erk_all,rolling_win)
    #akt_all=central_rolling_mean (akt_all,rolling_win)

    savePath = vps.dataSavePathPreProc.format(wellCur)
    saveNameH5 = vps.dataSaveNameRWh5.format(sigma, sigma, rolling_win)
    saveNameXDMF = vps.dataSaveNameRWxdmf.format(sigma,sigma, rolling_win)
    file = HDF5File(MPI.comm_world, savePath + saveNameH5, 'w')
    file.write(mesh,'/mesh')
    file_C = XDMFFile(MPI.comm_world, savePath + saveNameXDMF)
    #to remove
    """
    if oneD_flag:

        file = HDF5File(MPI.comm_world, '../results/PreProcess/cell_density_1D_'+str(sigma)+'_'+str(sigma)+'_rolling_win'+str(rolling_win)+'_refine4'+'.h5', 'w')
        file.write(mesh,'/mesh')
        file_C = XDMFFile(MPI.comm_world,'../results/PreProcess/density_1D_'+str(sigma)+'_'+str(sigma)+'_rolling_win'+str(rolling_win)+'_refine4'+'.xdmf');
        #file_erk = XDMFFile(MPI.comm_world,'../results/PreProcess/erk_1D_'+str(sigma)+'_'+str(sigma)+'_rolling_win'+str(rolling_win)+'_refine4'+'.xdmf');
        #file_akt = XDMFFile(MPI.comm_world,'../results/PreProcess/akt_1D_'+str(sigma)+'_'+str(sigma)+'_rolling_win'+str(rolling_win)+'_refine4'+'.xdmf');
    else:
        file = HDF5File(MPI.comm_world, '../results/PreProcess/'+ wellCur + '/cell_density_'+str(sigma)+'_'+str(sigma)+'_rolling_win'+str(rolling_win)+'.h5', 'w')
        file.write(mesh,'/mesh')
        file_C = XDMFFile(MPI.comm_world,'../results/PreProcess/' + wellCur + '/density_'+str(sigma)+'_'+str(sigma)+'_rolling_win'+str(rolling_win)+'.xdmf');
        #file_erk = XDMFFile(MPI.comm_world,'../results/PreProcess/erk_'+str(sigma)+'_'+str(sigma)+'_rolling_win'+str(rolling_win)+'.xdmf');
        #file_akt = XDMFFile(MPI.comm_world,'../results/PreProcess/akt_'+str(sigma)+'_'+str(sigma)+'_rolling_win'+str(rolling_win)+'.xdmf');
    """
    for t in range(total_time):
      if t%1==0:
        C_data=C_data_all[t]
        C_array=np.reshape(C_data,(-1))
        C.vector().set_local(C_array[d2v])
        C.vector().apply("insert")
        
        #akt_data=akt_all[t]
        #akt_array=np.reshape(akt_data,(-1))
        #akt.vector().set_local(akt_array[d2v])
        #akt.vector().apply("insert")
        
        #erk_data=erk_all[t]
        #erk_array=np.reshape(erk_data,(-1))
        #erk.vector().set_local(erk_array[d2v])
        #erk.vector().apply("insert")

      
        file_C.write(C,t);
        #file_akt.write(akt,t);
        #file_erk.write(erk,t);

        file.write(C,'/density_'+str(t))
        #file.write(akt,'/akt_'+str(t))
        #file.write(erk,'/erk_'+str(t))
print('~~~~~~~~~~~~~~~~~~~~~~~~~~~ ENDING READ WRITE ~~~~~~~~~~~~~~~~~~~~~~~~~~~')
