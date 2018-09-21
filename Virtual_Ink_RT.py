#Author: Symon Reza
#Virtual Ink Residence Time Comutation
import numpy
import vtk
import math
import sys
from scipy.stats import spearmanr
from scipy.stats import pearsonr
from vtk.util import numpy_support as VN


fieldname = 'f_20'
T_first = 0
T_delta = 1
T = 100
t_index= T_first-T_delta
input_filename =  '/Users/symon/Data/Research/Residence_Time/AAA/AAA_new_VI_'
output_root = '/Users/symon/Data/Research/Residence_Time/FinalResult/AAA/'
output_filename = output_root + 'VI_AAA_RT_test02_100.vtk'

################
#Read concentration


#Assign a constant to the nodes
#reader = vtk.vtkUnstructuredGridReader()
#reader.SetFileName(input_filename_2)
#reader.Update()
#data1= reader.GetOutput()
#n_points1 = data1.GetNumberOfPoints()

#gnID1  = VN.vtk_to_numpy(data1.GetPointData().GetArray('GlobalNodeID') )#Node_ID from ParaView
#Concentration1 = VN.vtk_to_numpy( data1.GetPointData().GetArray('concentration') )
#count2=0
#Concentration1.assign(interpolate(Constant(0.002), n_points))
#print ('Concentration1'), Concentration1

#Define time
#t_start= 0
#t_stop= 3.8
#n_tsteps= 1000*4
#tsteps_per_file= 4
#dt = (t_stop - t_start) / n_tsteps

n_points= 1143405
tag= numpy.zeros(n_points)
array_RT = numpy.zeros(n_points)

for t in xrange(T):
    print 't=', t
    t_index = t_index + T_delta
    input_filename2 = input_filename + str(t_index+T_delta) + '.vtk'
    print 'Loading', input_filename2
    # reader = vtk.vtkDataSetReader()
                
    #read the next t.s. to calculate time derivative
    #input_filename3 = input_filename + str(t_index+T_delta) + '.vtk'
    #print 'Loading', input_filename3
    reader = vtk.vtkUnstructuredGridReader()
    reader.SetFileName(input_filename2)
    reader.Update()
    data = reader.GetOutput()
    n_points = data.GetNumberOfPoints()
    #gnID  = VN.vtk_to_numpy( data.GetPointData().GetArray('GlobalNodeID') )#Node_ID from ParaView
    Concentration = VN.vtk_to_numpy( data.GetPointData().GetArray('f_20') )
    print ('Concentration'), Concentration
    #print ('Concentration'), Concentration

    for i in xrange (n_points):
        if (Concentration[i] < 0.02 and tag[i] ==0) :
           tag[i]=1
           array_RT[i] = (t+1.0) * 0.1053
print ('array_RT'), array_RT

output_vtk = VN.numpy_to_vtk(array_RT)
output_vtk.SetName('VI_RT')
data.GetPointData().AddArray(output_vtk)
#data.GetPointData().RemoveArray('f_20')
myoutput = vtk.vtkDataSetWriter()
myoutput.SetInputData(data)
myoutput.SetFileName(output_filename)
myoutput.Write()
    
 
