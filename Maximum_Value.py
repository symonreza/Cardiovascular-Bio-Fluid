#Author: Symon Reza
#Maximum Value of an arrey

#from dolfin import *
import numpy
import vtk
import math
import sys
from scipy.stats import spearmanr
from scipy.stats import pearsonr
from vtk.util import numpy_support as VN

root_dir = '/Users/symon/Data/Research/Residence_Time/FinalResult/AAA/'
input_filename_1 = root_dir+ 'Eulerian_Indicator_new.vtk'
#input_filename_3 = root_dir+ 'Surface_PRT.vtk'
#input_filename_4 = root_dir+ 'Surface_PRT.vtk'
#output_filename = root_dir + 'SC_PointRT_MET.vtk'


################
#Read data for 1st file
reader = vtk.vtkUnstructuredGridReader()
reader.SetFileName(input_filename_1)
reader.Update()
data = reader.GetOutput()
n_points = data.GetNumberOfPoints()

gnID  = VN.vtk_to_numpy( data.GetPointData().GetArray('GlobalNodeID'))#Node_ID from ParaView
RT1 = VN.vtk_to_numpy( data.GetPointData().GetArray('RT_mobility'))
print ('RT1'), RT1
print('Maximum Values:'),numpy.max(RT1)
print ('99.9 percentile'), numpy.percentile (RT1, 99.9)
print ('99.99 percentile'), numpy.percentile (RT1, 99.99)
#excluding wall
meshtosurf = vtk.vtkDataSetSurfaceFilter()
meshtosurf.SetInputData(data)
meshtosurf.Update()
data2 = meshtosurf.GetOutput()
n_points2 = data2.GetNumberOfPoints()
n_p_d=n_points-n_points2
gnID2  = VN.vtk_to_numpy(data2.GetPointData().GetArray('GlobalNodeID'))#Node_ID from ParaView
RT2 = VN.vtk_to_numpy( data2.GetPointData().GetArray('RT_mobility'))
print ('RT2'), RT2
print ('n_points2'), n_points2
tag= numpy.zeros(n_points)
arrey= numpy.zeros(n_points-n_points2)
C=-1
counter = 0
for i in xrange (n_points2):
     check=0
     for j in xrange (n_points):
         if gnID2[i] == gnID[j]:
            tag[j]=1
            check=1
            counter = counter + 1
         if check==1:
           break
print 'number of surface nodes detected', counter
for i in xrange (n_points):
        if tag[i]==0:
         C=C+1
         arrey[C]= RT1[i]
         #if C==n_p_d:
         # break
print ('arrey'),arrey
print('Maximum Values:'),numpy.max(arrey)
print ('99.9 percentile'), numpy.percentile(arrey, 99.9)
print ('99.99 percentile'), numpy.percentile (arrey, 99.99)
