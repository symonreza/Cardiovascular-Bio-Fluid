#Author: Symon Reza
#Spearman & Pearson Correlation between two arreys

#from dolfin import *
import numpy
import vtk
import math
import sys
from scipy.stats import spearmanr
from scipy.stats import pearsonr
from vtk.util import numpy_support as VN

root_dir = '/Users/symon/Data/Research/Residence_Time/FinalResult/AAA/'
input_filename_1 = root_dir+ 'Final_exp_Eu_MET0.vtk'
input_filename_2 = root_dir+ 'Eulerian_Indicator_new.vtk'
#input_filename_3 = root_dir+ 'Surface_PRT.vtk'
#input_filename_4 = root_dir+ 'Surface_PRT.vtk'
#output_filename = root_dir + 'SC_PointRT_MET.vtk'
flag=1

################
#Read data for 1st file
reader = vtk.vtkUnstructuredGridReader()
reader.SetFileName(input_filename_1)
reader.Update()
data = reader.GetOutput()
n_points = data.GetNumberOfPoints()

gnID  = VN.vtk_to_numpy( data.GetPointData().GetArray('GlobalNodeID'))#Node_ID from ParaView
RT1 = VN.vtk_to_numpy( data.GetPointData().GetArray('concentration'))
print ('RT1'), RT1

#Read data for 2nd file
reader = vtk.vtkUnstructuredGridReader()
reader.SetFileName(input_filename_2)
reader.Update()
data1= reader.GetOutput()
n_points1 = data1.GetNumberOfPoints()

gnID1  = VN.vtk_to_numpy(data1.GetPointData().GetArray('GlobalNodeID'))#Node_ID from ParaView
RT2 = VN.vtk_to_numpy( data1.GetPointData().GetArray('RT_mobility'))
print ('RT2'), RT2
count2=0


#Include a Flag (1 for Eulerian Indicator & Pointwise RT)
if (flag==1):
 meshtosurf = vtk.vtkDataSetSurfaceFilter()
 meshtosurf.SetInputData(data)
 meshtosurf.Update()
 data2 = meshtosurf.GetOutput()
 n_points2 = data2.GetNumberOfPoints()
 n_p_d=n_points-n_points2
 gnID2  = VN.vtk_to_numpy(data2.GetPointData().GetArray('GlobalNodeID'))#Node_ID from ParaView
 print ('gnID2'), gnID2
 RT3 = VN.vtk_to_numpy( data2.GetPointData().GetArray('concentration'))
 print ('RT3'), RT3
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


 meshtosurf = vtk.vtkDataSetSurfaceFilter()
 meshtosurf.SetInputData(data1)
 meshtosurf.Update()
 data3 = meshtosurf.GetOutput()
 n_points3 = data3.GetNumberOfPoints()
 n_p_d1=n_points1-n_points3
 gnID3  = VN.vtk_to_numpy(data3.GetPointData().GetArray('GlobalNodeID') )#Node_ID from ParaView
 print ('gnID3'), gnID3
 RT4 = VN.vtk_to_numpy( data3.GetPointData().GetArray('RT_mobility') )
 print ('RT4'), RT4
 print ('n_points3'), n_points3
 tag1= numpy.zeros(n_points1)
 arrey1= numpy.zeros(n_points1-n_points3)
 D= -1
 for i in xrange (n_points3):
    check=0
    for j in xrange (n_points1):
        if gnID3[i] == gnID1[j]:
            tag1[j]=1
            check=1
        if check==1:
            break
 for i in xrange (n_points1):
        if tag1[i]==0:
         D=D+1
         arrey1[D]= RT2[i]
        if D==n_p_d1:
          break
 print ('arrey1'),arrey1
 P= arrey
 Q=arrey1
else:
 #Mesh check
 n_points = data.GetNumberOfPoints()
 L_pt1 = numpy.zeros((n_points, 3))
 #sorting the node and elements wall mesh
 for i in xrange(n_points):
   clpts1 = data.GetPoint(i)
   L_pt1[i,0] = clpts1[0]
   L_pt1[i,1] = clpts1[1]
   L_pt1[i,2] = clpts1[2]
    
 n_points1 = data1.GetNumberOfPoints()
 L_pt2 = numpy.zeros((n_points, 3))
    #sorting the node and elements wall mesh
 for i in xrange(n_points1):
   clpts2 = data1.GetPoint(i)
   L_pt2[i,0] = clpts2[0]
   L_pt2[i,1] = clpts2[1]
   L_pt2[i,2] = clpts2[2]
    
    
 if (n_points != n_points1):
  print "Error: number of pts do not match"
  sys.exit()
                
 eps=1e-5
 for i in xrange(n_points):
  if ( abs(L_pt1[i,0] - L_pt2[i,0] ) > eps or abs(L_pt1[i,1] - L_pt2[i,1]) > eps or abs(L_pt1[i,2] -L_pt2[i,2]) > eps):
   print "Error: coordinates dont match"
   A= abs(L_pt1[i,0] - L_pt2[i,0] )
   B= abs(L_pt1[i,1] - L_pt2[i,1])
   C= abs(L_pt1[i,2] - L_pt2[i,2])
   print 'Difference-x', A
   print 'Difference-y', B
   print 'Difference-z', C
   sys.exit()
 P= RT1
 Q= RT2
rank_corr_S =spearmanr(P,Q)
print 'Spearmans rank corr value:', rank_corr_S

rank_corr_P =pearsonr(P,Q)
print 'Pearsons rank corr value:', rank_corr_P
 
