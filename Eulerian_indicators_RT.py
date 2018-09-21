#Author: Symon Reza
#Eulerian Indicator Residence Time computation
import numpy
import sys
import vtk
import glob
import re
import math
from vtk.util import numpy_support as VN

def Eulerian_Indicator (input_filename, output_filename, fieldname, T_first, T_delta, T ):

   t_index = T_first - T_delta
   for t in xrange(T):
     print 't=', t
     t_index = t_index + T_delta
     input_filename2 = input_filename + str(t_index) + '.vtu'
     print 'Loading', input_filename2
     reader = vtk.vtkXMLUnstructuredGridReader()
     reader.SetFileName(input_filename2)
     reader.Update()
     data = reader.GetOutput()
     n_points = data.GetNumberOfPoints()
     WSS = VN.vtk_to_numpy(data.GetPointData().GetArray('velocity'))

                # reader = vtk.vtkDataSetReader()
     
     #read the next t.s. to calculate time derivative
     input_filename3 = input_filename + str(t_index+T_delta) + '.vtu'
     print 'Loading', input_filename3
     reader = vtk.vtkXMLUnstructuredGridReader()
     reader.SetFileName(input_filename3)
     reader.Update()
     data2 = reader.GetOutput()
     WSS2 = VN.vtk_to_numpy(data2.GetPointData().GetArray('velocity'))

     if t==0:  #allocate once
       wss_mag = numpy.zeros((n_points, 1))
       wss_mag2 = numpy.zeros((n_points, 1))
       RROC = numpy.zeros((n_points, 1))
       v_fast_vec = numpy.zeros((n_points, 3))
       v_slow = numpy.zeros((n_points, 1))
       v_sum = numpy.zeros((n_points, 1))
       v_sum_multiple = numpy.zeros((n_points, 1))
       RT = numpy.zeros((n_points, 1))
       WSS_difference = numpy.zeros((n_points, 3))
      
     for i in xrange(n_points):
        if (i% 100000 ==0):
            print '%done:', float(i)/ float(n_points) * 100.
        wss_mag[i] = math.sqrt((WSS[i,0]**2)+(WSS[i,1]**2)+(WSS[i,2]**2))
        wss_mag2[i] = math.sqrt((WSS2[i,0]**2)+(WSS2[i,1]**2)+(WSS2[i,2]**2))
        #WSS time derivative (fwd difference)
        WSS_difference[i,0] = WSS2[i,0] - WSS[i,0]
        WSS_difference[i,1] = WSS2[i,1] - WSS[i,1]
        WSS_difference[i,2] = WSS2[i,2] - WSS[i,2]
        WSS_difference_mag = math.sqrt( WSS_difference[i,0]**2 + WSS_difference[i,1]**2 + WSS_difference[i,2]**2 )
        RROC[i] = RROC[i] + abs( WSS_difference_mag) / ( abs(wss_mag[i]) + abs(wss_mag2[i]) + 1e-16 )  #relative time derivative

        v_fast_vec[i,0] = v_fast_vec[i,0] + WSS[i,0]
        v_fast_vec[i,1] = v_fast_vec[i,1] + WSS[i,1]
        v_fast_vec[i,2] = v_fast_vec[i,2] + WSS[i,2]
        v_slow[i] = v_slow[i] + wss_mag[i]

    

   print 'Final processing...'
   #delta_t = T_cycle / (T -1)

   #Mobility
   for i in xrange(n_points):
     RROC[i] = RROC[i] / T
     v_fast = math.sqrt( v_fast_vec[i,0]**2 + v_fast_vec[i,1]**2 + v_fast_vec[i,2]**2  )  #Like TAWSS_vec
     #v_sum[i] =  ( 1 - math.exp(-RROC[i]* T_c ) )*v_fast /T  + math.exp( -RROC[i] * T_c )* v_slow[i]/T #!!!! Divided by T so that values not very large
     v_sum[i] =  ( 1 - math.exp(-RROC[i] ) )*v_fast / T + math.exp( -RROC[i] )* v_slow[i]  / T
     RT[i]  = 1.0 /  ( v_sum[i] + 1e-16)
  



   theta_vtk = VN.numpy_to_vtk(RT)
   theta_vtk.SetName('RT_mobility')
   data.GetPointData().AddArray(theta_vtk)

   myoutput = vtk.vtkDataSetWriter()
   myoutput.SetInputData(data)
   myoutput.SetFileName(output_filename)
   myoutput.Write()


   print 'Done!'


if __name__ == "__main__":
 


 fieldname = 'velocity'
 input_filename = '/Users/symon/Data/Research/Residence_Time/ICA42/Velocity/ICA42_vel_c_'
 output_filename = '/Users/symon/Data/Research/Residence_Time//FinalResult/ICA/Eulerian_Indicator_new_ICA.vtk'
 T_first = 0  #index of first file
 T_delta = 1
 T = 50
 #T_c = 0 #Characteristic time scale for significant transport
 #T_cycle = 0.95  #cardiac cycle
 Eulerian_Indicator(input_filename, output_filename, fieldname, T_first, T_delta, T )
  

