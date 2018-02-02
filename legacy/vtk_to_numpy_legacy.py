# -*- coding: utf-8 -*-
import vtk
from vtk.util import numpy_support as vtk_np
import numpy as np
import os

def import_space(filename):
	if(not os.path.isfile(filename)):
            print("Can't find file: "+filename)
	else:

		reader = vtk.vtkXMLImageDataReader()
		reader.SetFileName(filename)
		reader.Update()
		reader.GetNumberOfCells()
		
		data = reader.GetOutput()
		
		dim = np.asarray(data.GetDimensions())
		
		x = np.zeros(data.GetNumberOfPoints())
		y = x.copy()
		z = x.copy()
		
		for i in range(data.GetNumberOfPoints()):
				x[i],y[i],z[i] = data.GetPoint(i)
		
		x = x.reshape(dim,order='F')
		x = 0.5*(x[1:,0,0]+x[:-1,0,0])
		y = y.reshape(dim,order='F')
		y = 0.5*(y[0,1:,0]+y[0,:-1,0])
		z = z.reshape(dim,order='F')
		z = 0.5*(z[0,0,1:]+z[0,0,:-1])
		
		return x,y,z

def import_scalar(filename,varname):
	if(not os.path.isfile(filename)):
            print("Can't find file: "+filename)
	else:

		reader = vtk.vtkXMLImageDataReader()
		reader.SetFileName(filename)
		reader.Update()
		reader.GetNumberOfCells()
		
		data = reader.GetOutput()
		dim = data.GetDimensions()
		
		#vec = list(dim)
		vec = [int(i-1) for i in dim]
		
		v = vtk_np.vtk_to_numpy(data.GetCellData().GetArray(varname))
		v = v.reshape(vec,order='F')
		
		return v

def import_vector(filename,varname):
	if(not os.path.isfile(filename)):
            print("Can't find file: "+filename)
	else:

		reader = vtk.vtkXMLImageDataReader()
		reader.SetFileName(filename)
		reader.Update()
		reader.GetNumberOfCells()
		
		data = reader.GetOutput()
		dim = data.GetDimensions()
		
		#vec = list(dim)
		vec = [int(i-1) for i in dim]
		vec.append(3)
		
		v = vtk_np.vtk_to_numpy(data.GetCellData().GetArray(varname))
		v = v.reshape(vec,order='F')
		
		return v