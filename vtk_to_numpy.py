# -*- coding: utf-8 -*-
import vtk
from vtk.util.numpy_support import vtk_to_numpy
import numpy as np
import os

def VTK_import_files(filedir,ctime,varlist,space_import=True):
    if(space_import):
        x,y,z = VTK_import_space(filedir+varlist[0]+'-'+ctime+'.vti')

    D = {s: VTK_import_array(filedir+s+'-'+ctime+'.vti','x44_'+s) for s in varlist}
    
    if(space_import):
        return x,y,z,D
    else:
        return D
    
def VTK_import_file(filename,varname):
    if(os.path.isfile(filename)):
        print('File',filename,varname)
    else:
        print('File',filename + ' not found')
        return
    
    # Load file
    
    #varname = filename[:4] + varname

    reader = vtk.vtkXMLImageDataReader()
    reader.SetFileName(filename)
    reader.Update()
    reader.GetNumberOfCells()
    
    data = reader.GetOutput()
    
    ''' Import space '''

    dim = np.asarray(data.GetDimensions())
    c = np.asarray(data.GetOrigin())
    d = np.asarray(data.GetSpacing())

    x = np.arange(dim[0])*d[0]+c[0]
    y = np.arange(dim[1])*d[1]+c[1]
    z = np.arange(dim[2])*d[2]+c[2]

    x = 0.5*(x[1:]+x[:-1])
    y = 0.5*(y[1:]+y[:-1])
    z = 0.5*(z[1:]+z[:-1])
    
    # Import array
    
    cdata = data.GetCellData()
    N_comps = cdata.GetNumberOfComponents()
    
    v = vtk_to_numpy(data.GetCellData().GetArray(varname))
    
    vec = [int(i-1) for i in dim]
    if(N_comps>1):
        vec.append(N_comps)
    v = v.reshape(vec,order='F') 
    
    return x,y,z,v
        

def VTK_import_space(filename):
    if(os.path.isfile(filename)):
        print('Space',filename)
    else:
        print('Space',filename + ' not found')
        return
    
    reader = vtk.vtkXMLImageDataReader()
    reader.SetFileName(filename)
    reader.Update()
    reader.GetNumberOfCells()
    
    data = reader.GetOutput()
    
    ''' Import space '''

    dim = np.asarray(data.GetDimensions())
    c = np.asarray(data.GetOrigin())
    d = np.asarray(data.GetSpacing())

    x = np.arange(dim[0])*d[0]+c[0]
    y = np.arange(dim[1])*d[1]+c[1]
    z = np.arange(dim[2])*d[2]+c[2]

    x = 0.5*(x[1:]+x[:-1])
    y = 0.5*(y[1:]+y[:-1])
    z = 0.5*(z[1:]+z[:-1])
    
    return x,y,z


def VTK_import_array(filename,varname):
    if(os.path.isfile(filename)):
        print('File',filename,varname)
    else:
        print('File',filename + ' not found')
        return
    
    #varname = filename[:4] + varname

    reader = vtk.vtkXMLImageDataReader()
    reader.SetFileName(filename)
    reader.Update()
    reader.GetNumberOfCells()
    
    data = reader.GetOutput()
    dim = data.GetDimensions()
    
    cdata = data.GetCellData()
    N_comps = cdata.GetNumberOfComponents()
    
    v = vtk_to_numpy(data.GetCellData().GetArray(varname))
    
    vec = [int(i-1) for i in dim]
    if(N_comps>1):
        vec.append(N_comps)
    v = v.reshape(vec,order='F') 
    
    return v

def VTK_import_scalar_file(filename,varname):
    print('s',filename,varname)

    reader = vtk.vtkXMLImageDataReader()
    reader.SetFileName(filename)
    reader.Update()
    reader.GetNumberOfCells()
    
    data = reader.GetOutput()
    dim = data.GetDimensions()
    
    cdata = data.GetCellData()
    
    v = vtk_to_numpy(cdata.GetArray(varname))
    
    vec = [int(i-1) for i in dim]
    v = v.reshape(vec,order='F')
    
    return v

def VTK_import_vector_file(filename,varname):
    print('v',filename,varname)
    reader = vtk.vtkXMLImageDataReader()
    reader.SetFileName(filename)
    reader.Update()
    reader.GetNumberOfCells()
    
    data = reader.GetOutput()
    dim = data.GetDimensions()
    
    cdata = data.GetCellData()
    N_comps = cdata.GetNumberOfComponents()
    
    v = vtk_to_numpy(data.GetCellData().GetArray(varname))
    
    vec = [int(i-1) for i in dim]
    vec.append(N_comps)
    v = v.reshape(vec,order='F')   
    
    return v