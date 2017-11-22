from magnetosphere.fortran.null_finder import null_finder
import numpy as np
import vtk
from vtk.util import numpy_support as vtk_np

class null_points:
    def __init__(self, vec_arr):
        null_list, self.n_nulls = null_finder.find_nulls(vec_arr)
    
        self.index = null_list[:self.n_nulls]
        
        self.var = {'id': np.arange(self.n_nulls)}
        
    def cell_pos(self, sim):
    
        self.xs = np.asarray([[sim.x[ix]+0.5*sim.d[0], 
                               sim.y[iy]+0.5*sim.d[1], 
                               sim.z[iz]+0.5*sim.d[2]] for ix, iy, iz 
                                     in self.index])
    
    def write_vtp(self, fname):
        
        # Points
        pts = np.vstack([xi for xi in self.xs]).ravel()
        doubleArray = vtk_np.numpy_to_vtk(pts)
        doubleArray.SetNumberOfComponents(3)
        
        points = vtk.vtkPoints()
        points.SetData(doubleArray)
        
        # Pointdata
        
        point_arrays = self._vtk_pointData_arrays()
        
        # Polydata
        
        polyData = vtk.vtkPolyData()
        polyData.SetPoints(points)
        
        for Arr in point_arrays:
            polyData.GetPointData().AddArray(Arr)
        
        # Write to file
        
        writer = vtk.vtkXMLPolyDataWriter()
        writer.SetFileName(fname)
        writer.SetInputData(polyData)
        writer.Write()  
        
    def _vtk_pointData_arrays(self):
        
        def create_pointDataArrays(arr, name):
            
            if(len(arr[0].shape)==1):
                data = np.hstack([fi for fi in arr])
            else:
                data = np.vstack([fi for fi in arr]).ravel()
                
            data = data.astype(np.float64)
                
            doubleArray = vtk_np.numpy_to_vtk(data, deep=1)
            doubleArray.SetName(name)
            
            if(len(arr[0].shape)>1):
                doubleArray.SetNumberOfComponents(arr[0].shape[1])
        
            return doubleArray
        
        return [create_pointDataArrays(self.var[name], name) for name in self.var]