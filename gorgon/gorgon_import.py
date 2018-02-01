import glob
import os
import re
import vtk
import numpy as np
from vtk.util import numpy_support as vtk_np
''' gorgon_import class for importing pvti files generated by Gorgon 
    
    Any issues, feel free to email me at: lars.mejnertsen10@imperial.ac.uk
'''

class gorgon_sim:
    def __init__(self, data_directory=None, output=False):
        self.fdir = data_directory
        
        self.arr = {}
        
        if(data_directory is not None):
            files = glob.glob(self.fdir+r'/*.pvti')
            
            file_names = [f.split(os.sep)[-1] for f in files]
            
            self.index = np.array([re.split('[_]', f)[0] for f in file_names])
            self.index = np.unique(self.index)[0]
            i = len(self.index)
            
            out = np.array([re.split('[.-]', f[i+1:]) for f in file_names])
            
            self.arr_names = np.unique(out[:, 0])
            self.times = np.unique(out[1:, 1:-1])
            s = np.argsort([int(i) for i in self.times])
            self.times = self.times[s]
            
            for s in self.arr_names:
                if(s not in ['Avec', 'Evec']):
                    break
            
            file_name = self.create_filename(s, self.times[0])
            
            self.import_space(file_name)
            
            self.time = self.times[0]
            
            if(output):
                print(self.index)
                print(self.arr_names)
                print(self.times)
                print(self.time)
        
    def create_filename(self, arr_name, time=None, suffix='pvti', nodir=False):
       
        if(time is None):
            time = self.time
        
        if(type(time)!=type(str())):
            time = str(time)
        
        string = '{}{}_{}-{}.'+ suffix
        
        fdir = self.fdir+os.sep
        if(nodir):
            fdir=''
            
        return string.format(fdir,
                             self.index,
                             arr_name,
                             time)
        
    def import_space(self, filename):
        if(not os.path.isfile(filename)):
            print("Can't find file: " + filename)
        else:
    
            reader = vtk.vtkXMLPImageDataReader()
            reader.SetFileName(filename)
            reader.Update()
            
            data = reader.GetOutput()
    
            bounds = np.array(data.GetBounds())
            spacing = np.array(data.GetSpacing())
            
            self.d = spacing
            self.xc = -np.array([bounds[0], bounds[2], bounds[4]])
            
            self.x = np.arange(bounds[0], bounds[1], spacing[0])
            self.y = np.arange(bounds[2], bounds[3], spacing[1])
            self.z = np.arange(bounds[4], bounds[5], spacing[2])
            
    def import_timestep(self, i, arr_names=None, delete_existing=True):
        self.time = self.times[i]
        if(arr_names is None):
            arr_names = self.arr_names
            
        for s in self.arr_names:
            if(s in arr_names):
                filename = self.create_filename(s)
                self.import_vtk(filename, s)
            elif(delete_existing and s in self.arr):
                del self.arr[s]
        
    def import_vtk(self, filename, varname, log=False):
        if(not os.path.isfile(filename)):
            print("Can't find file: " + filename)
        else:
    
            reader = vtk.vtkXMLPImageDataReader()
            reader.SetFileName(filename)
            reader.Update()
            
            data = reader.GetOutput()
            dim = data.GetDimensions()
            
            v = vtk_np.vtk_to_numpy(data.GetCellData().GetArray(0))
            
            if(log):
                v = 10**v
            
            n_comp = data.GetCellData().GetArray(0).GetNumberOfComponents()
            
            vec = [int(i - 1) for i in dim]
            
            if(n_comp>1):
                vec.append(n_comp)
            
            self.arr[varname] = v.reshape(vec,order='F')
         
    def set_3D_axlim(self, ax):
        ax.set(xlim=[self.x[0], self.x[-1]],
               ylim=[self.y[0], self.y[-1]],
               zlim=[self.z[0], self.z[-1]])
        ax.set_aspect(1)

    def write_vti_scalar(self, arr, x=None, y=None, z=None, name=None, output_dir=''):
        ''' 
        Writes a 3D numpy array into a VTI file. Only works for scalar values. Must be uniform spacing.
        arr: the data you want to write.
            If str: writes from sim.arr
            If numpy array: writes the numpy array. Assumes the same size as sim unless x, y and z are specified
            If dict: writes multiple numpy arrays in dict to vti. Assumes the same size as sim unless x, y and z are specified
        x, y, z: the x, y and z position array,
        name: Name of the file/data array. Needed for np.array or dict arr input
        output_dir: Specify a different output folder.    
        '''
        from evtk.hl import imageToVTK 
        if(output_dir is not ''):
            output_dir += '/'
        
        if(x is None):
            x = self.x
            y = self.y
            z = self.z
        
        # If x, y, and z are different from the array
        d = np.array([x[1]-x[0], y[1]-y[0], z[1]-z[0]])
        xc = np.array([x[0], y[0], z[0]]) #-0.5*d
        
        if(type(arr) is type('')):        
            imageToVTK(output_dir+self.index+'_'+arr+'-'+self.time,
                       cellData={arr: self.arr[arr]},
                       origin=xc.tolist(), spacing=d.tolist())
            
        else:
            if name is None:
                print('For arr of type np.array or dict, you need to specify a name')
                return
            
            fname = output_dir+self.index+'_'+name+'-'+self.time
            
            if type(arr) is type(dict()):
                imageToVTK(fname,
                           cellData=arr,
                           origin=xc.tolist(), spacing=d.tolist())
            
            elif type(arr) is type(np.array([])):
                imageToVTK(fname,
                           cellData={name: arr},
                           origin=xc.tolist(), spacing=d.tolist())
                
            else:
                print('arr must be of type str, np.array or dict')
            