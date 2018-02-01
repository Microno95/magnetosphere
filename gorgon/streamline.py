import numpy as np
from magnetosphere.gorgon.fortran.streamtracer import streamtracer
from scipy.interpolate import RegularGridInterpolator as interpolate
import matplotlib.pyplot as plt
import vtk
from vtk.util import numpy_support as vtk_np

#%% Single Streamline
class streamline:
    def __init__(self, n_steps, step_size, direction=0):
        self.ns = n_steps # Number of steps
        self.ns0 = n_steps # Save original number
        self.ds = step_size # Integration step size
        self.dir = direction # Integration direction of streamline
        
        self.ROT = 0 # Reason of termination
        
        self._ROT_reasons = ['Uncalculated',
                            'Out of steps',
                            'Out of domain',
                            'Isnan']
        self._dir_str = {-1: 'Reverse',
                         0: 'Both',
                         1: 'Forward'}
        
        self.var = {}
        self.cell_data = {}
        
        # Preallocate some arrays
        if(direction==1 or direction==-1):
            self.xs = np.zeros([n_steps, 3])
            self.s = np.zeros(n_steps)
        else: #If direction is both (0), need twice the space
            self.xs = np.zeros([2*n_steps, 3])
            self.s = np.zeros(2*n_steps)
            
    def __str__(self):
        
        if(type(self.ROT) == type(int())):
            ROT = self._ROT_reasons[self.ROT]
        else:
            ROT = [self._ROT_reasons[i] for i in self.ROT]
            
        direction = self._dir_str[self.dir]
        
        var_list = [s+': '+str(self.var[s].shape) for s in self.var]
        
        return 'Streamline object:\n No. steps = {}\n Step Size = {}\n Direction = {}\n Reason of Termination = {}\n Variable list: {}'.format(
                self.ns, self.ds, direction, ROT, var_list)
        
    
    # Calculate the streamline from a vector array        
    def calc(self, x0, v, d, xc):
        
        self.x0 = x0
        
        if(self.dir==1 or self.dir==-1):
            self.xs, ROT, self.ns = streamtracer.streamline(x0, v, d, xc, 
                                                             self.dir, 
                                                             self.ns,
                                                             self.ds)
            
            self.xs = self.xs[:self.ns, :]
            
        elif(self.dir==0):
            xs_f, ROT_f, ns_f = streamtracer.streamline(x0, v, d, xc, 1, 
                                                         self.ns, self.ds)
            xs_r, ROT_r, ns_r = streamtracer.streamline(x0, v, d, xc, -1, 
                                                         self.ns, self.ds)
            
            self.ROT = np.array([ROT_f, ROT_r])
            
            self.xs = np.vstack([xs_r[ns_r:0:-1, :], xs_f[:ns_f, :]])
            self.ns = self.xs.shape[0]
            
            self.ROT = np.array([ROT_f, ROT_r])
            
        self.s = np.arange(self.ns)*self.ds
        
    def reset(self, ns=None, ds=None):
        if(ns is None):
            ns = self.ns
        if(ds is None): 
            ds = self.ds
        self.__init__(ns, ds)
        
    # Interpolate for other quantities
            
    def __interp_scalar(self, x, y, z, v):
        I = interpolate((x, y, z), v, bounds_error=False)
        return I(self.xs)
    
    def _interp_vector(self, x, y, z, v):
        
        Ix = interpolate((x, y, z), v[:,:,:,0], bounds_error=False)
        Iy = interpolate((x, y, z), v[:,:,:,1], bounds_error=False)
        Iz = interpolate((x, y, z), v[:,:,:,2], bounds_error=False)
        
        return np.array([Ix(self.xs), Iy(self.xs), Iz(self.xs)]).T
    
    def interp(self, x, y, z, v, varname):
        
        if(len(v.shape)==3):
            vi = self._interp_scalar(x, y, z, v)
        elif(len(v.shape)==4):
            vi = self._interp_vector(x, y, z, v)
            
        self.var[varname] = vi
        
    # Plotting functions
    def plot(self, ax=None):
        i = 0
        if ax is None:
            fig, ax = plt.subplots()
            i = 1
            
        ax.plot(self.s, self.xs)
        
        if(i==1):
            return fig, ax
            
    def plot3D(self, ax=None):
        from mpl_toolkits.mplot3d import Axes3D
        i = 0
        if ax is None:
            fig = plt.figure()
            ax = fig.add_subplot(111, projection='3d')
            i = 1
            
        ax.plot(self.xs[:,0], self.xs[:,1], self.xs[:,2])
        
        if(i==1):
            return fig, ax
        
#%% Streamline array
class streamline_array(streamline):
    def __init__(self, n_steps, step_size, direction=0, inner_boundary=True, r_IB=1.):
        self.ns = n_steps # Number of steps
        self.ns0 = n_steps # Save original number
        self.ds = step_size # Integration step size
        self.dir = direction # Integration direction of streamline
        
        streamtracer.inner_boundary = inner_boundary
        streamtracer.r_IB = 1.
        
        
        self._ROT_reasons = ['Uncalculated',
                            'Out of steps',
                            'Out of domain',
                            'Isnan']
        self._dir_str = {-1: 'Reverse',
                         0: 'Both',
                         1: 'Forward'}
        
        self.var = {}
        self.cell_data = {}
        
    def reset(self, ns=None, ds=None):
        
        del self.xs
        del self.ROT
        
        if(ns is None):
            ns = self.ns0
        if(ds is None): 
            ds = self.ds
            
        self.__init__(ns, ds)
    
    # Calculate the streamline from a vector array
    
    def calc(self, x0, v, d, xc):
        
        self.x0 = x0
        self.n_lines = x0.shape[0]
        
        if(self.dir==1 or self.dir==-1):
            # Calculate streamlines
            self.xs, ROT, self.ns = streamtracer.streamline_array(x0, 
                                                                   v, d, xc, 
                                                                   self.dir, 
                                                                   self.ns, 
                                                                   self.ds)
            
            # Reduce the size of the array
            self.xs =np.array([xi[:ni, :] for xi, ni in zip(self.xs, self.ns)])
            
            # Save the Reason of Termination
            self.ROT = ROT
            
        elif(self.dir==0):
            # Calculate forward streamline
            xs_f, ROT_f, ns_f = streamtracer.streamline_array(x0, v, d, xc, 
                                                               1, 
                                                               self.ns, 
                                                               self.ds)
            # Calculate backward streamline
            xs_r, ROT_r, ns_r = streamtracer.streamline_array(x0, v, d, xc, 
                                                               -1, 
                                                               self.ns, 
                                                               self.ds)
            
            # Reduce the size of the arrays, and flip the reverse streamline
            xs_f = np.array([xi[:ni, :] for xi, ni in zip(xs_f, ns_f)])
            xs_r = np.array([xi[ni:1:-1, :] for xi, ni in zip(xs_r, ns_r)])
            
            # Stack the forward and reverse arrays
            self.ns =  ns_f+ns_r-1
            self.xs = np.array([np.vstack([xri, xfi]) for xri, xfi, ns in zip(xs_r, xs_f, self.ns) if ns>0])
            
            self.ROT = np.vstack([ROT_f, ROT_r]).T
            
        # Remove streamlines with zero size
        el = self.ns>0
        self.ROT = self.ROT[el] #, :
        self.ns = self.ns[el]
        
        for s in self.cell_data:
            #print(s, self.cell_data[s].shape)
            self.cell_data[s] = self.cell_data[s][el]
        
    
        
    # Interpolate for other quantities
    
    def interp(self, x, y, z, v, var_name):
        
        if(len(v.shape)==3):
            vi = self._interp_scalar(x, y, z, v)
        elif(len(v.shape)==4):
            vi = self._interp_vector(x, y, z, v)
            
        self.var[var_name] = vi
        
        self.var_names = np.array([s for s in self.var])
            
    def _interp_scalar(self, x, y, z, f):
        
        I = interpolate((x, y, z), f,)# bounds_error=False)

        xI = np.vstack(self.xs)
        fI = I(xI)

        fI = np.array(np.split(fI, np.cumsum(self.ns)))[:-1]
        
        return fI
    
    def _interp_vector(self, x, y, z, v):
        
        Ix = interpolate((x, y, z), v[:,:,:,0])#, bounds_error=False)
        Iy = interpolate((x, y, z), v[:,:,:,1])#, bounds_error=False)
        Iz = interpolate((x, y, z), v[:,:,:,2])#, bounds_error=False)

        xI = np.vstack(self.xs)
        vI = np.array([Ix(xI), Iy(xI), Iz(xI)]).T

        vI = np.array(np.vsplit(vI, np.cumsum(self.ns)))[:-1]
        
        return vI
    
        
    # Plotting methods
    
    def plot_line(self, i, ax=None, plot_names=None):
        from mpl_toolkits.mplot3d import Axes3D
        
        if(plot_names is None):
            plot_names = self.var_names
            
        n_plots = len(plot_names)+1
        
        fig = plt.figure(figsize=(12, 6))
        
        ax_line = [plt.subplot2grid((n_plots, 5), (i, 0), colspan=3) 
                        for i in range(n_plots)]
                
        ax_3D = plt.subplot2grid((n_plots, 5), (0, 3), rowspan=n_plots,
                                 colspan=2, projection='3d')
        
        self.plot3D(ax_3D, i)
        ax_3D.plot([self.xs[i][0,0]], [self.xs[i][0,1]], [self.xs[i][0,2]], '.')
        
        self.plot_linevars(i, ax=ax_line)
        
        for axi in ax_line[-2::-1]:
            axi.set(xlim = ax_line[-1].get_xlim(),
                    xticklabels=[])
        
        fig.tight_layout()
        
        return fig, ax_line, ax_3D
    
    def plot_linevars(self, i, ax=None, plot_names=None):
        
        if(plot_names is None):
            plot_names = self.var_names
        
        ret = False
        if(ax is None):
            fig, ax = plt.subplots(len(plot_names)+1, 1, sharex=True)
            ret = True
            
        s = np.linspace(0, 1, self.ns[i])*self.ns[i]*self.ds
        
        ax[0].plot(s, self.xs[i])
        mag = np.sqrt(np.sum(self.xs[i]**2, axis=1))
        ax[0].plot(s, mag)
        
        for axi, n in zip(ax[1:], plot_names):
            axi.plot(s, self.var[n][i])
            if(self.var[n][i].shape[-1]==3):
                mag = np.sqrt(np.sum(self.var[n][i]**2, axis=1))
                axi.plot(s, mag)
            axi.set_ylabel(n)
        
        if(ret):
            return fig, ax
        
        
    def plot3D(self, ax=None, i=None, seed=False, max_lines=200):
        from mpl_toolkits.mplot3d import Axes3D
        
        ret = False
        if(ax is None):
            fig = plt.figure()
            ax = fig.add_subplot(111, projection='3d')
            ret = True
        
        step = int(np.ceil(self.n_lines/max_lines))
        
        if(i is None):
            if(seed): ax.plot(self.x0[:,0], self.x0[:,1], self.x0[:,2], '.')
            for xi in self.xs[::step]:
                ax.plot(xi[:,0], xi[:,1], xi[:,2], linewidth=1)
                
        elif(type(i)==type(int())):
            if(seed): ax.plot(self.x0[i,0], self.x0[i,1], self.x0[i,2], '.')
            ax.plot(self.xs[i][:,0], self.xs[i][:,1], self.xs[i][:,2], linewidth=1)
            
        else:
            if(seed): ax.plot(self.x0[i,0], self.x0[i,1], self.x0[i,2], '.')
            for xi in self.xs[i]:
                ax.plot(xi[:,0], xi[:,1], xi[:,2], linewidth=1)
                
        ax.set(xlabel='x', ylabel='y', zlabel='z', aspect=1)
        
        if(ret):
            return fig, ax
        
    def plot_seeds3D(self, ax=None):
        from mpl_toolkits.mplot3d import Axes3D
        
        ret = False
        if(ax is None):
            fig = plt.figure()
            ax = fig.add_subplot(111, projection='3d')
            ret = True
            
        ax.plot(self.x0[:,0], self.x0[:,1], self.x0[:,2], '.')
        ax.set(xlabel='x', ylabel='y', zlabel='z', aspect=1)
                
        if(ret):
            return fig, ax
        
        #%% Write to vtp

    def write_vtp(self, fname, pts_step=None):
        if(pts_step is None):
            pts_step = int(max(1, 1./self.ds))
        
        # Points
        pts = np.vstack([xi[::pts_step] for xi in self.xs]).ravel()
        doubleArray = vtk_np.numpy_to_vtk(pts)
        doubleArray.SetNumberOfComponents(3)
        
        points = vtk.vtkPoints()
        points.SetData(doubleArray)
        
        # Cells
        
        n_pts_in_cell = np.array([len(xi[::pts_step]) for xi in self.xs])
        
        i = np.arange(np.sum(n_pts_in_cell), dtype=np.int64)
        i = np.array(np.split(i, n_pts_in_cell.cumsum())[:-1])
        
        id_array = np.array([np.hstack([ni, ii]) for ni, ii in zip(n_pts_in_cell, i)])
        
        id_array = np.hstack([ii for ii in id_array])
        
        cellArray = vtk.vtkCellArray()
        idArray = vtk_np.numpy_to_vtkIdTypeArray(id_array)
        cellArray.SetCells(len(self.ns), idArray)
        
        # Pointdata
        
        point_arrays = self._vtk_pointData_arrays(pts_step)
        
        # Cell Data
        
        cell_arrays = self._vtk_cellData_arrays(pts_step)
        
        # Polydata
        
        polyData = vtk.vtkPolyData()
        polyData.SetPoints(points)
        polyData.SetLines(cellArray)
        
        for Arr in point_arrays:
            polyData.GetPointData().AddArray(Arr)
            
        for Arr in cell_arrays:
            polyData.GetCellData().AddArray(Arr)
        
        # Write to file
        
        writer = vtk.vtkXMLPolyDataWriter()
        writer.SetFileName(fname)
        writer.SetInputData(polyData)
        writer.Write()  
        
    def _vtk_pointData_arrays(self, pts_step):
        
        def create_pointDataArrays(arr, name):
            
            if(len(arr[0].shape)==1):
                data = np.hstack([fi[::pts_step] for fi in arr])
            else:
                data = np.vstack([fi[::pts_step] for fi in arr]).ravel()
                
            data = data.astype(np.float64)
                
            doubleArray = vtk_np.numpy_to_vtk(data, deep=1)
            doubleArray.SetName(name)
            
            if(len(arr[0].shape)>1):
                doubleArray.SetNumberOfComponents(arr[0].shape[1])
        
            return doubleArray
        
        return [create_pointDataArrays(self.var[name], name) for name in self.var]
        
    def _vtk_cellData_arrays(self, pts_step):
        
        def create_cellDataArrays(arr, name):
            
            data = arr.ravel()
                
            data = data.astype(np.float64)
                
            doubleArray = vtk_np.numpy_to_vtk(data, deep=1)
            doubleArray.SetName(name)
            
            if(len(arr.shape)>1):
                doubleArray.SetNumberOfComponents(arr.shape[1])
        
            return doubleArray
        
        cell_arrays = [self.ns, self.ROT]
        cell_names = ['ns', 'ROT']
        
        for s in self.cell_data:
            cell_arrays.append(self.cell_data[s])
            cell_names.append(s)
        
        return [create_cellDataArrays(arr, name) for arr, name in zip(cell_arrays, cell_names)]
        