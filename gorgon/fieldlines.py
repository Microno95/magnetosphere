# -*- coding: utf-8 -*-
"""
Created on Wed Sep  6 14:14:46 2017

@author: Lars
"""
import numpy as np
from magnetosphere.gorgon.streamline import streamline_array
import vtk
from vtk.util import numpy_support as vtk_np

class fieldlines(streamline_array):
    def __init__(self, n_steps, step_size):
        
        streamline_array.__init__(self, n_steps, step_size, direction=0)
		
		
    
    def calc_linkage(self, sim):
        
        def end_cat(xi, sim):
            
            
            ri = np.sqrt(np.sum(xi**2, axis=-1))
            
            i = 2
            el_x = np.logical_or(xi[:, 0]<sim.x[i], xi[:, 0]>sim.x[-i])
            el_y = np.logical_or(xi[:, 1]<sim.y[i], xi[:, 1]>sim.y[-i])
            el_z = np.logical_or(xi[:, 2]<sim.z[i], xi[:, 2]>sim.z[-i])
            
            el_IB = ri<1.+sim.d[0]
    
            # Check to see if streamline has gone out of bounds
            
            el_bounds = np.vstack([el_x, el_y, el_z]).any(axis=0)
            
#            for s, f in zip(['x', 'y', 'z', 'IB', 'bounds'],
#                            [el_x, el_y, el_z, el_IB, el_bounds]):
#                print(s, f)
            
            # Check to see if streamline has gone to the inner boundary
            
            Closed = np.logical_and(el_IB[0], el_IB[1])
            SW = np.logical_and(el_bounds[0], el_bounds[1])
            Open = np.logical_xor(el_bounds[0], el_bounds[1])
            
#            for s, f in zip(['Closed', 'SW', 'Open'], [Closed, SW, Open]):
#                print(s, f)
            
            if(Closed):
                link = 1
            elif(SW):
                link = 2
            elif(Open):
                
                i = np.argmin(ri)
                if(xi[i, 2] > 0.):
                    link = 3
                elif(xi[i, 2] < 0.):
                    link = 4
            else:
                link = 0
        
            return link
        
        # Place the end points into an array
        x_ends = np.array([ [xi[0, :], xi[-1, :]] for xi in self.xs if xi.size>0])
        
        self.cell_data['link'] = np.array([end_cat(xi, sim) for xi in x_ends]) 
#        self.link = end_cat(x_ends[0], sim) 
        
    def filter(self, i):
        self.xs = self.xs[i]
        self.ROT = self.ROT[i]
        
        self.ns = self.ns[i]
        
        for s in self.var:
            self.var[s] = self.var[s][i]
            
        #self.x0 = self.x0[i]
        
        for s in self.cell_data:
            self.cell_data[s] = self.cell_data[s][i]
        