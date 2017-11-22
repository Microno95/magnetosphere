import numpy as np
from magnetosphere.fortran.connectivity_tracer import connectivity_tracer
from magnetosphere.fieldlines import fieldlines
    

def calc_end_points(x0, arr, d, xc, ns=10000, ds=None):
    if ds is None:
        ds = 0.1*d[0]
    
    xf_f, ROT_f, ns_f = connectivity_tracer.streamline_array(x0, arr, d, xc, 1, ns, ds)
    xf_r, ROT_r, ns_r = connectivity_tracer.streamline_array(x0, arr, d, xc, -1, ns, ds)
    
    return np.stack([xf_f, xf_r]).transpose(1, 0, 2)

def calc_linkage(xf, sim):
        
        def end_cat(xi, sim):
            
            
            ri = np.sqrt(np.sum(xi**2, axis=-1))
            
            i = 1
            el_x = np.logical_or(xi[:, 0]<sim.x[i], xi[:, 0]>sim.x[-i-1])
            el_y = np.logical_or(xi[:, 1]<sim.y[i], xi[:, 1]>sim.y[-i-1])
            el_z = np.logical_or(xi[:, 2]<sim.z[i], xi[:, 2]>sim.z[-i-1])
            
            el_IB = ri<1.+np.sqrt(np.sum(sim.d**2))
    
            # Check to see if streamline has gone out of bounds
            
            el_bounds = np.vstack([el_x, el_y, el_z]).any(axis=1)
            
            # Check to see if streamline has gone to the inner boundary
            
            #Closed
            if(el_IB[0] and el_IB[1]):
                link = 1
            #SW
            elif(el_bounds[0] and el_bounds[1]):
                link = 2
            # North Open
            elif(el_IB[0] and el_bounds[1]):
                link = 3
            # South Open
            elif(el_bounds[0] and el_IB[1]):
                link = 4
            # Incomplete: reaches outer boundaries
            elif(el_bounds[0] or el_bounds[1]):
                link = 5
            # Incomplete: reaches North
            elif(el_IB[0]):
                link = 6
            # Incomplete: reaches South
            elif(el_IB[1]):
                link = 7
            # Incomplete: finishes within simulation domain
            else:
                link = 0
        
            return link
        
        return np.array([end_cat(xi, sim) for xi in xf])