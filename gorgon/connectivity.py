import numpy as np
from magnetosphere.gorgon.fortran.connectivity_tracer import connectivity_tracer
from magnetosphere.gorgon.fieldlines import fieldlines
    

def calc_connectivity(x0, arr, d, xc, ns=10000, ds=None):
    if ds is None:
        ds = 0.1*d[0]
        
    connectivity_tracer.ns = ns
    connectivity_tracer.ds = ds
    connectivity_tracer.r_IB = 1.
    connectivity_tracer.inner_boundary = 1.
    
    return connectivity_tracer.connectivity_array(x0, arr, d, xc)