import numpy as np
from magnetosphere.gorgon.fortran.connectivity_tracer import connectivity_tracer
from magnetosphere.gorgon.fieldlines import fieldlines
    

def calc_connectivity(x0, arr, d, xc, ns=10000, ds=None):
    if ds is None:
        ds = 0.1*d[0]
        
    connectivity_tracer.ns = ns
    connectivity_tracer.ds = ds
    
    return connectivity_tracer.streamline_array(x0, arr, d, xc)