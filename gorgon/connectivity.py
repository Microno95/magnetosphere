import numpy as np
#from magnetosphere.gorgon.fortran.connectivity_tracer import connectivity_tracer
from magnetosphere.gorgon.fortran.streamtracer import streamtracer
    

def calc_connectivity(x0, arr, d, xc, ns=10000, ds=None, n_threads=1):
    if ds is None:
        ds = 0.1*d[0]
        
    streamtracer.ns = ns
    streamtracer.ds = ds
    streamtracer.r_IB = 1.
    streamtracer.xc = xc
    streamtracer.inner_boundary = True
    
    return streamtracer.connectivity_array(x0+xc, arr, d, n_threads)