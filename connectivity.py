import numpy as np
from magnetosphere.fortran.connectivity_tracer import connectivity_tracer
from magnetosphere.fieldlines import fieldlines
    

def calc_end_points(x0, arr, d, xc, ns=10000, ds=None):
    if ds is None:
        ds = 0.1*d[0]
    
    return connectivity_tracer.streamline_array(x0, arr, d, xc, ns, ds)