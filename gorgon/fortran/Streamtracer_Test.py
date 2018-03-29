''' Streamtracer and Connectivity Tracer Test Script '''
import sys
mod = 'C:/Users/Lars/OneDrive - Imperial College London/Code/'
if(mod not in sys.path):
    sys.path.append(mod)
import numpy as np

run_name = 'new'

# Import data

from magnetosphere.gorgon.gorgon_import import gorgon_sim
sim_dir = 'C:/Users/Lars/OneDrive - Imperial College London/Code/Shared/sample_data'
sim = gorgon_sim(sim_dir)

sim.import_timestep(-1, 'Bvec_c')

# Run streamtracer

from magnetosphere.gorgon.fieldlines import fieldlines

flines = fieldlines(10000, 0.1*sim.d[0])

nx = 20
nz = 10
x = np.linspace(sim.x[0], sim.x[-1], nx)
z = np.linspace(sim.z[0], sim.z[-1], nz)
X, Z = np.meshgrid(x, z)
Y = np.zeros_like(X)

xs = np.stack([X.ravel(), Y.ravel(), Z.ravel()]).T

flines.calc(xs, sim.arr['Bvec_c'], sim.d, sim.xc)

flines.write_vtp(run_name+'.vtp')

del flines

# Run connectivity

from magnetosphere.gorgon.connectivity import calc_connectivity 

n = [10, 20, 20]
x = np.linspace(-10, -4, n[0])
y = np.linspace(-10, 10, n[1])
z = np.linspace(-10, 10, n[2])
dx = x[1]-x[0]; dy=y[1]-y[0]; dz=z[1]-z[0]

X, Y, Z = np.meshgrid(x, y, z)

xs = np.stack([X.ravel(), Y.ravel(), Z.ravel()]).T
print(xs.shape)

top = calc_connectivity(xs, sim.arr['Bvec_c'], sim.d, sim.xc)
top = top.reshape(X.shape)

sim.write_vti_scalar(top, x-dx/2, y-dy/2, z-dz/2, name='link')
