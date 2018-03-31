''' Streamtracer and Connectivity Tracer Test Script '''
import sys
mod = 'C:/Users/Lars/OneDrive - Imperial College London/Code/'
#mod = 'D:/Lars/OneDrive - Imperial College London/Code/'
if(mod not in sys.path):
    sys.path.append(mod)
import numpy as np
from datetime import datetime

run_name = 'new'
s = 'Bvec_c'

# Import data

from magnetosphere.gorgon.gorgon_import import gorgon_sim
sim_dir = 'C:/Users/Lars/OneDrive - Imperial College London/Code/Shared/sample_data'
#sim_dir = 'D:/Lars/BSM_9_restart'
sim = gorgon_sim(sim_dir)


startTime = datetime.now()
sim.import_timestep(-1, s)
print('Gorgon import took {}'.format(datetime.now()-startTime))

# Run streamtracer

from magnetosphere.gorgon.fieldlines import fieldlines

flines = fieldlines(10000, 0.1*sim.d[0])

step = 5
x = sim.x[::step]+0.5*sim.d[0]
y = sim.y[::step]+0.5*sim.d[1]
z = sim.z[::step]+0.5*sim.d[2]
X, Y, Z = np.meshgrid(x, y, z)

x = sim.x[100::-step]+0.5*sim.d[0]
y = sim.y[50::-step]+0.5*sim.d[1]
z = sim.z[50::-step]+0.5*sim.d[2]
Y, X, Z = np.meshgrid(y, x, z)

# x = sim.x+0.5*sim.d[0]
# X = x
# Y = np.zeros_like(X)
# Z = np.zeros_like(X)

xs = np.stack([X.ravel(), Y.ravel(), Z.ravel()]).T
print(xs.shape)

startTime = datetime.now()
flines.calc(xs, sim.arr[s], sim.d, sim.xc-0.5*sim.d)
print('Streamtracer took {}'.format(datetime.now()-startTime))
flines.calc_linkage(sim)

startTime = datetime.now()
flines.write_vtp(run_name+'.vtp')
print('vtp writing took {}'.format(datetime.now()-startTime))

del flines

# Run connectivity

from magnetosphere.gorgon.connectivity import calc_connectivity

x = sim.x[:100]+0.5*sim.d[0]
y = sim.y[:50]+0.5*sim.d[1]
z = sim.z[:50]+0.5*sim.d[2]
Y, X, Z = np.meshgrid(y, x, z)

xs = np.stack([X.ravel(), Y.ravel(), Z.ravel()]).T
print(xs.shape)

startTime = datetime.now()
top = calc_connectivity(xs, sim.arr[s], sim.d, sim.xc-0.5*sim.d)
print('Connectivity took {}'.format(datetime.now()-startTime))
top = top.reshape(X.shape)

print(X.shape, top.shape)

dx = x[1]-x[0]; dy=y[1]-y[0]; dz=z[1]-z[0]
sim.write_vti_scalar(top, x, y, z, name='link')
