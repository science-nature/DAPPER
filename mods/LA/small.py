# Smaller version.
# => convenient for debugging, e.g., scripts/test_iEnKS.py

from common import *

from mods.LA.core import Fmat, homogeneous_1D_cov
from mods.Lorenz95.liveplotting import LP_setup

tseq = Chronology(dt=1,dkObs=5,T=300,BurnIn=-1)

m = 100

# def step(x,t,dt):
  # return np.roll(x,1,axis=x.ndim-1)
Fm = Fmat(m,-1,1,tseq.dt)
def step(x,t,dt):
  assert dt == tseq.dt
  return x @ Fm.T

f = {
    'm': m,
    'model': step,
    'noise': 0
    }

X0 = GaussRV(C=homogeneous_1D_cov(m,m/8,kind='Gauss'))

p  = 4
jj = equi_spaced_integers(m,p)
h  = partial_direct_obs_setup(m,jj)
h['noise'] = 0.01

 
HMM = HiddenMarkovModel(f,h,tseq,X0,
    name = os.path.relpath(__file__,'mods/'),
    LP   = LP_setup(jj),
    )


####################
# Suggested tuning
####################
# cfgs += EnKF('PertObs',N=16 ,infl=1.02)
# cfgs += EnKF('Sqrt'   ,N=16 ,infl=1.0)

