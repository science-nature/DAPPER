# A mix of Evensen'2009 and Sakov'2008

# NB: Since there is no noise, and the system is stable,
#     the rmse's from this HMM go to zero as T-->infty.
#     => benchmarks largely depend on the initial error,
#     and so these absolute rmse values are not so useful
#     for quantatative evaluation of DA methods.
#     For that purpose, see mods/LA/raanes2015.py instead.

from common import *

from mods.LA.core import sinusoidal_sample, Fmat
from mods.Lorenz95.liveplotting import LP_setup

M = 1000
p = 4
jj = equi_spaced_integers(M,p)

tseq = Chronology(dt=1,dkObs=5,T=300,BurnIn=-1)

Fm = Fmat(M,c=-1,dx=1,dt=tseq.dt)
def step(x,t,dt):
  assert dt == tseq.dt
  return x @ Fm.T

# WITHOUT explicit matrix (assumes dt == dx/c):
# step = lambda x,t,dt: np.roll(x,1,axis=x.ndim-1)

Dyn = {
    'M'    : M,
    'model': step,
    'jacob': Fm,
    'noise': 0
    }

# In the animation, it can sometimes/somewhat occur
# that the truth is outside 3*sigma !!!
# Yet this is not so implausible because sinusoidal_sample()
# yields (multivariate) uniform (random numbers) -- not Gaussian.
wnum  = 25
X0 = RV(M=M, func = lambda N: sqrt(5)/10 * sinusoidal_sample(M,wnum,N))

Obs = partial_direct_obs_setup(M,jj)
Obs['noise'] = 0.01

HMM = HiddenMarkovModel(Dyn,Obs,tseq,X0,
    LP   = LP_setup(jj,conf_patch=True,conf_mult=1),
    name = os.path.relpath(__file__,'mods/'),
    )



####################
# Suggested tuning
####################
# config = EnKF('PertObs',N=100,infl=1.02)

