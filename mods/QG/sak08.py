# Reproduce results from Sakov and Oke "DEnKF" paper from 2008.

from common import *

from mods.QG.core import step, dt, nx, ny, m, square, sample_filename, generate_diags


# As specified in core.py: dt = 4*1.25 = 5.0.
#t = Chronology(dt=dt,dkObs=1,T=300,BurnIn=100) # Sakov: repeat 10 times
#t = Chronology(dt=dt,dkObs=1,T=600,BurnIn=20)
t = Chronology(dt=dt,dkObs=1,T=100,BurnIn=20)
#t = Chronology(dt=dt,dkObs=1,T=20,BurnIn=10)
#
# Considering that I have 8GB on the Mac, and the estimate:
# ≈ (8 bytes/float)*(129² float/stat)*(7 stat/k) * K,
# it should be possible to have an experiment with maximum K≈8000.


f = {
    'm'    : m,
    'model': step,
    'noise': 0,
    }

X0 = RV(m=m,file=sample_filename)


############################
# Observation settings
############################
"""
generate diagonal tracks:

Indications for one direction-tracks. Of course this double in case of cross-crossing observations (change the angle).

diags:2, lighten_factor:1, obs:131
diags:2, lighten_factor:2, obs:67
diags:2, lighten_factor:3, obs:43
diags:2, lighten_factor:4, obs:34
diags:2, lighten_factor:5, obs:26
diags:2, lighten_factor:6, obs:22
diags:4, lighten_factor:1, obs:261
diags:4, lighten_factor:2, obs:133
diags:4, lighten_factor:3, obs:86
diags:4, lighten_factor:4, obs:67
diags:4, lighten_factor:5, obs:52
diags:4, lighten_factor:6, obs:44
diags:8, lighten_factor:1, obs:521
diags:8, lighten_factor:2, obs:265
diags:8, lighten_factor:3, obs:173
diags:8, lighten_factor:4, obs:133
diags:8, lighten_factor:5, obs:104
diags:8, lighten_factor:6, obs:88
diags:16, lighten_factor:1, obs:1041
diags:16, lighten_factor:2, obs:529
diags:16, lighten_factor:3, obs:346
diags:16, lighten_factor:4, obs:265
diags:16, lighten_factor:5, obs:208
diags:16, lighten_factor:6, obs:176

"""
params={'m':16641,'diags':8,'angle':-1,'lighten_factor':2}
jj = generate_diags(**params)
p = len(jj)


def random_offset(t):
  rstream.seed(int(t/dt*100))
  u = rstream.rand()
  return int(floor(max_offset * u))

def obs_inds(t):
  return jj #+ random_offset(t)

@ens_compatible
def hmod(E,t):
  return E[obs_inds(t)]

from tools.localization import inds_and_coeffs, unravel
xIJ = unravel(arange(m), (ny,nx)) # 2-by-m
def locf(radius,direction,t,tag=None):
  """
  Prepare function:
  inds, coeffs = locf_at(state_or_obs_index)
  """
  yIJ = xIJ[:,obs_inds(t)] # 2-by-p
  def ic(cntr, domain):
    return inds_and_coeffs(cntr, domain, (ny,nx), radius, tag=tag)
  if direction is 'x2y':
    # Could pre-compute ic() for all xIJ,
    # but speed gain is not worth it.
    def locf_at(ind):
      return ic(xIJ[:,ind], yIJ)
  elif direction is 'y2x':
    def locf_at(ind):
      return ic(yIJ[:,ind], xIJ)
  else: raise KeyError
  return locf_at

h = {
    'm'    : p,
    'model': hmod,
    'noise': GaussRV(C=4*eye(p)),
    'loc_f': locf,
    }

setup = TwinSetup(f,h,t,X0)

setup.name = os.path.relpath(__file__,'mods/')

#Create R here

####################
# Suggested tuning
####################
# Reproduce Fig 7 from Sakov and Oke "DEnKF" paper from 2008.
#from mods.QG.sak08 import setup                                 # Expected RMSE_a:
#cfgs += LETKF(N=25,rot=True,infl=1.04,loc_rad=10,taper='Gauss') # 0.6



