# "Lorenz-63"  model. Classic exhibitor of chaos.
# Phase-plot looks like a butterfly.
# See demo.py for more info.

import numpy as np
from tools.math import with_rk4, is1d, ens_compatible, integrate_TLM

# Constants
sig = 10.0; rho = 28.0; beta = 8.0/3

@ens_compatible
def dxdt(x):
  d     = np.zeros_like(x)
  x,y,z = x
  d[0]  = sig*(y - x)
  d[1]  = rho*x - y - x*z
  d[2]  = x*y - beta*z
  return d

step = with_rk4(dxdt,autonom=True)


def TLM(x):
  """Tangent linear model"""
  assert is1d(x)
  x,y,z = x
  TLM=np.array(
      [[-sig , sig , 0],
      [rho-z , -1  , -x],
      [y     , x   , -beta]])
  return TLM

def dfdx(x,t,dt):
  """Integral of TLM. Jacobian of step."""
  return integrate_TLM(TLM(x),dt,method='approx')


