# "Lorenz-63"  model. Classic exhibitor of chaos.
# Phase-plot looks like a butterfly.
# See demo.py for more info.

import numpy as np
from tools.math import with_rk4, is1d, ens_compatible, integrate_TLM

# Constants
sig = 10.0; rho = 28.0; beta = 8.0/3

# Dynamics: time derivative.
@ens_compatible
def dxdt(x):
  d     = np.zeros_like(x)
  x,y,z = x
  d[0]  = sig*(y - x)
  d[1]  = rho*x - y - x*z
  d[2]  = x*y - beta*z
  return d

# Dynamics: time step integration.
step = with_rk4(dxdt,autonom=True)

# Time scale of system. Used as default plotting length.
Tplot = 4.0

# Example initial state.
# Specifics are usually not important coz system is chaotic,
# and we employ a BurnIn before averaging statistics.
# But it's often convenient to give a point on the attractor,
# or at least its basin, or at least ensure that it's "physical".
x0 = np.array([1.509, -1.531, 25.46])


################################################
# OPTIONAL (not necessary for EnKF or PartFilt):
################################################
def TLM(x):
  """Tangent linear model"""
  x,y,z = x
  A = np.array(
      [[-sig , sig , 0],
      [rho-z , -1  , -x],
      [y     , x   , -beta]])
  return A

def dfdx(x,t,dt):
  """Integral of TLM. Jacobian of step."""
  return integrate_TLM(TLM(x),dt,method='approx')


from mods.LotkaVolterra.liveplotting import sliding_marginals
from mods.Lorenz63     .liveplotting import phase3D

def LP(dt,jj=None):
  LP1 = sliding_marginals(jj, lag=int(Tplot/dt), labels='xyz')
  LP2 = phase3D(jj)
  return [LP1, LP2]
