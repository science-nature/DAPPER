from common import *

def simulate(setup,desc='Truth & Obs',reorder=False):
  """Generate synthetic truth and observations."""
  f,h,chrono,X0 = setup.f, setup.h, setup.t, setup.X0

  # Init
  xx    = zeros((chrono.K+1,f.m))
  xx[0] = X0.sample(1)
  yy    = zeros((chrono.KObs+1,h.m))

  # Loop
  for k,kObs,t,dt in progbar(chrono.forecast_range,desc):
    xx[k] = f(xx[k-1],t-dt,dt) + sqrt(dt)*f.noise.sample(1)
    if kObs is not None:
      bruit = h.noise.sample(1)

      if reorder:
        from mods.QG.sak08 import jj,params
        angle = params['angle']
        out = array([])
        if not hasattr(angle,'__iter__'):
          angle = [angle]
          bruit = [bruit.squeeze()]
        else:
          bruit = [bruit.squeeze()[:h.m//2],bruit.squeeze()[h.m//2:]]
        for (a,b) in zip(angle,bruit):
          d = int(f.m**0.5)
          b = sorted(zip(jj,b),key=lambda x: x[0]%d - a * (x[0]//d))
          b = map(lambda x: x[1], b)
          b = fromiter(b,dtype=float)

        out = hstack((out,b))
        bruit = expand_dims(out,axis = 0)

      yy[kObs] = h(xx[k],t) + bruit

  return xx,yy

def print_together(*args):
  "Print stacked 1D arrays."
  print(np.vstack(args).T)
