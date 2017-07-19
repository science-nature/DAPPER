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
      yy[kObs] = h(xx[k],t)
      bruit = h.noise.sample(1)

      if reorder:
        d = int(f.m**0.5)
        space = f.m//(h.m-1) #(hm-1) because starts at 0 and we want the #spaces
        bruit = map(lambda m:(m[1],(m[0]*space)//d+(m[0]*space)%d),enumerate(bruit.squeeze()))
        bruit = sorted(bruit,key=lambda x:x[1])
        bruit = map(lambda x:x[0], bruit)
        bruit = fromiter(bruit,dtype=int)
        bruit = expand_dims(bruit,axis=0)

      yy[kObs] = h(xx[k],t) + bruit

  return xx,yy

def print_together(*args):
  "Print stacked 1D arrays."
  print(np.vstack(args).T)
