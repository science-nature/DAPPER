# Explore filter performance as a function of the 
# nu (nX) parameter for EnKF_N.

from common import *

sd0 = seed(4)

##############################
# DA Configurations
##############################
cfgs  = List_of_Configs()

cfgs += Climatology()
# cfgs += OptInterp()
# cfgs += Var3D()

dk_range = arange(1,6)
from mods.Lorenz95.sak08 import HMM
cfgs += EnKF_N(N=20,nX=1.0)
cfgs += EnKF_N(N=20,nX=2)
cfgs += EnKF_N(N=20,nX=1,g=1)
cfgs += EnKF_N(N=20,nX=1,g=2)
cfgs += EnKF_N(N=20,nX='Chi2Fit')
cfgs += EnKF_N(N=20,nX='adapt')
cfgs += EnKF_N(N=20,nX='corr1')

##############################
# Assimilate
##############################
avrgs = np.empty((len(dk_range),len(cfgs)),dict)
stats = np.empty_like(avrgs)

for i,dk in enumerate(dk_range):
  seed(sd0)
  print_c('\ndkObs: ', dk)
  setup.t.dkObs = dk
  setup.t.KObs  = 2000
  xx,yy = simulate(setup)
  for k,config in enumerate(cfgs):
    seed(sd0)
    stats[i,k] = config.assimilate(setup,xx,yy)
    avrgs[i,k] = stats[i,k].average_in_time()
  print_averages(cfgs,avrgs[i])


infl_avrgs = {}
infl_avrgs['mean'] = []
infl_avrgs['var']  = []
infl_avrgs['90']   = []
infl_avrgs['nu']   = []
for k,dk in enumerate(dk_range):
  infls = stats[k,0].infl[setup.t.maskObs_BI]
  infl_avrgs['mean'] += [np.mean(infls)]
  infl_avrgs['var']  += [np.var(infls)]
  infl_avrgs['90']   += [np.percentile(infls,90)]
  i_min_rmse          = np.argmin([a['rmse_a'].val for a in avrgs[k]])
  infl_avrgs['nu']   += [cfgs[i_min_rmse].nu]

infl_avrgs = {key:array(infl_avrgs[key]) for key in infl_avrgs}
