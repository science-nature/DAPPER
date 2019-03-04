# Test graphics/plotting.
# This won't automatically verify if the plots are correct,
# only whether they cause errors or not.

from common import *

import tools.utils as utils
utils.disable_user_interaction = True # NB remember to set to True

sd0 = seed_init(3)

def test_L63():
  from mods.Lorenz63.sak12 import HMM
  HMM.t.BurnIn = HMM.t.dtObs
  HMM.t.KObs = 2

  cfgs  = List_of_Configs()
  cfgs += EnKF('Sqrt',   N=10 ,infl=1.02 ,rot=True)
  cfgs += PartFilt(      N=20 ,reg=2.4   ,NER=0.3)
  # cfgs += iEnKS('Sqrt',  N=10,  infl=1.02,rot=True)

  for iC,C in enumerate(cfgs):
    C.fail_gently=False
    C.store_u=True
    C.liveplotting="all"

  xx,yy = simulate(HMM)

  stats = []
  avrgs = []

  for ic,config in enumerate(cfgs):
    seed(sd0+2)

    stats += [ config.assimilate(HMM,xx,yy) ]
    avrgs += [ stats[ic].average_in_time() ]
    print_averages(config, avrgs[-1])

  print_averages(cfgs,avrgs,statkeys=['rmse_a'])

  for s in stats:
    replay(s,"all")
  replay(stats[-1], t2=1)
  replay(stats[-1], t2=0.0)
  replay(stats[-1], t2=0.3)
  replay(stats[-1], t2=0.8)
  replay(stats[-1], t2=0.8, t1=0.2)
  replay(stats[-1], t2=np.inf)
  replay(stats[-1], t2=np.inf, speed=1)
  replay(stats[-1], t2=np.inf, pause_a=0, pause_f=0)

  print(HMM); print(config); print(stats); print(avrgs)

  assert True
  return HMM, xx, yy, cfgs, stats, avrgs





# Non py.test runs:
# HMM, xx, yy, cfgs, stats, avrgs = test_L63()
