# Test graphics/plotting.
# This won't automatically verify if the plots are correct,
# only whether they cause errors or not.

from common import *

import tools.utils as utils
utils.disable_user_interaction = False # NB remember to set to True

sd0 = seed_init(3)

cfgs  = List_of_Configs()

# from mods.LotkaVolterra.dpr01 import HMM   ##################### Expected RMSE_a:
# HMM.t.T = 1000
# cfgs += EnKF_N(N=6,LP=[9])
# # cfgs += ExtKF(infl=1.2,LP=[9])


from mods.Lorenz63.sak12 import HMM
# shorten experiment
HMM.t.BurnIn = 0.1
HMM.t.T = 1
# HMM.t = Chronology(0.01,dtObs=0.24,T=4,BurnIn=0.5,Tplot=4)

# Specify a DA method configuration
cfgs += EnKF('Sqrt',   N=10 ,infl=1.02 ,rot=True, liveplotting=1)
# cfgs += EnKF_N(        N=10            ,rot=True, liveplotting=1)
# cfgs += PartFilt(      N=20 ,reg=2.4   ,NER=0.3 , liveplotting="all")
# cfgs += PFxN(xN=1000,  N=30  ,Qs=2     ,NER=0.2 , liveplotting=[1])
# cfgs += iEnKS('Sqrt',  N=10,  infl=1.02,rot=True, liveplotting=[1])

# Very quick experiment
# HMM.t.BurnIn = HMM.t.dtObs
# HMM.t.KObs = 2

for iC,C in enumerate(cfgs):
  cfgs[iC] = C.update_settings(fail_gently=False,store_u=False)


##############################
# Assimilate
##############################
xx,yy = simulate(HMM)

stats = []
avrgs = []

for ic,config in enumerate(cfgs):
  seed(sd0+2)

  stats += [ config.assimilate(HMM,xx,yy) ]
  avrgs += [ stats[ic].average_in_time() ]
  print_averages(config, avrgs[-1])

print_averages(cfgs,avrgs,statkeys=['rmse_a'])

# plot_time_series(stats[-1], t2=1)
# plot_time_series(stats[-1], t2=0.8, t1=0.2)
plot_time_series(stats[-1], t2=np.inf)
# plot_time_series(stats[-1], t2=0.2)
# plot_time_series(stats[-1], t2=1, t1=0.9)



# print(HMM)
# print(config)
# print(stats)
# print(avrgs)
