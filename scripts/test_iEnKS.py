from common import *

sd0 = seed_init(5)

cfgs  = List_of_Configs(unique=True)

from mods.LA.small import setup
setup.t.BurnIn=0
setup.t.KObs=10

cfgs +=  EnKF('Sqrt'        , N=20,                      )
cfgs +=  EnKF('PertObs'     , N=20,                      )
cfgs +=  EnKF('DEnKF'       , N=20,                      )
for Lag in [0,1,3]:
  cfgs +=  EnKS('Sqrt'      , N=20, Lag=Lag,             )
  cfgs +=  EnKS('PertObs'   , N=20, Lag=Lag,             )
  cfgs +=  EnKS('DEnKF'     , N=20, Lag=Lag,             )
  for nIter in [1,4]:
    for MDA in [False,True]:
      cfgs += iEnKS('Sqrt'    , N=20, Lag=Lag, nIter=nIter, MDA=MDA)
      cfgs += iEnKS('PertObs' , N=20, Lag=Lag, nIter=nIter, MDA=MDA)
      cfgs += iEnKS('Order1'  , N=20, Lag=Lag, nIter=nIter, MDA=MDA)


for iC,C in enumerate(cfgs):
  cfgs[iC] = C.update_settings(
      liveplotting=False,fail_gently=True,store_u=True)


##############################
# Assimilate
##############################
xx,yy = simulate(setup)

stats = []
avrgs = []

for ic,config in enumerate(cfgs):
  seed(sd0+2)

  stats += [ config.assimilate(setup,xx,yy) ]
  avrgs += [ stats[ic].average_in_time() ]

print_averages(cfgs,avrgs,statkeys=['rmse_a','rmse_f','rmse_u'])


# Example that can be used (after changes) with diff to unit-test developments.
# >>> sd0 = seed_init(5)
# >>> 
# >>> cfgs  = List_of_Configs(unique=True)
# >>> 
# >>> from mods.LA.small import setup
# >>> setup.t.BurnIn=0
# >>> setup.t.KObs=10
# >>> 
# >>> cfgs +=  EnKF('Sqrt'        , N=20,                      )
# >>> cfgs +=  EnKF('PertObs'     , N=20,                      )
# >>> cfgs +=  EnKF('DEnKF'       , N=20,                      )
# >>> for Lag in [0,1,3]:
# >>>   cfgs +=  EnKS('Sqrt'      , N=20, Lag=Lag,             )
# >>>   cfgs +=  EnKS('PertObs'   , N=20, Lag=Lag,             )
# >>>   cfgs +=  EnKS('DEnKF'     , N=20, Lag=Lag,             )
# >>>   for nIter in [1,4]:
# >>>     for MDA in [False,True]:
# >>>       cfgs += iEnKS('Sqrt'    , N=20, Lag=Lag, nIter=nIter, MDA=MDA)
# >>>       cfgs += iEnKS('PertObs' , N=20, Lag=Lag, nIter=nIter, MDA=MDA)
# >>>       cfgs += iEnKS('Order1'  , N=20, Lag=Lag, nIter=nIter, MDA=MDA)
# 
# Yielded output:
#       da_method  upd_a    Lag  MDA    nIter  |  rmse_a ±     rmse_f ±      rmse_u ±
# ----  ---------  -------  ---  -----  -----  -  -----------  -----------  -------------
# [0]   EnKF       Sqrt                        |  0.1294 0.06   0.174 0.07   0.1651 0.09
# [1]   EnKF       PertObs                     |  0.1375 0.07  0.1814 0.08   0.1726 0.09
# [2]   EnKF       DEnKF                       |  0.1614 0.08  0.2049 0.09   0.1962 0.1
# [3]   EnKS       Sqrt       0                |  0.1294 0.06   0.174 0.07   0.1651 0.09
# [4]   EnKS       PertObs    0                |  0.1375 0.07  0.1814 0.08   0.1726 0.09
# [5]   EnKS       DEnKF      0                |  0.1614 0.08  0.2049 0.09   0.1962 0.1
# [6]   iEnKS      Sqrt       0  False      1  |  0.1294 0.06   0.174 0.07   0.1651 0.09
# [7]   iEnKS      PertObs    0  False      1  |  0.1375 0.07  0.1814 0.08   0.1726 0.09
# [8]   iEnKS      Order1     0  False      1  |  0.1614 0.08  0.2049 0.09   0.1962 0.1
# [9]   iEnKS      Sqrt       0  True       1  |  0.1294 0.06   0.174 0.07   0.1651 0.09
# [10]  iEnKS      PertObs    0  True       1  |  0.1375 0.07  0.1814 0.08   0.1726 0.09
# [11]  iEnKS      Order1     0  True       1  |  0.1614 0.08  0.2049 0.09   0.1962 0.1
# [12]  iEnKS      Sqrt       0  False      4  |  0.1294 0.06   0.174 0.07   0.1651 0.09
# [13]  iEnKS      PertObs    0  False      4  |  0.1375 0.07  0.1814 0.08   0.1726 0.09
# [14]  iEnKS      Order1     0  False      4  |  0.1309 0.06  0.1749 0.06   0.1661 0.09
# [15]  iEnKS      Sqrt       0  True       4  |  0.1294 0.06   0.174 0.07   0.1651 0.09
# [16]  iEnKS      PertObs    0  True       4  |   0.188 0.04  0.2269 0.07   0.2191 0.08
# [17]  iEnKS      Order1     0  True       4  |  0.1318 0.07  0.1764 0.08   0.1674 0.09
# [18]  EnKS       Sqrt       1                |  0.1294 0.06   0.174 0.07   0.1213 0.04
# [19]  EnKS       PertObs    1                |  0.1375 0.07  0.1814 0.08   0.1295 0.04
# [20]  EnKS       DEnKF      1                |  0.1614 0.08  0.2049 0.09   0.1534 0.06
# [21]  iEnKS      Sqrt       1  False      1  |  0.1294 0.06   0.174 0.07   0.1213 0.04
# [22]  iEnKS      PertObs    1  False      1  |  0.1375 0.07  0.1814 0.08   0.1295 0.04
# [23]  iEnKS      Order1     1  False      1  |  0.1614 0.08  0.2049 0.09   0.1534 0.06
# [24]  iEnKS      Sqrt       1  True       1  |  0.1294 0.06   0.174 0.07   0.1213 0.04
# [25]  iEnKS      PertObs    1  True       1  |  0.1375 0.07  0.1814 0.08   0.1295 0.04
# [26]  iEnKS      Order1     1  True       1  |  0.1614 0.08  0.2049 0.09   0.1534 0.06
# [27]  iEnKS      Sqrt       1  False      4  |  0.1294 0.06   0.174 0.07   0.1213 0.04
# [28]  iEnKS      PertObs    1  False      4  |  0.1375 0.07  0.1814 0.08   0.1295 0.04
# [29]  iEnKS      Order1     1  False      4  |  0.1309 0.06  0.1749 0.06   0.1229 0.04
# [30]  iEnKS      Sqrt       1  True       4  |  0.1294 0.06   0.174 0.07   0.1213 0.04
# [31]  iEnKS      PertObs    1  True       4  |   0.188 0.04  0.2269 0.07   0.1807 0.04
# [32]  iEnKS      Order1     1  True       4  |  0.1318 0.07  0.1764 0.08   0.1235 0.04
# [33]  EnKS       Sqrt       3                |  0.1294 0.06   0.174 0.07   0.0727 0.01
# [34]  EnKS       PertObs    3                |  0.1375 0.07  0.1814 0.08  0.08216 0.01
# [35]  EnKS       DEnKF      3                |  0.1614 0.08  0.2049 0.09  0.09731 0.02
# [36]  iEnKS      Sqrt       3  False      1  |  0.1294 0.06   0.174 0.07   0.0727 0.01
# [37]  iEnKS      PertObs    3  False      1  |  0.1375 0.07  0.1814 0.08  0.08216 0.01
# [38]  iEnKS      Order1     3  False      1  |  0.1614 0.08  0.2049 0.09  0.09731 0.02
# [39]  iEnKS      Sqrt       3  True       1  |  0.1294 0.06   0.174 0.07   0.0727 0.01
# [40]  iEnKS      PertObs    3  True       1  |  0.1375 0.07  0.1814 0.08  0.08216 0.01
# [41]  iEnKS      Order1     3  True       1  |  0.1614 0.08  0.2049 0.09  0.09731 0.02
# [42]  iEnKS      Sqrt       3  False      4  |  0.1294 0.06   0.174 0.07   0.0727 0.01
# [43]  iEnKS      PertObs    3  False      4  |  0.1375 0.07  0.1814 0.08  0.08216 0.01
# [44]  iEnKS      Order1     3  False      4  |  0.1309 0.06  0.1749 0.06  0.07585 0.009
# [45]  iEnKS      Sqrt       3  True       4  |  0.1294 0.06   0.174 0.07   0.0727 0.01
# [46]  iEnKS      PertObs    3  True       4  |   0.188 0.04  0.2269 0.07   0.1378 0.01
# [47]  iEnKS      Order1     3  True       4  |  0.1318 0.07  0.1764 0.08  0.07394 0.01

# For proper unit-testing, this should of course be made seed-independent.
# This requires testing (instead of the absolute values, as with the above)
# which values match, and that these values also match after the changes you've made. 

# The values that match should be interpreted. For example, in the linear case (LA):
#  - Non-iterative methods yield equal f/a stats as iterative methods,
#    and this is independent of the Lag used.
#    An exception to the above rule is that nIter affects the f/a stats for
#    * PertObs with         MDA
#    * Order1   with/without MDA
#    though they remain equal to non-iterative methods for nIter==1.
# NB: For iEnKS PertObs (not MDA) to be independent of Lag and nIter,
#     and reproduce the non-iterative filter, one must ensure that
#     the perturbations are drawn in the same order, and use the same sign.
#  - EnKS and iEnKS yield the same 'u' stats.
#  - If Lag==0, then this also equals the 'u' stats of EnKF (excpet with Order1).
#
# For nonlinear dynamics, the (non-iterative) EnKF (f/a/u stats)
# are reproduced by the iEnKS with Lag=0 (and nIter==1 if h is also nonlin).
# However, the 'u' stats of the non-iterative EnKS(Lag>0) are not reproduced.
# Re-use cfgs and test with:
# from mods.Lorenz95.sak08 import setup
# setup.t.KObs=100 # Here, must use >100 to avoid indistinguishable rmse stats.




