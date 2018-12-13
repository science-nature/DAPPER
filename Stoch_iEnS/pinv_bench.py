# Test if EnKF using pinv(barB) is any good.
# Conclusion: it ain't.
#
# To run script, first insert this snippet in da_methods.py:EnKF_analysis():
#     if 'PertObs' in upd_a:
#         ...
#         E  = E + dE
#     elif 'pinv' in upd_a:
#         # Uses tinv(Bb) to compute Pb
#         D  = mean0(hnoise.sample(N))
#         H  = stats.setup.h.jacob(np.nan, np.nan)
#         Bb = A.T @ A / N1
#         P  = inv( tinv(Bb) + R.inv )
#         if   'KG' in upd_a:
#           KG = P @ H.T @ R.inv
#           E  = E + (KG @ ( y + D - hE ).T).T
#         elif 'Hess' in upd_a:
#           E  = ( P @ ( tinv(Bb)@E.T + R.inv@(y + D).T ) ).T
#     elif 'Sqrt' in upd_a:
#         ...


from common import *

sd0 = seed_init(4)

from mods.Lorenz95.sak08 import HMM
setup.t = Chronology(dt=0.05, dkObs=1, T=20, BurnIn=5)

##############################
# DA Configurations
##############################
cfgs  = List_of_Configs()
cfgs += OptInterp()

for UPD_A in ['PertObs','pinv KG','pinv Hess']:
  for INFL in [1.05, 1.10, 1.20, 1.30, 1.50, 1.70]:
    cfgs += EnKF(UPD_A,N=22,infl=INFL)


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
  print_averages(config, avrgs[-1])
print_averages(cfgs,avrgs,[],['rmse_a','rmv_a'])




