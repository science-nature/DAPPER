# Just stupidly compare the full table.
# Seed-dependent test.
# ==> Test cannot be run with different seeds or computers.

from common import *

from da_methods.admin import _print_averages
def pa(cfgs,avrgs): return _print_averages(cfgs,avrgs,statkeys=['rmse_a','rmse_f','rmse_u'])


##############################
# L63
##############################
from mods.Lorenz63.sak12 import HMM
HMM.t.BurnIn=0
HMM.t.KObs=10
sd0 = seed(9)

# Cfgs
cfgs  = List_of_Configs()
cfgs += Climatology()
cfgs += OptInterp()
cfgs += Var3D(infl=0.9)
cfgs += ExtKF(infl=90)
cfgs += EnKF('Sqrt',    N=3 ,  infl=1.30)
cfgs += EnKF('Sqrt',    N=10,  infl=1.02,rot=True)
cfgs += EnKF('PertObs', N=500, infl=0.95,rot=False)
cfgs += EnKF_N(         N=10,            rot=True)
cfgs += iEnKS('Sqrt',   N=10,  infl=1.02,rot=True)
cfgs += PartFilt(       N=100 ,reg=2.4  ,NER=0.3)
cfgs += PartFilt(       N=800 ,reg=0.9  ,NER=0.2)
cfgs += PartFilt(       N=4000,reg=0.7  ,NER=0.05)
cfgs += PFxN(xN=1000,   N=30  ,Qs=2     ,NER=0.2)

# Run
xx,yy = simulate(HMM)
avrgs = []
for ic,config in enumerate(cfgs):
  config.store_u = True
  seed(sd0+2)
  stats = config.assimilate(HMM,xx,yy)
  avrgs += [ stats.average_in_time() ]

table = pa(cfgs,avrgs)
old = """
      da_method       N  upd_a     infl  rot     NER  Qs  reg    xN  |  rmse_a ±    rmse_f ±    rmse_u ±
----  -----------  ----  -------  -----  -----  ----  --  ---  ----  -  ----------  ----------  ----------
[0]   Climatology                                                    |   8.265 1     8.265 1      7.93 3
[1]   OptInterp                                                      |   1.265 0.2   8.265 1     1.325 0.1
[2]   Var3D                        0.9                               |   1.212 0.2   2.569 1     1.584 0.3
[3]   ExtKF                       90                                 |  0.9001 0.2    1.86 0.8   1.152 0.2
[4]   EnKF            3  Sqrt      1.3                               |  0.8368 0.2   1.514 0.5  0.9236 0.1
[5]   EnKF           10  Sqrt      1.02  True                        |  0.7257 0.3     1.3 0.3  0.8452 0.2
[6]   EnKF          500  PertObs   0.95  False                       |   0.713 0.3   1.289 0.3  0.8216 0.2
[7]   EnKF_N         10                  True                        |  0.7471 0.2   1.317 0.3  0.8661 0.2
[8]   iEnKS          10  Sqrt      1.02  True                        |  0.8153 0.3  0.9145 0.3  0.4997 0.3
[9]   PartFilt      100                         0.3       2.4        |  0.6182 0.2   1.401 0.7  0.8489 0.1
[10]  PartFilt      800                         0.2       0.9        |  0.4684 0.2  0.9491 0.4   0.587 0.1
[11]  PartFilt     4000                         0.05      0.7        |  0.4317 0.1   0.815 0.3  0.5209 0.1
[12]  PFxN           30                         0.2    2       1000  |   0.862 0.2   1.809 0.5   1.155 0.2
"""[1:-1]

def test_len():
  assert len(old)==len(table)

table = table.split('\n')
old   = old  .split('\n')

import pytest
@pytest.mark.parametrize(('lineno'),arange(len(table)))
def test_tables(lineno):
    assert table[lineno] == old[lineno]



##############################
# L95
##############################
from mods.Lorenz95.sak08 import HMM
HMM.t.BurnIn=0
HMM.t.KObs=10
sd0 = seed(9)

# Cfgs
cfgs  = List_of_Configs()
cfgs += Climatology()
cfgs += OptInterp()
cfgs += Var3D(infl=1.05)
cfgs += ExtKF(infl=6)
cfgs += EnKF('PertObs'        ,N=40,infl=1.06)
cfgs += EnKF('Sqrt'           ,N=28,infl=1.02,rot=True)

cfgs += EnKF_N(N=24,rot=True)
cfgs += EnKF_N(N=24,rot=True,xN=2)
cfgs += iEnKS('Sqrt',N=40,infl=1.01,rot=True)

cfgs += LETKF(         N=7,rot=True,infl=1.04,loc_rad=4)
cfgs += SL_EAKF(       N=7,rot=True,infl=1.07,loc_rad=6)


# Run
xx,yy = simulate(HMM)
avrgs = []
for ic,config in enumerate(cfgs):
  config.store_u = True
  seed(sd0+2)
  stats = config.assimilate(HMM,xx,yy)
  avrgs += [ stats.average_in_time() ]

table = pa(cfgs,avrgs)
old = """
      da_method     N  upd_a    infl  rot   loc_rad  xN  |   rmse_a ±        rmse_f ±        rmse_u ±
----  -----------  --  -------  ----  ----  -------  --  -  --------------  --------------  --------------
[0]   Climatology                                        |   0.8347 0.2      0.8347 0.2      0.8347 0.2
[1]   OptInterp                                          |   0.1275 0.03     0.8347 0.2      0.1275 0.03
[2]   Var3D                     1.05                     |  0.08993 0.04     0.2405 0.2     0.08993 0.04
[3]   ExtKF                     6                        |  0.02986 0.0008  0.02987 0.0009  0.02986 0.0008
[4]   EnKF         40  PertObs  1.06                     |  0.03031 0.001   0.03018 0.0009  0.03031 0.001
[5]   EnKF         28  Sqrt     1.02  True               |     0.03 0.001   0.02991 0.001      0.03 0.001
[6]   EnKF_N       24                 True               |   0.0304 0.001   0.03032 0.001    0.0304 0.001
[7]   EnKF_N       24                 True            2  |  0.03039 0.001   0.03031 0.001   0.03039 0.001
[8]   iEnKS        40  Sqrt     1.01  True               |  0.03017 0.001   0.03011 0.0009  0.03017 0.001
[9]   LETKF         7           1.04  True        4      |  0.03027 0.001   0.03022 0.001   0.03027 0.001
[10]  SL_EAKF       7           1.07  True        6      |  0.03067 0.001    0.0305 0.001   0.03067 0.001
"""[1:-1]

table = table.split('\n')
old   = old  .split('\n')

def test_len():
  assert len(old)==len(table)

import pytest
@pytest.mark.parametrize(('lineno'),arange(len(table)))
def test_tables(lineno):
    assert table[lineno] == old[lineno]


