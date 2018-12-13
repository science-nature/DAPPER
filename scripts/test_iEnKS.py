# Tests with the LA model, which are very useful for testing
# boundary cases of the iEnKS (e.g. nIter=1, Lag=0).

from common import *

sd0 = seed_init(3)

cfgs  = List_of_Configs(unique=True)

from mods.LA.small import HMM
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


##############################
# Test aux functions
##############################
inds = cfgs.inds # lookup inds

def allsame(xx): 
  return np.allclose(xx, xx[0], atol=1e-10, rtol=0) # rtol=0 => only atol matters.

def gr(ii,f_a_u='a'): # grab_rmses
  return [avrgs[i]['rmse_'+f_a_u].val for i in ii]


##############################
# Tests
##############################
# These don't evaluate the actual numbers, just that they match with the results from other
# methods that are supposed to be equivalent.
# * Understanding why the matches should (and shouldn't) arise is very instructive.
# * Tests can be run with any initial seed.

# For Sqrt flavour:
# - All filter/analysis (f/a) stats should be equal for all EnKF/EnKS/iEnKS(nIter/Lag/MDA).
ii = inds(strict=1,upd_a='Sqrt')
def test_EnKF_Sqrt():                                      assert allsame(gr(ii,'a')) and allsame(gr(ii,'f'))
# - However, the smoother stats (u) depend on the Lag, of course:
ii = inds(strict=0,upd_a='Sqrt',Lag=0) # Also has the non-iter filters
def test_EnKF_Sqrt_u():                                    assert allsame(gr(ii,'u'))
ii = inds(strict=1,upd_a='Sqrt',Lag=1) # Also has the non-iter smoothers
def test_EnKF_Sqrt_Lag1_u():                               assert allsame(gr(ii,'u'))
ii = inds(strict=1,upd_a='Sqrt',Lag=3) # Also has the non-iter smoothers
def test_EnKF_Sqrt_Lag3_u():                               assert allsame(gr(ii,'u'))

# For PertObs flavour::
# - f/a stats all equal except for MDA with nIter>1.
ii = inds(strict=0,upd_a='PertObs',MDA=0) +\
     inds(strict=1,upd_a='PertObs',MDA=1,nIter=1)
def test_EnKF_PertObs():                                   assert allsame(gr(ii,'a')) and allsame(gr(ii,'f'))
# - Still with nIter=4, f/a stats of MDA does not depend on Lag.
ii = inds(strict=1,upd_a='PertObs',nIter=4,MDA=1)
def test_EnKF_PertObs_nIter4():                            assert allsame(gr(ii,'a')) and allsame(gr(ii,'f'))
# - u stats equal for filter and (iter/non-iter) smoothers with Lag=0, except MDA with nIter>1:
ii = inds(strict=0,upd_a='PertObs',MDA=0,Lag=0) +\
     inds(strict=1,upd_a='PertObs',MDA=1,Lag=0,nIter=1)
def test_EnKF_PertObs_u():                                 assert allsame(gr(ii,'u'))
# - u stats equal for            (iter/non-iter) smoothers with Lag=1, except MDA with nIter>1:
ii =      inds(upd_a='PertObs',Lag=1)
ii.remove(inds(upd_a='PertObs',Lag=1,MDA=1,nIter=4)[0])
def test_EnKF_PertObs_Lag1_u():                            assert allsame(gr(ii,'u'))
# - u stats equal for            (iter/non-iter) smoothers with Lag=3, except MDA with nIter>1:
ii =      inds(upd_a='PertObs',Lag=3)
ii.remove(inds(upd_a='PertObs',Lag=3,MDA=1,nIter=4)[0])
def test_EnKF_PertObs_Lag3_u():                            assert allsame(gr(ii,'u'))

# For Order1 (DEnKF) flavour:
# - f/a stats all equal except for nIter>1:
ii = inds(upd_a='DEnKF') +\
     inds(upd_a='Order1',nIter=1)
def test_EnKF_Order1():                                    assert allsame(gr(ii,'a')) and allsame(gr(ii,'f'))
# - f/a stats independent of Lag for non-MDA and a given nIter:
ii = inds(upd_a='Order1',nIter=4,MDA=0)
def test_EnKF_Order1_nIter4_MDA0():                        assert allsame(gr(ii,'a')) and allsame(gr(ii,'f'))
# - f/a stats independent of Lag for     MDA and a given nIter:
ii = inds(upd_a='Order1',nIter=4,MDA=1)
def test_EnKF_Order1_nIter4_MDA1():                        assert allsame(gr(ii,'a')) and allsame(gr(ii,'f'))
# - u   stats equal for EnKS/iEnKS(nIter=1) for a given Lag:
ii = inds(strict=0,upd_a='DEnKF' ,Lag=0) +\
     inds(         upd_a='Order1',Lag=0,nIter=1)
def test_EnKF_Order1_nIter1_Lag0_u():                      assert allsame(gr(ii,'u'))
ii = inds(da=EnKS, upd_a='DEnKF' ,Lag=1) +\
     inds(         upd_a='Order1',Lag=1,nIter=1)
def test_EnKF_Order1_nIter1_Lag1_u():                      assert allsame(gr(ii,'u'))
ii = inds(da=EnKS, upd_a='DEnKF' ,Lag=3) +\
     inds(         upd_a='Order1',Lag=3,nIter=1)
def test_EnKF_Order1_nIter1_Lag3_u():                      assert allsame(gr(ii,'u'))

# A table perspective of the tests:
# Each benchmark result has been replaced by a letter,
# but the letter is repeated for equal results.
# I.e. the tests check that the letters occur everywhere they should.
# The 'X#' stand for results that do not repeat anywhere.

#       da_method  upd_a    Lag  MDA    nIter  |  rmse_a  rmse_f   rmse_u
# ----  ---------  -------  ---  -----  -----  -  ------  ------  --------
# [0]   EnKF       Sqrt                        |       A       B        C
# [1]   EnKF       PertObs                     |       F       G        J
# [2]   EnKF       DEnKF                       |       M       N        S


# [3]   EnKS       Sqrt       0                |       A       B        C
# [4]   EnKS       PertObs    0                |       F       G        J
# [5]   EnKS       DEnKF      0                |       M       N        S

# [6]   iEnKS      Sqrt       0  False      1  |       A       B        C
# [7]   iEnKS      PertObs    0  False      1  |       F       G        J
# [8]   iEnKS      Order1     0  False      1  |       M       N        S
# [9]   iEnKS      Sqrt       0  True       1  |       A       B        C
# [10]  iEnKS      PertObs    0  True       1  |       F       G        J
# [11]  iEnKS      Order1     0  True       1  |       M       N        S

# [12]  iEnKS      Sqrt       0  False      4  |       A       B        C
# [13]  iEnKS      PertObs    0  False      4  |       F       G        J
# [14]  iEnKS      Order1     0  False      4  |       O       P       X4
# [15]  iEnKS      Sqrt       0  True       4  |       A       B        C
# [16]  iEnKS      PertObs    0  True       4  |       H       I       X1
# [17]  iEnKS      Order1     0  True       4  |       Q       R       X5


# [18]  EnKS       Sqrt       1                |       A       B        D
# [19]  EnKS       PertObs    1                |       F       G        K
# [20]  EnKS       DEnKF      1                |       M       N        T

# [21]  iEnKS      Sqrt       1  False      1  |       A       B        D
# [22]  iEnKS      PertObs    1  False      1  |       F       G        K
# [23]  iEnKS      Order1     1  False      1  |       M       N        T
# [24]  iEnKS      Sqrt       1  True       1  |       A       B        D
# [25]  iEnKS      PertObs    1  True       1  |       F       G        K
# [26]  iEnKS      Order1     1  True       1  |       M       N        T

# [27]  iEnKS      Sqrt       1  False      4  |       A       B        D
# [28]  iEnKS      PertObs    1  False      4  |       F       G        K
# [29]  iEnKS      Order1     1  False      4  |       O       P       X6
# [30]  iEnKS      Sqrt       1  True       4  |       A       B        D
# [31]  iEnKS      PertObs    1  True       4  |       H       I       X2
# [32]  iEnKS      Order1     1  True       4  |       Q       R       X7


# [33]  EnKS       Sqrt       3                |       A       B        E
# [34]  EnKS       PertObs    3                |       F       G        L
# [35]  EnKS       DEnKF      3                |       M       N        U

# [36]  iEnKS      Sqrt       3  False      1  |       A       B        E
# [37]  iEnKS      PertObs    3  False      1  |       F       G        L
# [38]  iEnKS      Order1     3  False      1  |       M       N        U
# [39]  iEnKS      Sqrt       3  True       1  |       A       B        E
# [40]  iEnKS      PertObs    3  True       1  |       F       G        L
# [41]  iEnKS      Order1     3  True       1  |       M       N        U

# [42]  iEnKS      Sqrt       3  False      4  |       A       B        E
# [43]  iEnKS      PertObs    3  False      4  |       F       G        L
# [44]  iEnKS      Order1     3  False      4  |       O       P       X8
# [45]  iEnKS      Sqrt       3  True       4  |       A       B        E
# [46]  iEnKS      PertObs    3  True       4  |       H       I       X3
# [47]  iEnKS      Order1     3  True       4  |       Q       R       X9


##############################
# Seed-dependent test
##############################
# Just stupidly compare the full table.
# ==> Test cannot be run with different seeds or computers.
# Seed used: sd0 = seed_init(3). Test run on my Mac.

def test_a():
  assert gr(arange(len(cfgs)),'a') == \
      [0.11696645982504372 , 0.13707255797444515 , 0.16013057281282431 ,
       0.11696645982504382 , 0.13707255797444512 , 0.16013057281282431 ,
       0.11696645982504426 , 0.13707255797444598 , 0.16013057281282492 ,
       0.11696645982504426 , 0.13707255797444598 , 0.16013057281282492 ,
       0.1169664598250439  , 0.13707255797444512 , 0.1186316809275947  ,
       0.11696645982504399 , 0.15353644282849052 , 0.11857822272148509 ,
       0.11696645982504382 , 0.13707255797444512 , 0.16013057281282431 ,
       0.11696645982504422 , 0.13707255797444584 , 0.16013057281282467 ,
       0.11696645982504422 , 0.13707255797444584 , 0.16013057281282467 ,
       0.1169664598250439  , 0.13707255797444509 , 0.1186316809275947  ,
       0.116966459825044   , 0.15353644282849072 , 0.11857822272148516 ,
       0.11696645982504382 , 0.13707255797444512 , 0.16013057281282431 ,
       0.11696645982504425 , 0.13707255797444587 , 0.16013057281282467 ,
       0.11696645982504425 , 0.13707255797444587 , 0.16013057281282467 ,
       0.11696645982504389 , 0.13707255797444509 , 0.11863168092759468 ,
       0.11696645982504399 , 0.15353644282849077 , 0.11857822272148516 ]

def test_f():
  assert gr(arange(len(cfgs)),'f') == \
      [0.19684477973528794 , 0.21429229233638106 , 0.23934387451007286 ,
       0.19684477973528802 , 0.214292292336381   , 0.23934387451007291 ,
       0.19684477973528852 , 0.21429229233638189 , 0.23934387451007347 ,
       0.19684477973528852 , 0.21429229233638189 , 0.23934387451007347 ,
       0.19684477973528816 , 0.21429229233638106 , 0.19828766346772814 ,
       0.19684477973528819 , 0.23049048795386004 , 0.19858885018530187 ,
       0.19684477973528802 , 0.214292292336381   , 0.23934387451007286 ,
       0.19684477973528844 , 0.21429229233638178 , 0.23934387451007322 ,
       0.19684477973528844 , 0.21429229233638178 , 0.23934387451007322 ,
       0.19684477973528811 , 0.214292292336381   , 0.19828766346772814 ,
       0.19684477973528824 , 0.23049048795386023 , 0.19858885018530195 ,
       0.19684477973528802 , 0.214292292336381   , 0.23934387451007286 ,
       0.19684477973528847 , 0.21429229233638178 , 0.23934387451007322 ,
       0.19684477973528847 , 0.21429229233638178 , 0.23934387451007322 ,
       0.19684477973528811 , 0.214292292336381   , 0.19828766346772814 ,
       0.19684477973528819 , 0.23049048795386026 , 0.19858885018530195 ]

def test_u():
  assert gr(arange(len(cfgs)),'u') == \
      [0.18086911575323916  , 0.19884834546399385  , 0.22350121417062313  ,
       0.18086911575323922  , 0.19884834546399388  , 0.22350121417062319  ,
       0.18086911575323966  , 0.19884834546399466  , 0.22350121417062382  ,
       0.18086911575323966  , 0.19884834546399466  , 0.22350121417062382  ,
       0.18086911575323933  , 0.1988483454639938   , 0.18235646695970145  ,
       0.18086911575323938  , 0.21509967892878615  , 0.18258672469253867  ,
       0.10809005405428855  , 0.12872786931335145  , 0.15138717068466812  ,
       0.10809005405428891  , 0.12872786931335212  , 0.15138717068466845  ,
       0.10809005405428891  , 0.12872786931335212  , 0.15138717068466845  ,
       0.10809005405428862  , 0.12872786931335142  , 0.10979974263086155  ,
       0.10809005405428863  , 0.14479974441156154  , 0.10969822797890079  ,
       0.05884180343829841  , 0.087850672189102169 , 0.083764782614186017 ,
       0.058841803438298215 , 0.08785067218910271  , 0.083764782614185893 ,
       0.058841803438298215 , 0.08785067218910271  , 0.083764782614185893 ,
       0.058841803438298222 , 0.087850672189102114 , 0.061641461112443896 ,
       0.058841803438298201 , 0.096727617207417904 , 0.059231682718777934 ]




# For nonlinear dynamics, the (non-iterative) EnKF (f/a/u stats)
# are reproduced by the iEnKS with Lag=0 (and nIter==1 if h is also nonlin).
# However, the 'u' stats of the non-iterative EnKS(Lag>0) are not reproduced.
# Re-use cfgs and test with:
# from mods.Lorenz95.sak08 import HMM
# setup.t.KObs=100 # Here, must use >100 to avoid indistinguishable rmse stats.




