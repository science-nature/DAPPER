# Settings from
# Pajonk, Oliver, et al. 
#   "A deterministic filter for non-Gaussian Bayesian estimation—applications to dynamical system estimation with noisy measurements."
#   Physica D: Nonlinear Phenomena 241.7 (2012): 775-788.
#
# More interesting settings: mods.Lorenz84.harder

from common import *

from mods.Lorenz84.core import step, dfdx
from mods.Lorenz63.liveplotting import LP_setup

Nx = 3
Ny = Nx


day = 0.05/6 * 24 # coz dt=0.05 <--> 6h in "model time scale"
t = Chronology(0.05,dkObs=1,T=200*day,BurnIn=10*day)

Nx = 3
Dyn = {
    'M'    : Nx,
    'model': step,
    'jacob': dfdx,
    'noise': 0
    }

X0  = GaussRV(C=0.01,M=Nx) # Decreased from Pajonk's C=1.

Obs = {
    'M'    : Ny,
    'model': Id_op(),
    'jacob': Id_mat(Ny),
    'noise': 0.1,
    }

other = {'name': os.path.relpath(__file__,'mods/')}

HMM = HiddenMarkovModel(Dyn,Obs,t,X0,**other)

HMM.liveplotting = LP_setup(arange(Nx))

####################
# Suggested tuning
####################
# cfgs += ExtKF(infl=2)
# cfgs += EnKF('Sqrt',N=3,infl=1.01)
# cfgs += PartFilt(reg=1.0, N=100, NER=0.4) # add reg!
# cfgs += PartFilt(reg=1.0, N=1000, NER=0.1) # add reg!
