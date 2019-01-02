# Moderate dtObs and non-0 Q.
# Named 'm30' in Datum.


from common import *

from mods.Lorenz63.sak12 import *

t.dkObs = 15
Dyn['noise'] = 2
X0.C = CovMat(0.5*ones(M))

other = {'name': os.path.relpath(__file__,'mods/')}
HMM = HiddenMarkovModel(Dyn,Obs,t,X0,**other)


####################
# Suggested tuning
####################
