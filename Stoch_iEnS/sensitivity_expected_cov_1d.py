# Question: does the LS linear regression estimate,
#           barCyx/barCx, have expectation Cyx/Cx ?
# Note: this question does not even pre-suppose
#       the existence of an operator (such as "h") linking y and x,
#       but it's easier to test that way.
# Answer: Only if the "h" is linear or 2nd-O polynom.

##
from common import *
# sd0 = seed_init(3)

def c(xx,axis=None): # centering
  return xx - xx.mean(axis,keepdims=True)

K = 10**4 # num of repetitions of experiment.
N = 3     # size of sample that is used for estimator

xx = randn((N,K))
h  = lambda x: 0*x**3 + x**2 + 3*x + 4
yy = h(xx)

## True covs
print("*************")
Cyx  = ( c(yy)*c(xx) ).sum()       / (N*K-1)
Cx   = ( c(xx)*c(xx) ).sum()       / (N*K-1)
spell_out(Cyx / Cx)


## Mean of estimated covs
print("*************")
barCyx = ( c(yy)*c(xx) ).sum(axis=0) / (N  -1)
barCx  = ( c(xx)*c(xx) ).sum(axis=0) / (N  -1)
spell_out(mean(barCyx / barCx))

##

