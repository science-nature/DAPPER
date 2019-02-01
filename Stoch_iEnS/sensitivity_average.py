# Test (with K repetitions) if
#  - Expected projection  matrix == identity. Answer: yes, if N>Nx.
#  - Expected sensitivity matrix == Cyx/Cx.
#    Answer: Only if N>Nx and the obs operator is linear or 2nd-O polynom.
#            See 1d (quicker) tests in sensitivity_expected_cov_1d.
#
# Also recall: Cyx/Cx == E(jacob) if x is Gaussian.
#

## Preamble
from common import *

K  = 10**5   # Num of experiments
Nx = 5       # State length
Ny = 4       # Obs length
N  = 5       # Ens size

c = lambda E: mean0(E,axis=1) # centering

## Random, but repeatable experiment settings
seed(6)

# B = randcov(Nx)              # yields non-diag Pi_infty. Why? 
B = np.diag(1+np.arange(Nx)) # yields diagonal Pi_infty.
B12 = sqrtm(B)

# Linear (part of) obs operator
H = np.round(10*rand((Ny,Nx))) 

# Nonlinear obs operator
h = lambda x: (H@x)**2 + 3*(H@x) + 4                           # A 2nd-O polynomial.
# dhdx = lambda x: 2*diag(H@x)           @H + 3*H              # Works for 1d x, only.
dhdx =   lambda x: 2*    (H@x)[None,:].T *H + 3*H              # Works for 1d x and col-wise E.
# x0 = randn(Nx); dhdx(x0) - approx_jacob(h, x0, colwise=True) # validate dhdx


## One "infinite" experiment (can use K up to 1e7).
seed()

# h_infty == Expect(dhdx(x)).
E_infty = B12@randn((Nx,K))
h_infty = dhdx(E_infty).mean(axis=0)                                            # v1
# h_infty = c(h(E_infty)) @ tinv(c(E_infty))                                    # v2
# h_infty = ( c(h(E_infty)) @ c(E_infty).T ) @ inv(c(E_infty) @ c(E_infty).T)   # v3
# Versions v1 and v2 converge as K-->infty. Of course, v2==v3.


## # Many experiments (loop over K=1e5 takes 1 min)
bPi = zeros((Nx,Nx))
bH  = zeros((Ny,Nx))
bh  = zeros((Ny,Nx))
for k in progbar(range(K)):
  E = B12@randn((Nx,N))
  invE = tinv(c(E), threshold=1-1e-4)

  bPi += c(   E ) @ invE / K
  bH  += c( H@E ) @ invE / K
  bh  += c( h(E)) @ invE / K

##
print("\n=================\nLinear experiment:")
print("****H\n"               ,        H)
print("****H bPi\n"           , round2(H @ bPi           , 0.01))
print("****mean bH\n"         , round2(bH                , 0.01))

print("\n=================\nLinearization experiment:")
print("****h_infty\n"         , round2(h_infty           , 0.01))
print("****h_infty bPi\n"     , round2(h_infty @ bPi     , 0.01))
print("****mean bh\n"         , round2(bh                , 0.01))

print("\n=================\nProj matrix:")
print("****mean bPi\n"        , round2(bPi               , 0.01))
print("****mean trace(bPi)\n" , trace(bPi))

##


