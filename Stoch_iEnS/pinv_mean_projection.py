##

# Test if expected
#  - projection matrix is identity  (answer: only if N>=Nx)
#  - sensitivity matrix is unbiased (answer: only if N>=Nx -- not Ny)

## Preamble
from common import *

K  = 10**4   # Num of experiments
Nx = 3       # State length
Ny = 4       # Obs length
N  = 6       # Ens size

cntr = lambda E: mean0(E,axis=1)

seed(6)
#B = randcov(Nx)             # yields non-diag Pi_infty. Why? 
B = np.diag(1+np.arange(Nx)) # yields diagonal Pi_infty.
H = np.round(10*rand((Ny,Nx)))
Obs = lambda x: (H@x)**2 + 3*(H@x) + 4 # 2nd-D polynomial
seed()

# "infty" is exaggerated. There's still noticeable sampling error here.
E_infty = sqrtm(B)@randn((Nx,K))
h_infty = cntr(Obs(E_infty)) @ tinv(cntr(E_infty))


P_av = zeros((Nx,Nx))
H_av = zeros((Ny,Nx))
h_av = zeros((Ny,Nx))
for k in range(K):
  E  = sqrtm(B)@randn((Nx,N))
  iE = tinv(cntr(E), threshold=1-1e-3)

  P_av +=   E  @ iE / K
  H_av += H@E  @ iE / K
  h_av += Obs(E) @ iE / K

##
print("****trace(P_av)\n", trace(P_av))
print("****P_av\n", round2(P_av  , 0.01))

print("****H\n"    , H)
print("****H_av\n" , round2(H_av , 0.01))

print("****h_infty\n" , round2(h_infty , 0.01))
print("****h_av\n"    , round2(h_av    , 0.01))

##


