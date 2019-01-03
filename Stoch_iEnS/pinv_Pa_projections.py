# Explore relationship between (KG-form / ens-subspace-form) P_bar
#     X @ inv(N1*eye(N) + Y.T @ inv(R) @ Y) @ X.T
# and the state-form of P_bar
#     inv( tinv(X @ X.T)/N1 + H.T @ inv(R) @ H )


####################################
# Preamble
####################################
from common import *

sd0 = seed_init(3)

def mean1(E): return mean(E,axis=1,keepdims=True)

Nx = 4       # State length
Ny = 4       # Obs length
N  = 3       # Ens size      # NB must have N1<Nx (want inv(B_e) undefined!)
N1 = N-1

Pi1   = ones((N,N))/N
PiAN  = eye(N) - Pi1

#B = randcov(Nx)                       # non-diag
B = np.diag(1+np.arange(Nx))          # diagonal
R = randcov(Ny)                       # non-diag
H = np.round(10*rand((Ny,Nx)))
#Obs = lambda x: (H@x)**2 + 3*(H@x) + 4  # 2nd-D polynomial
Obs = lambda x: H@x                     # Linear

seed(sd0)


E   = randn((Nx,N))
X   = E - mean1(E)
B_e = X@X.T / N1
Eo  = Obs(E)
Y   = Eo - mean1(Eo)

PiXT = tinv(X) @ X; PiXTC = eye(N)  - PiXT
PiX  = X @ tinv(X); PiXC  = eye(Nx) - PiX


####################################
# KG-form / ens-subspace-form of P_bar (the "right" way)
####################################
#KG = X   @ Y.T @ inv( Y @ Y.T    +    R * N1 )
KG = B_e @ H.T @ inv( H @ B_e @ H.T + R )
#P0 = X @ X.T / N1 - KG@Y@X.T / N1
P0 = X @ X.T / N1 - X@Y.T@inv( Y @ Y.T + R * N1 )@Y@X.T / N1
#P0 = X @ inv(N1*eye(N) + Y.T @ inv(R) @ Y) @ X.T
#P0 = (eye(Nx) - KG @ H) @ B_e
print("\nP0: KG-form / ens-subspace-form of P_bar (the 'right' way)")
spell_out(P0)

print("\nShould be zero:")
spell_out(PiXC @ B_e)
spell_out(PiXC @ P0)


####################################
# state-form of P_bar
####################################

print("\nP2 will differ from P0 because of tinv()...")
#P2 = inv( tinv(B_e)        + H.T @ inv(R) @ H )
P2 = inv( tinv(X @ X.T)*N1 + H.T @ inv(R) @ H )
print("Notably the analysis has **added** uncertainty outside of the ensemble subspace!:")
spell_out(PiXC @ P2 @ PiXC)
print("Projecting into ensemble subspace does not make it equal P0:")
spell_out(PiX @ P2 @ PiX)


####################################
# Try to make a KG-updated P_bar that equals P2
####################################
# Result: fail





