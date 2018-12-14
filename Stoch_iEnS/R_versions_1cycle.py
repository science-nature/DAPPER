# One-cycle EnKF PertObs update with/without projecting (full) R.
# Conclusion: Only if R is a*eye(P) will projecting R be equivalent to classic EnKF.
# Conversely: Projecting is quite unnecessary, also with nonlinearity.
#             Does this hold also for (iterative) EnRML?
#

from common import *
seed(2)

M = 5
R = eye(M)
# R = diag(1+arange(M))
# R = randcov(M)

N  = 3
N1 = N-1
E  = randn((M,N))

# Nonlinerity makes no difference to conclusions.
hE = E + 3 # Linear
# hE = sin(E) + 3 # Non-Lin
# fig, ax = plt.subplots()
# lh = ax.plot(E, hE, 'o')

X, xb = center(E ,1)
Y, yb = center(hE,1)

P = Y @ tinv(Y)
y = 4*ones((M,1))
Innov = y - hE - sqrtm(R) @ randn((M,N))

# These become equal when R is a*eye(P):
V1 = Y.T @ ( Y@Y.T + N1*    R    )
V2 = Y.T @ ( Y@Y.T + N1* P@ R @P )
# The (non-)equality of these is then just a corollary:
E1 = E + X@ V1 @ Innov # Classic
E2 = E + X@ V2 @ Innov # Projected

spell_out(E1-E2)

