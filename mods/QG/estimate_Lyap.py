# Estimtate the Lyapunov spectrum of the QG model,
# using a limited (rank-N) ensemble.
# Inspired by EmblAUS/Lor95_Lyap.py
from common import *
from mods.QG.core import shape, step, sample_filename, dt
import mods.QG.core as mod

# NB: "Sometimes" multiprocessing does not work here.
# This may be nn ipython bug (stackoverflow.com/a/45720872).
# Solutions: 1) run script from outside of ipython,
#         or 2) Turn it off:
mod.mp = True

sd0 = seed(5)

eps = 0.01 # ensemble rescaling

T  = 1000.0
K  = round(T/dt)
tt = linspace(dt,T,K)

m  = np.prod(shape) # ndim

########################
# Main loop
########################

x = np.load(sample_filename)['sample'][-1]

# Init U
N = 300
U = eye(m)[:N]
E = x + eps*U

LL_exp = zeros((K,N))

for k,t in enumerate(tt):
  if t%10.0==0: print(t)

  x         = step(x,t,dt)
  E         = step(E,t,dt)

  E         = (E-x).T/eps
  [Q, R]    = sla.qr(E,mode='economic')
  E         = x + eps*Q.T
  LL_exp[k] = log(abs(diag(R)))


# Running averages
running_LS = ( tp(1/tt) * np.cumsum(LL_exp,axis=0) )
LS         = running_LS[-1]
print('Lyap spectrum estimate at t=T:')
with printoptions(precision=4): print(LS)
n0 = sum(LS >= 0)
print('n0: ', n0)


#########################
## Plot
#########################
plt.clf()
plt.plot(tt,running_LS,lw=1,alpha=0.4)
plt.title('Lyapunov Exponent estimates')
plt.xlabel('Time')
plt.ylabel('Exponent value')





