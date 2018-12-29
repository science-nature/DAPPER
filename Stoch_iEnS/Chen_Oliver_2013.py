##

# Consider the univariate problem from Chen, Yan and Oliver, Dean S.
# "Levenberg-Marquardt forms of the iES for efficient history matching and UQ"
# alias: CO13
#
# This is here generalized to M dimensions
# by repeating the problem independently for each dim.
# However, sampling error will couple the dimensions somewhat.
#
# Setting M = P = 1 can be used to reproduce the results of the paper.
#
# This script has been used extensively to test the equivalence of various
# forms of the matrix Y_k = H_k X_0.

## Preamble ===================
from common import *
from scipy.stats import norm

#plt.xkcd()

sd0 = seed_init(3)

def mean1(E):
  "Enables correct broadcasting for MxN shaped ensemble matrices"
  return mean(E,axis=1,keepdims=True)

M  = 2       # State length
P  = 2       # Obs length
N  = 4       # Ens size
N1 = N-1     #
nbins = 150  # Histogram bins


## Prior ===================
b  = -2*ones((M,1))
B  = 1*eye(M)
E0 = b + sqrtm(B)@randn((M,N))


## Obs  ===================
jj = arange(P)
y  = 48*ones((P,1))
R  = 16*eye(P)
#R  = 1*eye(P)
if True:
  # == Non-Lin H == (as in CO13):
  def h1(x): return 7/12*x*x*x - 7/2*x*x + 8*x # A scalar fun...
  def  Obs(x): return h1(x) [jj,:]               # ...duplicated to obs dims.
  def hp(x):                                   # The tangent.
    H = zeros((P,M))
    for i,j in enumerate(jj):
      H[i,j] = 7/4*x[j]**2 - 7*x[j] + 8
    return H
else:
  # == Linear H == (not in CO13):
  def h1(x): return 5*x
  def  Obs(x): return h1(x) [jj,:]
  def hp(x): return 5*ones((P,M))


## PDFs ===================
xx = linspace(-5,10,2001)
dx = xx[1]-xx[0]
def normlz(pp):
  return pp / sum(pp) / dx

prior_xx = norm.pdf(xx,b[0],sqrt(B[0,0]))
lklhd_xx = normlz( norm.pdf(y[0],h1(xx),sqrt(R[0,0])) )
postr_xx = normlz( prior_xx * lklhd_xx )


## Plotting PDFs ==========
fig = plt.figure(1)
fig.clear()
fig, (ax1, ax2) = plt.subplots(2,1,sharex=True,num=1, gridspec_kw = {'height_ratios':[5, 1]})

ax1.hist(E0[0,:],bins=linspace(*xx[[0,-1]],nbins),
    normed=True, label='$E_0$',alpha=0.6)

ax1.plot(xx,prior_xx, 'b-' , label='$p(x)$')
ax1.plot(xx,lklhd_xx, 'g-' , label='$p(y='+str(y[0,0])+'|x)$')
ax1.plot(xx,postr_xx, 'r--', label='$p(x|y)$')

with np.errstate(divide='ignore'):
  ax2.plot(xx,-log( prior_xx ), 'b-' )
  ax2.plot(xx,-log( lklhd_xx ), 'g-' )
  ax2.plot(xx,-log( postr_xx ), 'r--')
  ax2.set_ylim(-2,70)

ax1.set_ylabel('p')
ax2.set_ylabel('-log p')
ax1.legend()
ax2.set_xlabel('x')


## Ensemble setup ===========

# Initialize ensemble
E   = E0
x0  = mean1(E0)
A0  = E0 - x0
B0  = A0@A0.T/N1

# Proj matrices
Pi0  = tinv(A0)@A0
Pi0C = eye(N) - Pi0
Pi1  = ones((N,N))/N
AN   = eye(N) - Pi1

# Obs perturbations
D = sqrtm(R)@randn((P,N))
D = mean0(D.T).T

# Sqrt initialization matrix
w      = zeros((N,1))
W      = eye(N)
T      = eye(N)
# Geir-Evensen formulation
We     = W - eye(N)
Om     = eye(N) + We@AN # = Pi1 + W@AN


## Algorithm  ===============
#FORM='MDA'                   # Ensemble Multiple Data Assimilation
#FORM='RML-GN'                # RML with exact, local gradients. Gauss-Newton.

#FORM='iEnS-Det-GN'           # Sqrt, determin, iter EnS

#FORM='EnRML-GN-obs'          # EnRML, Gauss-Newton
#FORM='EnRML-GN-state'        # Idem. Formulated in state space
#FORM='EnRML-GN-ens'          # Idem. Formulated in ensemble space

FORM='iEnS-GN'               # Sqrt, stoch, iter EnS. Equals EnRML-GN ? 


## Assimilation =============
nIter = 3
for k in range(nIter):
  A  = E - mean1(E)
  Eo = Obs(E)
  Z  = Eo - mean1(Eo)
  H  = Z@tinv(A)

  if FORM=='RML-GN':
    dLkl = zeros((M,N))
    dPri = zeros((M,N))
    for n in range(N): 
      Hn        = hp(tp(E[:,n]))
      Pn        = inv( inv(B0) + Hn.T@inv(R)@Hn )
      dLkl[:,n] = Pn@Hn.T@inv(R)@(y.ravel()-D[:,n]-Eo[:,n])
      dPri[:,n] = Pn@inv(B0)@(E0[:,n]-E[:,n])

  elif FORM=='MDA':
    D     = sqrtm(R)@randn((P,N))
    D    -= mean1(D)
    D    *= sqrt(N/N1)
    K     = A@Z.T@inv(Z@Z.T + nIter*N1*R)
    dLkl  = K@(y-sqrt(nIter)*D-Eo)
    dPri  = 0

  elif FORM=='iEnS-Det-GN':
    #Y     = H@A0
    Y     = Z @ inv(T)
    # Mean update
    Pw    = inv(Y.T@inv(R)@Y + N1*eye(N))
    dLkl  = Pw@Y.T@inv(R)@(y-mean1(Eo))
    dPri  = Pw@           (0-w)*N1
    w    += dLkl + dPri
    # CVar to ensemble space
    dLkl  = A0@dLkl
    dPri  = A0@dPri
    # Anomalies update (just add to dLkl)
    T     = funm_psd(Pw, sqrt)*sqrt(N1)
    dLkl  = dLkl + A0@T - A


  elif FORM=='EnRML-GN-obs':
    Y     = H@A0
    C     = Y@Y.T + N1*R
    K     = A0@Y.T@inv(C)
    dLkl  = K@(y-D-Eo)
    dPri  = (eye(M) - K@H)@(E0-E)
  elif FORM=='EnRML-GN-ens':
    Y     = H@A0
    Pw    = inv(Y.T@inv(R)@Y + N1*eye(N))
    K     = A0@Pw@Y.T@inv(R)
    dLkl  = K@(y-D-Eo)
    dPri  = (eye(M) - K@H)@(E0-E)
  elif FORM=='EnRML-GN-state':
    assert nla.matrix_rank(B0) == M # Not well-done by inv()
    P     = inv( inv(B0) + H.T@inv(R)@H )
    dLkl  = P@H.T@inv(R)@(y-D-Eo)
    dPri  = P@inv(B0)@(E0-E)

  elif FORM=='iEnS-GN':
    # Could re-construct W from E...
    # W     = tinv(A0) @ (E - x0) + Pi0C @ W 
    # ... but, of course, it's better to just use W to hold the state.

    # Regression_E: Y = Z tinv(A) A0:
    Y     = H @ A0                          # 
    # Y     = Z @ tinv( A ) @ A0              # 
    # Y     = Z @ tinv( A0 @ W @ AN ) @ A0    # 
    # Y     = Z @ tinv( Pi0 @ W @ AN )        # 

    # Regression_W: Y = Z tinv(T)   -- is it better?
    # Y     = Z  @ tinv( AN  @ W @ AN )       # 
    # Y     = Z  @ tinv( W @ AN )             # 
    # Y     = Eo @ tinv( W @ AN )             # 
    # Y     = Eo @ AN @ tinv( W ) @ AN        # 
    # Y     = Z  @ tinv( W ) @ AN             # 

    # Geir-Evensen forms -- equivalent to Regression_W forms
    # Y     = Z @ tinv( A ) @ A @ inv(Om)     # 
    # Y     = Z @ inv(Om)                     # 
    # Y     = Eo @ AN @ inv(Om)               # 

    # The rest of the algo
    Pw    = inv(Y.T@inv(R)@Y + N1*eye(N))
    dLkl  = Pw@Y.T@inv(R)@(y-D-Eo)
    dPri  = Pw@(eye(N)-W)*N1
    W    += dLkl + dPri
    We    = W - eye(N)
    Om    = eye(N) + We@AN
    # CVar to ensemble space
    dLkl  = A0@dLkl
    dPri  = A0@dPri


  E = E + dLkl + dPri

  # Animation
  if not (k+1)%nIter:
  #if not k%1:
    ax1.set_title(FORM+', k = '+str(k)+'. Plotting for state index [0]')
    if 'hist' in locals():
      for patch in hist:
        patch.remove()
    _,_,hist = ax1.hist(E[0,:],bins=linspace(*xx[[0,-1]],nbins),
        normed=True, label='$E_0$',color='r',alpha=0.6)
    plt.pause(0.5)


## Print result =================
# For (comparing methods) debugging
print("%15.15s, E_k state[0...5]:"%FORM, E[0,:5])


