# Investigate what goes wrong in iEnS when using
# Y = ME @ tinv(Wc, threshold=1-eps)
# and the threshold isn't minuscule (but still reasonably small).

# Curiosity: with M=2, N>3, and an unreasonably small threshold (e.g. 0.3),
# the posteriors of even-numbered iterations equal the prior!
# Explanation: It yields Y=0 ==> grad_y=0 and grad_b that exactly undoes W.
# This also makes intuitive sense: zero relationship ==> zero update.

# Note: This is not an asymptotic problem. It does not get worse as nIter-->infty.

# Note: Using actual threshold (between 0 and 1) confuses the question.
#       It is easier to study using a fixed (wrt iterations) trunc (number of s. vals).

# Conclusions:
# - For M==2<N1:
#   - If trunc=N1: all good, even for large N (slow ~1min to compute for N>1e3)
#   - If trunc<N1: 
#      - 'New' version yields bad posteriors for some iter>1,
#        sometimes in an alternating (between iterations) fashion.
#        It generally sucks for any trunc<N1.
#      - 'Old' version: question does not make sense,
#        since any trunc>M still yields truncation at M.
# - For 2<M<trunc<N1: the alternations are more spaced out, but still happen (unless trunc=N1).
# - For trunc<N1<M: Both 'Old' and 'New' versions yield non-convergent posteriors.

##
from common import *
sd0 = seed_init(3)

def mean(x):
  return np.mean(x,axis=1,keepdims=True)

##
nIter = 10    # Num of iterations
M     = 2    # Param size
N     = 3   # Ens size
N1    = N-1

# EnRML version: 'New' or 'Old'
version = 'New'

# Truncation control
trunc = N-2 # by rank -- try N-1, N-2, 1, etc
# trunc = 0.5 # by threshold -- try 1-eps for eps small or large

# These are left so as to recognize the formulae, but are not of interest here.
Mod = lambda x: x
y   = zeros((M,1))
D   = zeros((M,N))

##
# Init
E  = randn((M,N))
E0 = E.copy()
xb = mean(E)
X  = E-xb
W  = eye(N)

# Plotting
fig, axs = freshfig(1,None,*nrowcol(nIter+1),sharex=True,sharey=True)
def splot(i,E):
  nP = min(N,100)
  ax = axs.flatten()[i]
  if M>1:
    ax.plot(E[0,:nP],E[1,:nP],'*')
  elif M==1:
    ax.plot(E[0,:nP],zeros(nP),'*')
    ax.set_ylim(-1,1)
  ax.set_title(i)
  plt.pause(0.01)
splot(0,E)


# Algorithm
for i in 1+arange(nIter):

  ME = Mod(E)

  if version=='New':
    Wc      = W - mean(W)
    # iWc   = tinv(Wc, trunc)
    U,ss,VT = tsvd(Wc, trunc)
    # ProjW   = U@U.T
    iWc     = (VT.T / ss) @ U.T
    Y       = ME @ iWc
    grad_y  = Y.T @ (y + D - ME)
    grad_b  = (eye(N) - W)*N1
    # grad_b  = ProjW @ (eye(N) - W)*N1 # safe
    iCw     = Y.T@Y + N1*eye(N)
    W       = W + nla.solve(iCw, grad_y + grad_b)
    Wc      = W - mean(W) # just for debugging
    E       = xb + X@W

  elif version=='Old':
    Xk     = E - mean(E)
    H      = ME @ tinv(Xk, trunc)
    Y      = H @ X
    C      = Y@Y.T + N1*eye(M)
    K      = X@Y.T@inv(C)
    dLkl   = K@(y-D-ME)
    dPri   = (eye(M) - K@H)@(E0-E)
    E      = E + dLkl + dPri

  # print("Iter", i, "Mean", mean(E), "Cov:", np.cov(E.T,ddof=1), sep='\n')
  splot(i,E)


##
##


