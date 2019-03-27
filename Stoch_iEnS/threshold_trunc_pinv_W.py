# Investigate what goes wrong in iEnS when using
# Y = ME @ tinv(Wc, threshold=1-eps), and eps isn't minuscule.

# Conclusions (regarding truncation level and statistical errors):
# - For trunc<min(N1,M): Both 'Old' and 'New' versions yield non-convergent posteriors.
# - For M<trunc<N1: 'New' yields bad posteriors for iter>1. For 'Old' the question does not make sense.
# Note: Using actual threshold (between 0 and 1) confuses the question.
#       It is easier to study using a fixed (wrt iterations) trunc (number of s. vals).
# Curiosity: with M=2, N>3, and an unreasonably small threshold (e.g. 0.3),
#   the posteriors of even-numbered iterations equal the prior!
#   Explanation: It yields Y=0 ==> grad_y=0 and grad_b that exactly undoes W.
#   This also makes intuitive sense: zero relationship ==> zero update.

# Conclusions (regarding the iterative growth of numerical error):
# For numerical stability, need to center Y following M(E)@inv(T).
# For numerical+statistical stability (e.g. in case W isn't invertible,
# or you've truncated it), need to pre-multiply prior increment with
# by the projection undergone by the lklhd increment, U@U.T or Wc@tinv(Wc).


##
from common import *
sd0 = seed_init(3)

def mean(x):
  return np.mean(x,axis=1,keepdims=True)

##
nIter = 10   # Num of iterations
M     = 10   # Param size
N     = 30   # Ens size
N1    = N-1

# EnRML version
version = 'New'

# Truncation control
trunc = N-1 # by rank -- Try N-1, N-2, 1, etc for some strangeness.
# trunc = 0.5 # by threshold -- try 1-eps for eps small or large for some strangeness

# These are left so as to recognize the formulae, but are not of interest here.
Mod = lambda x: 10*x
y   = zeros((M,1))
D   = zeros((M,N))

##
# Init
E  = randn((M,N))
E0 = E.copy()
xb = mean(E)
X  = E-xb
W  = eye(N)
AN = eye(N) - 1/N

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
  ax.text(0.01, 0.7,i,transform=ax.transAxes,fontsize=16,color='r')
  plt.pause(0.01)
splot(0,E)


# Algorithm
for i in 1+arange(nIter):

  ME = Mod(E)

  if version=='New':
    # Y       = ME @ tinv(W - mean(W), trunc) # v1
    Y       = nla.solve(W.T, ME.T).T        # v2
    Y       = Y - mean(Y) # Required for v2. Also helps v1 (stability).
    grad_y  = Y.T @ (y + D - ME)
    grad_b  = (eye(N) - W)*N1
    iCw     = Y.T@Y + N1*eye(N)
    W       = W + nla.solve(iCw, grad_y + grad_b)
    E       = xb + X@W

  if version=='GE':
    const   = 1 # impacts numerical errors => rate of convergence
    Wc      = W - mean(W) + const # coz Omeg = I + (W-I)@AN = W@AN + 1
    Y       = nla.solve(Wc.T, ME.T).T
    grad_y  = Y.T @ (y + D - ME)
    grad_b  = (eye(N) - W)*N1
    iCw     = Y.T@Y + N1*eye(N)
    W       = W + nla.solve(iCw, grad_y + grad_b)
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

  # with printoptions(precision=10):
    # print("\n",i)
    # am = lambda x: abs(x).max()
    # print(am(grad_y + grad_b))
    # # print(am(nla.solve(iCw, grad_y + grad_b)))

  # print("Iter", i, "Mean", mean(E), "Cov:", np.cov(E.T,ddof=1), sep='\n')
  splot(i,E)


##
##


