# Test if the "sensitivity" (based on exact covs) equals the mean gradient
# - Gaussian distributions;   answer: Yes [by Stein's identity]
# - Non-Gaussian, with the nonlinear Obs.mod being a polynomial of degree:
#   - 2: answer: Only if skew=0. But the sensitivity *at* the mean is the same anyways.
#   - 3: answer: No -- even with skew = 0

##

from common import *
import seaborn as sns

#sd0 = seed_init(2)

## Prior
b = 0
B = 2

# For this study, N should be "practically infinite"
N  = 10**7
N1 = N-1


## h
R  = 0
# def h0(x): return -0.1*x + x**2
# def h1(x): return -0.1 + 2*x
def h0(x): return -0.1*x + x**2 - 0.2 *x**3
def h1(x): return -0.1 + 2*x    - 0.6 *x**2
# def h0(x): return 4*x
# def h1(x): return 4


## Ensemble
TriMod = False                
if TriMod: # Non-Gaussian
  wN = [int(w*N) for w in [0.5, 0.1876]]; wN.append(N-sum(w))
  xx = ccat(    xx[:,:wN[0]],
     3 + 0.5*randn((1,wN[1])),
    -4 + 0.5*randn((1,wN[2])),
    axis=1)
else: # Gaussian
  xx = b + sqrt(B)*randn((1,N)) 

prnt = lambda s,v: print('%13.13s: %.5f'%(s,v))
prnt("mean(xx)", mean(xx) )
prnt(" var(xx)", np.var(xx,ddof=1) )
prnt("skew(xx)", ss.skew(xx.T,bias=False) )


## Stats
yy = h0(xx) + sqrt(R)*randn((1,N))

Cyx = np.cov(yy,xx,ddof=1)[0,1]
Cx  = np.var(xx,ddof=1)
H   = Cyx/Cx

prnt("\nH",H)
prnt("mean(h1(xx))", mean(h1(xx)) )
prnt("h1(mean(xx))", h1(mean(xx)) )


## Plotting
df = pd.DataFrame(ccat(xx,yy).T[np.random.choice(N,10**4)], columns=['x','y'])
sns.jointplot(x="x", y="y", data=df)

