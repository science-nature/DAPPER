# Test if chain rule holds for LAD: least absolute deviation.
# Other ideas:
#  - min-max: should hold for monotonic functions.
#  - least-cubes, etc...

from common import *
from scipy.stats import norm

sd0 = seed_init(5)

def mean1(E): return mean(E,axis=1,keepdims=True)
def     a(E): return E - mean1(E)

M  = 3   # State length
N  = 10   # Ens size
N1 = N-1 #

E0 = 1 + randn((M,N))
xX = array([-4,4])

## Funcs ===================
def f1(x): return 0.4*x**2
def f2(x): return sin(x)
def  f(x): return f2(f1(x))

E1 = f1(E0)
E2 = f2(E1)

## Least squares ===================
LS_1 = a(E1) @ tinv(a(E0))
LS_2 = a(E2) @ tinv(a(E1))
LS_c = LS_2 @ LS_1

LS_d = a(E2) @ tinv(a(E0))

spell_out(LS_d)
spell_out(LS_c)

## Least squares with regularization ===================

eps = 0.1
b   = zeros(M)

def LS_reg(Y,X):
  X = a(X)
  Y = a(Y)
  return inv( X@X.T + eps ) @ (X@Y.T + eps*b )

LS_reg_d = LS_reg(E2, E0)
LS_reg_c = LS_reg(E2, E1) @ LS_reg(E1, E0)

spell_out(LS_reg_d)
spell_out(LS_reg_c)

## Least absolute deviations ===================
# Trying to follow statsmodels.org/dev/examples/notebooks/generated/quantile_regression.html
# but what I got so far does not seem promising...
# import statsmodels.formula.api as smf
# import pandas
# EE_marg = np.vstack( [E0[0,:], E2[0,:]] ).T
# df = pd.DataFrame(data=EE_marg, columns=['x','y'])
# mod = smf.quantreg('y ~ x', df)
# res = mod.fit(q=0.5)
# res.summary()

## Plotting ===================
fig = plt.figure(1)
fig.clear()
fig, (ax1, ax2) = plt.subplots(2,1,sharex=True,num=1, gridspec_kw = {'height_ratios':[5, 1]})

xx = linspace(*xX,1001)
ax1.plot(xx,f(xx), label=f.__name__)
ax1.plot(E0[0,:],f(E0)[0,:], 'o')
ax1.legend()
ax1.set_title('Dim 0')
ax2.hist(E0[0,:],bins=linspace(*xX,40),normed=True)
