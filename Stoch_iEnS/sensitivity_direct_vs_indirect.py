# Try to explain which is better of:
#   - explicit regression (to obtain H from state to obs space)
#   - ensemble-estimation of HBH and BH, which, via Woodbury, can be shown to be
#     direct regression to the ensemble.
# Recall that these only differ in the case of nonlinear h() and N-1 > M.
# Context: (single-cycle update of non-iterative) EnKF.

from common import *
from scipy.stats import norm
from scipy.optimize import fsolve

sd0 = seed_init(4)

N     = 5              # Ens size
N1    = N-1               #
bw    = N**-0.2
nbins = int(12/bw) # Histogram bins
Pi1   = ones((N,N))/N
PiAN  = eye(N) - Pi1



## Prior ===================
b  = 1.5
B  = 3
E  = b + sqrt(B)*randn(N)


# NB NB NB
E.sort()


hcx = 1/sqrt(2*pi*B)


## Obs  ===================
y = 40
R = 16
hcy = 1/sqrt(2*pi*R)
# == Non-Lin H ==:
def  h(x): return x**2
def hp(x): return 2*x


#################################
# Plot setup, scatter plot
#################################
plt.figure(1).clear()

ax_s, ax_x, ax_y = axes_with_marginals(4, 1)
ax_x.set_yticklabels([])
ax_y.set_xticklabels([])
ax_s.scatter(E, h(E), 14**2, color='b', label='E',zorder=4)

xX = [-5, 5] # ax_s.get_xlim()
yY = [-5, 10] # ax_s.get_ylim()
#xX = [-4, 10] # ax_s.get_xlim()
#yY = [-5, 70] # ax_s.get_ylim()
xx = linspace(*xX,2001)
yy = linspace(*yY,2001)

ax_s.plot(   xx, h(xx),       'k', label='h')

ax_s.set_xlim(xX)
ax_s.set_ylim(yY)

## PDFs ===================
dx = xx[1]-xx[0]
def normlz(pp):
  return pp / sum(pp) / dx

prior_xx = norm.pdf(xx,b,sqrt(B))
ax_x.plot(xx,prior_xx, 'b' )

## Histograms ===================

#ax_x.hist(  E ,bins=nbins, normed=True)
#ax_y.hist(h(E),bins=nbins, normed=True, orientation="horizontal")


## Stem plotting ===================
def stem_pieces(xx):
  N = len(xx)
  N1 = N-1

  # Create piece-wise line
  pw_x = array([ [xx[k], xx[k], xx[k+1] ] for k in range(N1)]).ravel()
  pw_y = array([ [0    , 1    , nan     ] for k in range(N1)]).ravel()
  # Add the last stem
  pw_x = ccat( pw_x, xx[N1], xx[N1] )
  pw_y = ccat( pw_y, 0     , 1      )

  return pw_x, pw_y

pw_x, pw_y = stem_pieces(E)
pw_y = 0.7 * hcx * pw_y
ax_x.plot(pw_x, pw_y ,alpha=bw, lw=1)

pw_x, pw_y = stem_pieces(h(E))
pw_y = 0.7 * hcy * pw_y
ax_y.plot(pw_y, pw_x ,alpha=bw, lw=1)


## X ===================

def arrow(x0,y0,x1,y1,**arrowprops):
  ax_s.annotate("",arrowprops=dict(arrowstyle="->",**arrowprops),
      xy=(x1, y1), xytext=(x0, y0),
      annotation_clip=False)

def norm_with_AR(xy):
  "Norm of vector scaled by plot aspect_ratio"
  return sqrt( xy[0]**2 + (xy[1]*AR)**2 )

# Plot ensemble mean
#ax_s.plot(E.mean(), h(E).mean(), 'r*',ms=12)

# Graphics length scales and aspect ratio
wx = xX[1]-xX[0]
wy = yY[1]-yY[0]
AR = wx/wy

# Ensemble (average) gradient
X   = E    @ PiAN
b_e = E.mean()
Y   = h(E) @ PiAN
y_e = h(E).mean()
B_e = X@X / N1
PiX = tinv(X[None,:]) @ X[None,:]
H   = Y@X / (X@X) # Calc ensemble gradient
x0  = fsolve(lambda x: hp(x) - H, E.mean())  # Find x0 such that hp(x0) == H
arrow(x0 - 0.3*wx, h(x0) - 0.3*wx*H, x0 + 0.3*wx,h(x0) + 0.3*wx*H, lw=2) # Draw

for n in range(N):
  xy0 = array([E[n] , h(E[n])])
  dxy = array([1    , hp(E[n])])
  dxy = dxy * 0.2*wx / norm_with_AR(dxy)
  xy1 = xy0 + dxy
  arrow(*xy0, *xy1, color='g',lw=2)



## X ===================

LinXY = lambda x: y_e + H*x
LinWY = lambda w: y_e + Y@w
LinWX = lambda w: b_e + X@w


N2 = 20000
W2 = randn((N,N2))/sqrt(N1)
ax_s.scatter(LinWX(W2), LinXY(LinWX(W2)), color='y', alpha=N2**-0.2, label='E2 - LinX', zorder=2)
ax_s.scatter(LinWX(W2), LinWY(W2)       , color='r', alpha=N2**-0.2, label='E2 - LinW', zorder=1)

ax_x.hist(LinWX(W2)        ,bins=nbins, normed=True, color='c')
ax_y.hist(LinXY(LinWX(W2)) ,bins=nbins, normed=True, color='y', alpha=0.7, orientation="horizontal")
ax_y.hist(LinWY(W2)        ,bins=nbins, normed=True, color='r', alpha=0.4, orientation="horizontal")


bb = ax_y.boxplot([LinXY(LinWX(W2)), LinWY(W2)],
    notch=0,sym='',vert=True,
    whis=[5, 95],medianprops={'linewidth':0},
    showmeans=True,meanline=True,meanprops={'color':'k','linestyle':'-'},
    patch_artist=True,widths=0.7)



for (c, patch) in zip(['r','y'], bb['boxes'][::-1]):
  patch.set(facecolor=c,alpha=0.8)







##
