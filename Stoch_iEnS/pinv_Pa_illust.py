# Illustrate the result of using the state-form KF update in the EnKF with pinv.
# Also see theoretical note deriving these results on paper.

from common import *

sd0 = seed_init()

GUI = False # else make animation

# N=2 gives the interesting tinv case.
# Use N>2 for validation purposes.
N  = 2 
N1 = N-1

xlims = array([-2,7])

y = array([[3],[4]])

# Setup figure
plt.figure(1,figsize=(6,5)).clear()
fig, ax = plt.subplots(num=1,subplot_kw={})

####################
# Slider controls
####################
if GUI:
  fig.subplots_adjust(bottom=0.2,top=1-0.05)
  from matplotlib.widgets import Slider
  x1 = Slider(fig.add_axes([1/8 , 1.0/9, 2/8 , 0.6/9]), '$E^f_1$ offset', -3, 7, y[0])
  x2 = Slider(fig.add_axes([1/8 , 0.2/9, 2/8 , 0.6/9]), '$E^f_2$ offset', -2, 2, 0)
  th = Slider(fig.add_axes([6/10, 0.2/9, 3/10, 0.6/9]), '$\\theta$', 0, 180, 0 )
else:
  x1 = Bunch(val=y[0])
  x2 = Bunch(val=0)
  th = Bunch(val=0)

hh = dict()
output = [None]

def foo(_):
  # ------ Obs ------
  t   = np.radians(th.val)
  U   = array([[ cos(t),  -sin(t)],
               [ sin(t),   cos(t)]])
  r1  = 6
  r2  = 1
  R   = U @ diag([r1,r2]) @ U.T
  iR  = inv(R)

  # ------ Prior ------
  # Prior ensemble and moments
  E = np.vstack([array([-1,1]) + x1.val,
                 array([0,0])  + x2.val])
  if N>2:
    # Total randomness -- no slider control!
    E = 1 + randn((2,N))     

  xb  = E.mean(axis=1,keepdims=True)
  X   = E - xb
  Bb  = X @ X.T / N1

  # ------ Posterior ------
  G   = inv(N1*eye(N) + X.T@iR@X)
  Xa  = X @ sqrtm(G) * sqrt(N1)
  Pb  = Xa @ Xa.T / N1
  KG  = X @ G @ X.T@iR
  xab = xb + KG@(y - xb)
  Ea  = xab + Xa

  # ------ Pinv (tinv) version ------
  Pt    = inv( tinv(Bb) + iR )
  xt    = Pt @ (tinv(Bb)@xb + iR@y) # i.e.   xt = xb + KG@(y - xb),   with KG = Pt @ iR
  # Trigonometric version
  Be    = Bb[0,0] # Marg
  R11   = R[0,0]  # Marg
  K11   = Be / (Be + R11) # Marginal KG. Not quite the same as KG[0,0].
  xt_1  = xb[0] + K11 * (y[0] - xb[0]) # Compare to xab[0].
  xt_2  = y[1] + (y[0] - xb[0]) * (r2-r1) / (Be + R11) * sin(2*t)/2
  xt    = np.array([xt_1, xt_2])
  # A choice must be made for the ensemble (coz it does not have enough DoF to "span" Pt)...
  Et    = xt + CovMat(Pt).sym_sqrt @ CovMat(Bb).sym_sqrt_inv @ X # ... One possible choice.
  

  theta2slope = lambda h: tan(pi * h.angle / 180)
  slope_pinv  = lambda t: (r1-r2)/R11*sin(2*t)/2

  def contour1(mu, C,col,**kwargs):
    return cov_ellipse(ax, mu, C, edgecolor=col,lw=4, facecolor='none',**kwargs)

  if len(hh):
    hh['E_y'].set_offsets(y.ravel())
    hh['E_f'].set_offsets(E.T)
    hh['E_a'].set_offsets(Ea.T)
    hh['E_t'].set_offsets(Et.T)


    hh['c_y'].remove(); hh['c_y'] = contour1(y  , R ,'y')
    hh['c_f'].remove(); hh['c_f'] = contour1(xb , Bb,'b')
    hh['c_a'].remove(); hh['c_a'] = contour1(xab, Pb,'r',zorder=3)
    hh['c_t'].remove(); hh['c_t'] = contour1(xt , Pt,'m',zorder=4)

    hh['l_y'].set_ydata( y[1,0] + theta2slope(hh['c_y']) *(xlims- y[0,0]))
    hh['l_t'].set_ydata(xt[1,0] + theta2slope(hh['c_t']) *(xlims-xt[0,0]))
    # hh['l_k'].set_ydata(xt[1,0] + slope_pinv(t)          *(xlims-xt[0,0]))

  else:
    symbolic = True
    hh['E_y'] = ax.scatter(y[0], y[1] ,15**2,'y','*', label='$\mathcal{N}(y|x)$' if symbolic else 'Likelihood')
    hh['E_f'] = ax.scatter(E[0] ,E[1] ,15**2,'b',     label='$\\bar{P}^f$'       if symbolic else 'Prior (ens and contour)')
    hh['E_a'] = ax.scatter(Ea[0],Ea[1],15**2,'r',     label='$\\bar{P}^a$'       if symbolic else 'Updated using K-gain.')
    hh['E_t'] = ax.scatter(Et[0],Et[1],15**2,'m',     label='$\\bar{P}^a_+$'     if symbolic else 'Updated using p-inv.')

    hh['c_y'] = contour1(y  , R , 'y'         )
    hh['c_f'] = contour1(xb , Bb, 'b'         )
    hh['c_a'] = contour1(xab, Pb, 'r',zorder=3)
    hh['c_t'] = contour1(xt , Pt, 'm',zorder=4)

    hh['l_y'], = ax.plot(xlims,  y[1,0] + theta2slope(hh['c_y']) *(xlims- y[0,0]), 'y',alpha=0.7)
    hh['l_t'], = ax.plot(xlims, xt[1,0] + theta2slope(hh['c_t']) *(xlims-xt[0,0]), 'm',alpha=0.7)
    # hh['l_k'], = ax.plot(xlims, xt[1,0] + slope_pinv(t)          *(xlims-xt[0,0]), 'k',lw=0.7)

  output[0] = t, U, r1, r2, R, iR, E, xb, X, Bb, G, Xa, Pb, KG, xab, Ea, Pt, xt, Et, xt_1, xt_2

foo(None)

# Adjust
ax.set_aspect('equal')
ax.set_xlim(xlims)
ax.set_ylim(xlims+.5)
ax.set_xlabel('$x_1$')
ax.set_ylabel('$x_2$')
ax.set_xticklabels([])
ax.set_yticklabels([])
ax.set_title('EnKF analysis update illustration')
ax.legend(fontsize=12,loc='upper left')

if GUI:
  # Activate listeners
  x1.on_changed(foo)
  x2.on_changed(foo)
  th.on_changed(foo)
else:
  savefig_n.index = 1
  def save(param,val):
    param.val = val
    foo(None)
    savefig_n('data/Stoch_iEnS/figs/prez/pinv_Pa_illust_n')
  ax.set_title('')

  mpl.rcParams['savefig.bbox'] = 'tight'
  for val in ccat(arange(3,-1,-.5),arange(-1,6,1)):                save(x1,val)
  for val in ccat(arange(0,-2,-.5),arange(-2,2,1),arange(2,0,-1)): save(x2,val)
  for val in ccat(arange(6,3,-1)):                                 save(x1,val)
  th.val = 30
  for val in ccat(arange(3,-1,-.5),arange(-1,6,1)):                save(x1,val)



