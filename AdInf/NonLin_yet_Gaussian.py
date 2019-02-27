# Experiments with a non-linear transformation (NLT)
# as the dynamical model that has the property that
# Y = NLT(X) is Gaussian if X is Gaussian.
# 
# Shows that, if using a Sqrt EnKF:
# - Linearity abates sampling error:
#    - barB CV to B after initial transitory regime,
#    - and the ensemble itself also CVes.
#    - Both of these are analytical results of Bocquet/Carrassi/Grudzien!
# - Conversely, non-linearity creates sampling error,
#   as observed by the everpresent jitteriness in barB and barP,
#   even if Gaussianity is maintained.
# - Thus, sampling errors cause problems on their own,
#   and can be regarded as a separate issue from non-Gaussianity.
# - Sampling error causes jitteriness (in B itself) and P and xa.
#   Adapt.Infl can remedy jitteriness, since y also contains info on B.
# - Sampling error also fuels the negative bias (in P).
#   Adapt.Infl should remedy this too.
# - The "sampling error" due to the deterministic NLT behaves just like random noise:
#   for example, barP is negatively biased.
#   This may be ascertained by observing that NLT(sqrt(eN)*normalized(Ens)) has
#   expected sample var of 2 exactly (Ex(barB) = 2) and yet Ex(barP)<1,
#   however, the normalization is turned off for design simplicity.
# - More applications of model render the sample more noisy.
#   Thus, a longer DAW (e.g.) will cause a greater need for infl.
#
# Todo: M-dim case: NLT can be made for radius of multivariate Gaussian,
#       which is chi2(ddof=M).


###########################
# Preamble
###########################

from common import *
plt.style.use('AdInf/paper.mplstyle')

seed(88) # 4,59,88,89 are good ones


# The Non-Lin-Transf is a composite of familiar func's
NLT  = lambda x: sqrt(2)*ss.norm.ppf(ss.chi2.cdf(x**2,1))
Lin  = lambda x: sqrt(2)*x
aprx = lambda x: sqrt(2)*(-0.40 + 0.88*abs(zz) + 0.23*log(zz**2))
# aprx found by fitting a + b|x| + c log x²

def normalize(xx):
  N  = len(xx)
  mu = sum(xx)/N
  xx = xx - mu
  s2 = sum(xx**2)/(N-1)
  xx = xx/sqrt(s2)
  return xx

xlim  = 3.1
ylim  = xlim*sqrt(2)


###########################
# Histograms -- one axes
###########################
nBins = 30
nHist = int(1e5)

xx = randn(nHist)
yyNLT = NLT(xx)
yyLin = Lin(xx)

plt.figure(1,figsize=(5,sqrt(2)*5)).clear()
fig, ax0 = plt.subplots(num=1)
axX = ax0.twinx()
axY = ax0.twiny()
axX.set_yticks([])
axY.set_xticks([])
axY.grid(None)
axX.grid(None)

hhx, bbx, ppx    = axX.hist(xx   ,bins=nBins,range=(-xlim,xlim))
hhy, bby, ppy    = axY.hist(yyNLT,bins=nBins,range=(-ylim,ylim),orientation="horizontal")
_, _, ppyLin     = axY.hist(yyLin,bins=bby                     ,orientation="horizontal")

# Coloring
cm = plt.cm.get_cmap('jet')
customize = lambda x: x**0.8 # shifts cmaps
# Indices of inverse mapping
mids = 0.5 * (bbx[1:] + bbx[:-1])
inds = sorted(range(nBins), key=lambda k: NLT(mids[k]) + 0.01*np.sign(mids[k]))
# Bin coloring
for i,pr in enumerate(ppx   ): pr.set_color(cm( customize(i      /(nBins-1)) )); pr.set_ec('k'); pr.set_lw(0.6)
for i,pr in enumerate(ppy   ): pr.set_color(cm( customize(inds[i]/(nBins-1)) )); pr.set_ec('k'); pr.set_lw(0.6)
for i,pr in enumerate(ppyLin): pr.set_color(cm( customize(i      /(nBins-1)) )); pr.set_ec('k'); pr.set_lw(0.6)

axX.set_ylim(top  =5*axX.get_ylim()[1]) # downscale histogram
axY.set_xlim(right=5*axY.get_xlim()[1]) # downscale histogram
axY.invert_xaxis()

ax0.yaxis.tick_right()
ax0.set_xlim(-xlim,xlim)
ax0.set_ylim(-ylim,ylim)
zz  = linspace(-xlim,xlim,20001)
zz  = zz[zz!=0]
line_Lin, = ax0.plot(zz,Lin (zz), 'k--', lw=3, label=r'$\mathcal{M}_{\rm Lin}$')
line_NLT, = ax0.plot(zz,NLT (zz), 'k-' , lw=3, label=r'$\mathcal{M}_{\rm NonLin}$')
# line_aprx, = ax0.plot(zz,aprx(zz), 'k:' , lw=3, label=r'$\mathcal{M}_{\rm Approx}$')
ax0.legend(loc="upper center")
ax0.set_zorder(99)
ax0.patch.set_alpha(0)
# ax0.grid("on",color='k')
ax0.set_xticks(np.array([-1, 0, 1]).astype(int))
ax0.set_yticks(np.array([-sqrt(2), 0, sqrt(2)]))
ax0.set_yticklabels(["$-\sqrt{2}$","$0$","$+\sqrt{2}$"])
ax0.grid(None)

# Animate
if True:
  savefig_n.index = 1
  def save():
    savefig_n('data/AdInf/NonLin_yet_Gaussian/models_n', ext=".eps")

  def leg(*handles):
    ax0.legend(handles, [x.get_label() for x in handles], loc="upper center")

  # Hide all
  toggle_viz(*(ppx + ppyLin + ppy), line_Lin, line_NLT, pause=False)
  # Play/save
  toggle_viz(line_Lin);  leg(line_Lin)          ; save()
  toggle_viz(*ppx)    ;  pass                   ; save()
  toggle_viz(*ppyLin) ;  pass                   ; save()
  toggle_viz(line_NLT);  leg(line_Lin, line_NLT); save()
  toggle_viz(*ppyLin)
  toggle_viz(*ppy)    ;  pass                   ; save()


###########################
# EnKF experiment
###########################
R  = 2
P  = 1
y  = 0

N  = 40
eN = 1 + 1/N
E0 = NLT(randn(N))
nT = 100
#nT = 10000

# NLT
###########################
E  = E0.copy()
PP = zeros(nT)
BB = zeros(nT)
bb = zeros(nT)
for k in arange(nT):
  b = mean(E)
  A = E-b
  B = sum(A**2)/(N-1)
  K = B/(B+R)
  T = sqrt(1 - K)
  xa= b + K*(y - b)
  E = xa + T*A
  bb[k] = xa
  BB[k] = B
  PP[k] = T**2*B
  #E = sqrt(eN) * normalize(E) # Yields Ex(barB) = 2 (exactly)
  E = NLT(E)   # DEFAULT
  #E = sqrt(2)*randn(N) # To estimate bias with fully random ensemble
  
plt.figure(2)
plt.clf()

#matplotlib.colors.ColorConverter.colors['mc1'] = (0.976,0.333,0.518)
plt.plot(BB,color='mlo',lw=1.5)
plt.plot(PP,color='mlb',lw=1.5)
plt.plot(bb,color='k'  ,lw=1.5)

plt.ylim((-0.3,3))
plt.xlim(xmin=0,xmax=nT)

print('mean B:', mean(BB))
print('mean P:', mean(PP))

B_true          = 2 # known analytically
# I believe the following two would have been equal
# if not for the filtering removing some sampling error in B.
expected_err_B  = 2/N * 2**2 # Wick's theorem (Boc'2011)
empiricl_err_B  = np.var(BB-B_true,ddof=1)

# Lin
###########################
#plt.figure(3)
#plt.clf()
plt.figure(2)

E  = E0.copy()
PP = zeros(nT)
BB = zeros(nT)
bb = zeros(nT)
for k in arange(nT):
  b = mean(E)
  A = E-b
  B = sum(A**2)/(N-1)
  K = B/(B+R)
  T = sqrt(1 - K) # Yields the same as left-multiplying sym sqrt
  Aa= T*A
  xa= b + K*(y - b)
  E = xa + Aa
  bb[k] = xa
  BB[k] = B
  PP[k] = (1-K)*B
  E = Lin(E)

plt.plot(BB,color='mlo')
plt.plot(PP,color='mlb')
plt.plot(bb,color='k')


# Lin - PertObs
###########################
# E  = E0.copy()
# PP = zeros(nT)
# BB = zeros(nT)
# bb = zeros(nT)
# for k in arange(nT):
#   b = mean(E)
#   A = E-b
#   B = sum(A**2)/(N-1)
#   K = B/(B+R)
#   E = E + K*(y+sqrt(R)*randn(N)-E)
#   bb[k] = mean(E)
#   BB[k] = B
#   PP[k] = np.var(E,ddof=1)
#   E = Lin(E)
# 
# plt.plot(BB,color='k',alpha=0.3)
# plt.plot(PP,color='k',alpha=0.3)
# plt.plot(bb,color='k',alpha=0.3)



# Tune figure
###########################
plt.yticks(arange(0,4),('0.0','1.0','2.0','3.0'))
plt.xticks([0,nT])
plt.text(1.01*nT,0-0.09,'$\\bar{\mathbf{x}}^a$' )
plt.text(1.01*nT,1-0.09,'$\\bar{\mathrm{P}}^a$' )
plt.text(1.01*nT,2-0.09,'$\\bar{\mathrm{B}}$'   )
plt.text(0.5,-0.05,'DA cycle (i.e. time) index',
    fontsize=mpl.rcParams['axes.labelsize'],
    horizontalalignment='center',
    transform=plt.gca().transAxes
    )

ax = plt.gca()
ax.yaxis.tick_left()


##
#fname='NonLin_yet_Gaussian'
#plt.savefig('data/AdInf/figs/'+fname+'.eps')
#plt.savefig('data/AdInf/figs/'+fname+'.pdf')


###########################
# QUIT
###########################
sys.exit(0)


###########################
# Mini experiment 1
###########################
# Show that NLT(normalize(E))
# is slightly biased unless compensated for by sqrt(eN).

N  = 10
eN = 1 + 1/N
B  = 0
P  = 0
nT = 10000
for k in range(nT):
  E = randn(N)
  E = sqrt(eN) * normalize(E)
  b = sum(E)/N
  A = E-b
  B += sum(A**2)/(N-1) / nT

  E = NLT(E)
  b = sum(E)/N
  A = E-b
  P += sum(A**2)/(N-1) / nT
print()
print('mean var before NLT:', B)
print('mean var after  NLT:', P)


###########################
# Mini experiment 2
###########################
# How much "sampling error" does NLT cause?
# I.e. How inexact (as measured by its variance)
#      is NLT(E) if E is an exact sqrt.?

def NLT_recursive(x,level=1):
  if level == 0:
    return x
  x = NLT_recursive(x,level-1)
  x = ss.norm.ppf(ss.chi2.cdf(x**2,1)) # don't use sqrt(2)
  return x


PP = zeros(nT)
for k in range(nT):
  E = randn(N)
  E = sqrt(eN) * normalize(E)    # Remove "sampling error"
  E = NLT_recursive(E,2)
  b = sum(E)/N
  A = E-b
  PP[k] = sum(A**2)/(N-1)

print()
print('mean varance',mean(PP)) # Should be approx 1.
print('var of var  ',np.var(PP)) 
# This should increase with the num of NLT applied.
# Tops out at 0.22 (for recursion level=0 or >5). Min=0.15 (for level=1).



