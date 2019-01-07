from common import *

#from mpl_toolkits.mplot3d import Axes3D
from mpl_toolkits.mplot3d.art3d import juggle_axes

from mpl_toolkits.axes_grid1 import make_axes_locatable
from mpl_toolkits.axes_grid1.inset_locator import inset_axes

from matplotlib import colors
from matplotlib.ticker import MaxNLocator


# TODO:
# - figure number management
# - for loops
class LivePlot:
  """
  Live plotting functionality.
  """
  def __init__(self,stats,key,E=None,P=None,only=False,**kwargs):
    """
    Initialize plots.
     - only: (possible) fignums to plot.
    """

    if isinstance(only,bool):
      only=range(99)
    #else: assume only is a list fo fignums

    # config = stats.config
    HMM = stats.HMM
    Nx  = HMM.Nx
    N   = len(E) if E is not None else None

    # Store
    self.HMM   = HMM
    self.stats = stats

    # Set up prompts
    self.is_on     = True
    self.is_paused = False
    print('Initializing liveplotting...')
    print('Press <Enter> to toggle live plot OFF/ON.')
    print('Press <Space> and then <Enter> to pause.')

    #ens_props = {} yields rainbow
    ens_props = {'color': 0.7*RGBs['w'],'alpha':0.3}

    # Diagnostics
    if 1 in only:
      self.sliding_diagnostics = sliding_diagnostics_LP(1,HMM.t,N,stats)
      set_figpos('2312')

    # Correlation plot
    if 2 in only and Nx<1001:
      self.correlations = correlations_LP(2,Nx,E,P)
      set_figpos('2311')

    # Spectral error plot
    if 4 in only:
      self.spectral_errors = spectral_errors_LP(4,stats)
      set_figpos('2311')

    # Weighted histogram
    if 4 in only and E is not None and stats._has_w:
      self.weight_histogram = weight_histogram_LP(4,stats.w[0])
      set_figpos('2321')

    # User-defined state
    if 9 in only and hasattr(HMM,'liveplotting'):
      LPs = HMM.liveplotting
      LPs = LPs if hasattr(LPs,'__len__') else [LPs]
      self.custom = [LP(stats,key,E,P) for LP in LPs]
      for i, (fig, _) in enumerate(self.custom):
        win_title(fig,"Custom plot %d"%i)
      set_figpos('2322')
      plot_pause(0.01)

    plot_pause(0.01)


  def update(self,key,E=None,P=None,**kwargs):
    """Update liveplots"""
    if self.skip_plotting(): return

    k,kObs,f_a_u = key

    # Diagnostics
    if hasattr(self, 'sliding_diagnostics'):
      fignum = self.sliding_diagnostics.fig.number
      if plt.fignum_exists(fignum):
        plt.figure(fignum)
        self.sliding_diagnostics.update(key,self.stats)
        plot_pause(0.01)

    # Correlation plot
    if hasattr(self, 'correlations'):
      fignum = self.correlations.fig.number
      if plt.fignum_exists(fignum):
        plt.figure(fignum)
        self.correlations.update(E,P)
        plot_pause(0.01)

    # Spectral error plot
    if hasattr(self, 'spectral_errors'):
      fignum = self.spectral_errors.fig.number
      if plt.fignum_exists(fignum) and self.spectral_errors.do_spectral_error:
        plt.figure(fignum)
        self.spectral_errors.update(k)
        plot_pause(0.01)

    # Weight histogram
    if kObs and hasattr(self, 'weight_histogram'):
      fignum = self.weight_histogram.fig.number
      if plt.fignum_exists(fignum):
        plt.figure(fignum)
        self.weight_histogram.update(self.stats.w[k])
        plot_pause(0.01)

    # User-defined state
    if hasattr(self,'custom'):
      for fig, updator in self.custom:
        if plt.fignum_exists(fig.number):
          plt.figure(fig.number)
          updator(key,E,P)
          plot_pause(0.01)


  def skip_plotting(self):
    """
    Poll user for keypresses.
    Decide on toggling pause/step/plot:
    """
    open_figns = plt.get_fignums()
    if open_figns == []:
      return True

    if self.is_paused:
      # If paused
      ch = getch()
      # Wait for <space> or <enter>
      while ch not in [' ','\r']:
        ch = getch()
      # If <enter>, turn off pause
      if '\r' in ch:
        self.is_paused = False

    key = poll_input() # =None if <space> was pressed above
    if key is not None:
      if key == '\n':
        # If <enter> 
        self.is_on = not self.is_on # toggle plotting on/off
      elif key == ' \n':
        # If <space>+<enter> 
        self.is_on = True # turn on plotting
        self.is_paused = not self.is_paused # toggle pause
        print("Press <Space> to step. Press <Enter> to resume.")
    return not self.is_on


# TODO:
# - rm estimate_good_plot_length in favour of Tplot
# - clean up HMM, stats, passing-around
# - iEnKS diagnostics don't work at all when store_u=False
# - mv label and shape to stats.py
# - re-use settings with plot_time_series
# - make set_limits for t-axis the same here and in LV
star = "${}^*$"
class sliding_diagnostics_LP:

  def __init__(self,fignum,chrono,N,stats):
      GS = {'left':0.125,'right':0.76}
      self.fig, (self.ax1, self.ax2) = \
          freshfig(fignum, (5,3.5), nrows=2, sharex=True, gridspec_kw=GS)

      self.legend_not_yet_created = True
      win_title(self.fig,"Scalar diagnostics")
      self.ax1.set_ylabel('RMS')
      self.ax2.set_ylabel('Values') 
      self.ax2.set_xlabel('Time (t)')

      def lin(a,b): return lambda x: a + b*x
      def divN()  : return lambda x: x/N
      def Id(x)   : return x

      # RMS
      d1 = {
          'rmse'    : [Id          , None   , dict(c='k'      , label='Error'            )],
          'rmv'     : [Id          , None   , dict(c='b'      , label='Spread', alpha=0.6)],
        }

      # OTHER         transf       , shape  , plt kwargs
      d2 = OrderedDict([
          ('skew'   , [Id          , None   , dict(c=     'g' , label=star+'Skew/$\sigma^3$'        )]),
          ('kurt'   , [Id          , None   , dict(c=     'r' , label=star+'Kurt$/\sigma^4$'        )]),
          ('trHK'   , [Id          , None   , dict(c=     'k' , label=star+'HK'                     )]),
          ('infl'   , [lin(-10,10) , 'step' , dict(c=     'c' , label='10(infl-1)'                  )]),
          ('N_eff'  , [divN()      , 'dirac', dict(c=RGBs['y'], label='N_eff/N'             ,lw=3   )]),
          ('iters'  , [lin(0,.1)   , 'dirac', dict(c=     'm' , label='iters/10'                    )]),
          ('resmpl' , [Id          , 'dirac', dict(c=     'k' , label='resampled?'                  )]),
        ])

      def init_ax(ax,style_table):
        plotted_lines = OrderedDict()
        for name in style_table:

            # SKIP -- if stats[name] is not (1) in existence, or (2) active.
            try: stat = getattr(stats,name) # (1)
            except AttributeError: continue
            try: val0 = stat[0] # (2)
            except KeyError: continue
            # PS: recall (from series.py) that even if store_u is false, stat[k] is
            # still present if liveplotting=True via the k_tmp functionality.

            ln = {}
            ln['transf'] = style_table[name][0]
            ln['shape']  = style_table[name][1]
            ln['plt']    = style_table[name][2]

            # Create series
            u_lag = K_lag if isinstance(stat,FAU_series) else 0
            ln['data'] = RollingArray(u_lag + a_lag)
            ln['tt']   = RollingArray(u_lag + a_lag)

            # Plot (init)
            ln['handle'], = ax.plot(ln['tt'],ln['data'],**ln['plt'])

            plotted_lines[name] = ln
        return plotted_lines


      t = chrono
      K_lag = estimate_good_plot_length(stats.xx,t,mult = 80)
      a_lag = int(K_lag/ t.dkObs)

      self.dt_margin = (t.tt[-1] - t.tt[-min(K_lag,len(t.tt))]) / 20
      self.chrono    = t

      self.d1        = init_ax(self.ax1, d1);
      self.d2        = init_ax(self.ax2, d2);
      self.update((0,None,'u'), stats)


  def update(self,key,stats):
      k, kObs, f_a_u = key

      def update_arrays(plotted_lines):
        for name, ln in plotted_lines.items():
          stat = getattr(stats,name)
          if isinstance(stat,FAU_series):
            # The stat series is an FAU_series, defined for each t in tt.
            # ln['data'] (a RollingArray) *will* include duplicate instances for f/a times.
            ln['tt']  .update(k   , t.tt[k])
            ln['data'].update(k   , ln['transf'](stat[k]))
          elif 'a' in f_a_u:
            # The stat series is an ndarray, defined for each t in ttObs.
            # ln['data'] (a RollingArray) *will not* have duplicates, coz only 'a' is input.
            ln['tt']  .update(kObs, t.ttObs[kObs])
            ln['data'].update(kObs, ln['transf'](stat[kObs]))


      def update_plots(ax,plotted_lines):

        def bend_into(shape, xx, yy):
          # Get arrays. Repeat (to use for intermediate nodes). 
          yy = yy.array.repeat(3)
          xx = xx.array.repeat(3)
          if shape == 'step':
            yy = np.hstack([yy[1:], nan]) # roll leftward
          elif shape == 'dirac':
            nonlocal nDirac
            axW      = np.diff(ax.get_xlim())
            yy[0::3] = False          # set "curve" to 0
            xx[2::3] = nan            # make "curve" disappear
            xx      += nDirac*axW/100 # offset curve horizontally
            nDirac  +=1
          return xx, yy

        nDirac = 0
        for name, ln in plotted_lines.items():
          ln['handle'].set_data(*bend_into(ln['shape'], ln['tt'], ln['data']))


      def adjust_ax(ax,reference_line):
          t1 = reference_line['tt'].leftmost()
          t2 = t.tt[k] + self.dt_margin
          ax.set_xlim(t1, t2)


      def do_once(ax,plotted_lines):
          # Rm lines that only contain NaNs
          for name in list(plotted_lines):
            ln = plotted_lines[name]
            if not np.any(np.isfinite(ln['data'])):
              ln['handle'].remove()
              del plotted_lines[name]
          # Add legends
          if plotted_lines:
            ax.legend(loc='upper left', bbox_to_anchor=(1.01, 1),borderaxespad=0)


      t   = self.chrono
      ax1 = self.ax1
      ax2 = self.ax2

      update_arrays(    self.d1)
      update_arrays(    self.d2)
      update_plots(ax1, self.d1)
      update_plots(ax2, self.d2)

      # Set x-limits (time) 
      adjust_ax(ax1, self.d1['rmse'])
      # Set y-limits
      data1 = [ln['data'].array for ln in self.d1.values()]
      data2 = [ln['data'].array for ln in self.d2.values()]
      ax1.set_ylim(0, d_ylim(data1, ax1,                cC=0.2,cE=0.9)[1])
      ax2.set_ylim(  *d_ylim(data2, ax2, Max=4, Min=-4, cC=0.3,cE=0.9))

      # Init legend. Rm nan lines. 
      if  self.legend_not_yet_created and k>t.kkObs[0]: # Don't use == (fails if user skipped)
          self.legend_not_yet_created = False
          do_once(ax1,self.d1)
          do_once(ax2,self.d2)

          ax2.annotate(star+": mean of\nmarginals", xy=(0,-1.5/len(self.d2)),
              xycoords=ax2.get_legend().get_frame(),
              bbox=dict(alpha=0.0), fontsize='small')



class weight_histogram_LP:

  def __init__(self,fignum,w0):

    fig, ax = freshfig(fignum, (6,3), gridspec_kw={'bottom':.15})
    win_title(fig,"Weight histogram")
    ax.set_xscale('log')
    ax.set_xlabel('weigth [× N]')
    ax.set_ylabel('count')
    if len(w0)<10001:
      hist   = ax.hist(w0)[2]
      N      = len(w0)
      xticks = 1/N * 10**arange(-4,log10(N)+1)
      xtlbls = array(['$10^{'+ str(int(log10(w*N))) + '}$' for w in xticks])
      xtlbls[xticks==1/N] = '1'
      ax.set_xticks(xticks)
      ax.set_xticklabels(xtlbls)
      self.fig  = fig
      self.ax   = ax
      self.hist = hist
    else:
      not_available_text(ax,'Not computed (N > threshold)')

    #TODO:
    # self.update(0)

  def update(self,w):
    ax        = self.ax
    _         = [b.remove() for b in self.hist]
    N         = len(w)
    wmax      = w.max()
    bins      = exp(linspace(log(1e-5/N), log(1), int(N/20)))
    counted   = w>bins[0]
    nC        = np.sum(counted)
    nn,_,pp   = ax.hist(w[counted], bins=bins, color='b')
    self.hist = pp
    #thresh   = '#(w<$10^{'+ str(int(log10(bins[0]*N))) + '}/N$ )'
    ax.set_title('N: {:d}.   N_eff: {:.4g}.   Not shown: {:d}. '.\
        format(N, 1/(w@w), N-nC))
    ax.set_ylim(*d_ylim([nn]))


class spectral_errors_LP:

  def __init__(self,fignum,stats):
    fig, ax = freshfig(fignum, (6,3))
    win_title(fig,"Spectral view")
    ax.set_xlabel('Sing. value index')
    ax.set_yscale('log')
    ax.set_ylim(bottom=1e-5)
    #ax.set_ylim([1e-3,1e1])

    try:
      self.msft = stats.umisf
      self.sprd = stats.svals
    except KeyError:
      do_spectral_error = False
      not_available_text(ax, "Spectral stats not being computed")
    else:
      if np.any(np.isinf(self.msft[0])):
        not_available_text(ax, "Spectral stats not finite")
        do_spectral_error = False
      else:
        do_spectral_error = True

    if do_spectral_error:
      M = len(self.msft[0])
      self.line_msft, = ax.plot(arange(M),ones(M),'k',lw=2,label='Error')
      self.line_sprd, = ax.plot(arange(M),ones(M),'b',lw=2,label='Spread',alpha=0.9)
      ax.get_xaxis().set_major_locator(MaxNLocator(integer=True))
      ax.legend()

    self.do_spectral_error = do_spectral_error
    self.fig = fig
    self.ax  = ax
    self.update(0)

  def update(self,k):
    msft = abs(self.msft[k])
    sprd =     self.sprd[k]
    self.line_sprd.set_ydata(sprd)
    self.line_msft.set_ydata(msft)
    self.ax.set_ylim(*d_ylim(msft))


class correlations_LP:

  # Whether to show half/full (symmetric) corr matrix.
  half = True

  def __init__(self,fignum,Nx,E,P):
    GS = {'height_ratios':[4, 1],'hspace':0.09,'top':0.95}
    fig, (ax,ax2) = freshfig(fignum, (5,6), nrows=2, gridspec_kw=GS)
    win_title(fig, "Correlations")

    if Nx<=1003:
      C = eye(Nx)
      # Mask half
      mask = np.zeros_like(C, dtype=np.bool)
      mask[np.tril_indices_from(mask)] = True
      # Make colormap. Log-transform cmap, but not internally in matplotlib,
      # so as to avoid transforming the colorbar too.
      cmap = plt.get_cmap('RdBu')
      trfm = colors.SymLogNorm(linthresh=0.2,linscale=0.2,vmin=-1, vmax=1)
      cmap = cmap(trfm(linspace(-0.6,0.6,cmap.N)))
      cmap = colors.ListedColormap(cmap)
      #
      VM   = 1.0 # abs(np.percentile(C,[1,99])).max()
      im   = ax.imshow(C,cmap=cmap,vmin=-VM,vmax=VM)
      # Colorbar
      cax = ax.figure.colorbar(im,ax=ax,shrink=0.8)
      # Tune plot
      plt.box(False)
      ax.set_facecolor('w') 
      ax.grid(False)
      ax.set_title("State correlation matrix:", y=1.07)
      ax.xaxis.tick_top()
      
      # ax2 = inset_axes(ax,width="30%",height="60%",loc=3)
      line_AC, = ax2.plot(arange(Nx), ones(Nx), label='Correlation')
      line_AA, = ax2.plot(arange(Nx), ones(Nx), label='Abs. corr.')
      _        = ax2.hlines(0,0,Nx-1,'k','dotted',lw=1)
      # Align ax2 with ax
      bb_AC = ax2.get_position()
      bb_C  = ax.get_position()
      ax2.set_position([bb_C.x0, bb_AC.y0, bb_C.width, bb_AC.height])
      # Tune plot
      ax2.set_title("Auto-correlation:")
      ax2.set_ylabel("Mean value")
      ax2.set_xlabel("Distance (in state indices)")
      ax2.set_xticklabels([])
      ax2.set_yticks([0,1] + list(ax2.get_yticks()[[0,-1]]))
      ax2.set_ylim(top=1)
      ax2.legend(frameon=True,facecolor='w',
          bbox_to_anchor=(1, 1), loc='upper left', borderaxespad=0.02)

      self.fig     = fig
      self.ax      = ax
      self.ax2     = ax2
      self.im      = im
      self.line_AC = line_AC
      self.line_AA = line_AA
      self.mask    = mask

      self.update(E,P)
    else:
      not_available_text(ax)

  def update(self,E,P):
    # Get cov matrix
    if E is not None:
      C = np.cov(E,rowvar=False)
    else:
      assert P is not None
      C = P.full if isinstance(P,CovMat) else P
      C = C.copy()
    # Compute corr from cov
    std = sqrt(diag(C))
    C  /= std[:,None]
    C  /= std[None,:]
    # Mask
    if self.half:
      C = np.ma.masked_where(self.mask, C)
    # Plot
    self.im.set_data(C)
    # Auto-corr function
    ACF = circulant_ACF(C)
    AAF = circulant_ACF(C,do_abs=True)
    self.line_AC.set_ydata(ACF)
    self.line_AA.set_ydata(AAF)


def plot_pause(duration):
  """
  plt.pause is not supported by jupyter notebook.
  Provide fallback that does work.
  stackoverflow.com/q/34486642
  """
  try:
    plt.pause(duration)
  except:
    fig = plt.gcf()
    fig.canvas.draw()
    time.sleep(0.1)


def setup_wrapping(M,periodic=True):
  """
  Make periodic indices and a corresponding function
  (that works for ensemble input).
  """

  if periodic:
    ii = np.hstack([-0.5, arange(M), M-0.5])
    def wrap(E):
      midpoint = (E[[0],...] + E[[-1],...])/2
      return ccat(midpoint,E,midpoint)

  else:
    ii = arange(M)
    wrap = lambda x: x

  return ii, wrap
  
def adjust_position(ax,adjust_extent=False,**kwargs):
  """
  Adjust values (add) to get_position().
  kwarg must be one of 'x0','y0','width','height'.
  """
  # Load get_position into d
  pos = ax.get_position()
  d   = OrderedDict()
  for key in ['x0','y0','width','height']:
    d[key] = getattr(pos,key)
  # Make adjustments
  for key,item in kwargs.items():
    d[key] += item
    if adjust_extent:
      if key=='x0': d['width']  -= item
      if key=='y0': d['height'] -= item
  # Set
  ax.set_position(d.values())

def span(xx,axis=None):
  a = xx.min(axis)
  b = xx.max(axis)
  return a, b

def stretch(a,b,factor=1,int=False):
  """
  Stretch distance a-b by factor.
  Return a,b.
  If int: floor(a) and ceil(b)
  """
  c = (a+b)/2
  a = c + factor*(a-c) 
  b = c + factor*(b-c) 
  if int:
    a = floor(a)
    b = ceil(b)
  return a, b



def d_ylim(data,ax=None,cC=0,cE=1,pp=(1,99),Min=-1e20,Max=+1e20):
  """
  Provide new ylim's intelligently,
  computed from percentiles of the data.
  - data: iterable of arrays for computing percentiles.
  - pp: percentiles

  - ax: If present, then the delta_zoom in/out is also considered.
    - cE: exansion (widenting) rate ∈ [0,1].
        Default: 1, which immediately expands to percentile.
    - cC: compression (narrowing) rate ∈ [0,1].
        Default: 0, which does not allow compression.
  
  - Min/Max: bounds

  Despite being a little involved,
  the cost of this subroutine is typically not substantial
  because there's usually not that much data to sort through.
  """

  # Find "reasonable" limits (by percentiles), looping over data
  maxv = minv = -np.inf # init
  for d in data:
    d = d[np.isfinite(d)]
    if len(d):
      minv, maxv = np.maximum([minv, maxv], \
          array([-1, 1]) * np.percentile(d,pp))
  minv *= -1

  # Pry apart equal values
  if np.isclose(minv,maxv):
    maxv += 0.5
    minv -= 0.5

  # Make the zooming transition smooth
  if ax is not None:
    current = ax.get_ylim()
    # Set rate factor as compress or expand factor. 
    c0 = cC if minv>current[0] else cE
    c1 = cC if maxv<current[1] else cE
    # Adjust
    minv = np.interp(c0, (0,1), (current[0], minv))
    maxv = np.interp(c1, (0,1), (current[1], maxv))

  # Bounds 
  maxv = min(Max,maxv)
  minv = max(Min,minv)

  # Set (if anything's changed)
  def worth_updating(a,b,curr):
    # Note: should depend on cC and cE
    d = abs(curr[1]-curr[0])
    lower = abs(a-curr[0]) > 0.002*d
    upper = abs(b-curr[1]) > 0.002*d
    return lower and upper
  #if worth_updating(minv,maxv,current):
    #ax.set_ylim(minv,maxv)

  return minv, maxv


def set_ilim(ax,i,Min=None,Max=None):
  """Set bounds on axis i.""" 
  if i is 0: ax.set_xlim(Min,Max)
  if i is 1: ax.set_ylim(Min,Max)
  if i is 2: ax.set_zlim(Min,Max)

def estimate_good_plot_length(xx,chrono=None,mult=100):
  """
  Estimate good length for plotting stuff
  from the time scale of the system.
  Provide sensible fall-backs (better if chrono is supplied).
  """
  if xx.ndim == 2:
    # If mult-dim, then average over dims (by ravel)....
    # But for inhomogeneous variables, it is important
    # to subtract the mean first!
    xx = xx - mean(xx,axis=0)
    xx = xx.ravel(order='F')

  try:
    K = mult * estimate_corr_length(xx)
  except ValueError:
    K = 0

  if chrono != None:
    t = chrono
    K = int(min(max(K, t.dkObs), t.K))
    T = round2sigfig(t.tt[K],2) # Could return T; T>tt[-1]
    K = find_1st_ind(t.tt >= T)
    if K: return K
    else: return t.K
  else:
    K = int(min(max(K, 1), len(xx)))
    T = round2sigfig(K,2)
    return K

def get_plot_inds(xx,chrono,K=None,T=None,**kwargs):
  """
  Def subset of kk for plotting, from one of
   - K
   - T
   - mult * auto-correlation length of xx
  """
  t = chrono
  if K is None:
    if T: K = find_1st_ind(t.tt >= min(T,t.T))
    else: K = estimate_good_plot_length(xx,chrono=t,**kwargs)
  plot_kk    = t.kk[:K+1]
  plot_kkObs = t.kkObs[t.kkObs<=K]
  return plot_kk, plot_kkObs


def plot_3D_trajectory(stats,dims=0,**kwargs):
  """
  Plot 3D phase-space trajectory.
  kwargs forwarded to get_plot_inds().
  """
  if is_int(dims):
    dims = dims + arange(3)
  assert len(dims)==3

  xx     = stats.xx
  mu     = stats.mu
  chrono = stats.HMM.t

  kk,kkA = get_plot_inds(xx,chrono,mult=100,**kwargs)

  if mu.store_u:
    xx = xx[kk]
    mu = mu[kk]
    T  = chrono.tt[kk[-1]]
  else:
    xx = xx[kkA]
    mu = mu.a[:len(kkA)]
    T  = chrono.tt[kkA[-1]]

  plt.figure(14).clf()
  set_figpos('3321 mac')
  ax3 = plt.subplot(111, projection='3d')

  xx = xx.T[dims]
  mu = mu.T[dims]

  ax3.plot   (*xx      ,c='k',label='Truth')
  ax3.plot   (*mu      ,c='b',label='DA estim.')
  ax3.scatter(*xx[:, 0],c='g',s=40)
  ax3.scatter(*xx[:,-1],c='r',s=40)

  ax3.set_title('Phase space trajectory up to t={:<5.2f}'.format(T))
  ax3.set_xlabel('dim ' + str(dims[0]))
  ax3.set_ylabel('dim ' + str(dims[1]))
  ax3.set_zlabel('dim ' + str(dims[2]))
  ax3.legend(frameon=False)
  ax3.set_facecolor('w')


def plot_time_series(stats,**kwargs):
  """
  Plot time series of various statistics.
  kwargs forwarded to get_plot_inds().
  """

  # Figure, axes
  fg = plt.figure(12,figsize=(5,3.5))
  fg.clf()
  set_figpos('1313 mac')
  fg, (ax_e,ax_K) = plt.subplots(2,1,sharex=True,num=12)

  # Time
  chrono = stats.HMM.t
  xx     = stats.xx
  Nx     = xx.shape[1]
  dims   = equi_spaced_integers(Nx, min(Nx, 10))
  kk,kkA = get_plot_inds(xx[:,dims],chrono,mult=80,**kwargs)
  tt,ttA = chrono.tt[kk], chrono.tt[kkA]
  KA     = len(kkA)      

  # Stats
  s = stats
  if s.mu.store_u:
    tt_  = tt 
    rmse = s.rmse[kk]
    rmv  = s.rmv [kk]
  else:
    tt_  = ttA
    rmse = s.rmse.a[:KA]
    rmv  = s.rmv .a[:KA]

  trKH   = s.trHK  [:KA]
  skew   = s.skew.a[:KA]
  kurt   = s.kurt.a[:KA]

  ax_e.plot(        tt_, rmse,'k',lw=2 ,label='Error')
  ax_e.fill_between(tt_, rmv ,alpha=0.7,label='Spread') 
  ax_e.set_ylim(0, 1.1*max(np.percentile(rmse,99), rmv.max()) )
  ax_e.set_ylabel('RMS')
  ax_e.legend()

  ax_K.plot(ttA, trKH,'k',lw=2,label='HK')
  ax_K.plot(ttA, skew,'g',lw=2,label='Skew')
  ax_K.plot(ttA, kurt,'r',lw=2,label='Kurt')
  ax_K.set_xlabel('Time (t)')
  ax_K.set_ylabel('Mean of marginal\n $\sigma$-normalized values',
      fontsize='small', labelpad=0)
  ax_K.legend()


def plot_hovmoller(xx,chrono=None,**kwargs):
  """
  Plot Hovmöller diagram.
  kwargs forwarded to get_plot_inds().
  """
  #cm = mpl.colors.ListedColormap(sns.color_palette("BrBG", 256)) # RdBu_r
  #cm = plt.get_cmap('BrBG')
  fig, ax = freshfig(16,(4,3.5))
  set_figpos('3311 mac')

  Nx = xx.shape[1]

  if chrono!=None:
    kk,_ = get_plot_inds(xx,chrono,mult=40,**kwargs)
    tt   = chrono.tt[kk]
    ax.set_ylabel('Time (t)')
  else:
    K    = estimate_good_plot_length(xx,mult=40)
    tt   = arange(K)
    ax.set_ylabel('Time indices (k)')

  plt.contourf(arange(Nx),tt,xx[kk],25)
  plt.colorbar()
  ax.set_position([0.125, 0.20, 0.62, 0.70])
  ax.set_title("Hovmoller diagram (of 'Truth')")
  ax.set_xlabel('Dimension index (i)')
  add_endpoint_xtick(ax)


def add_endpoint_xtick(ax):
  """Useful when xlim(right) is e.g. 39 (instead of 40)."""
  xF = ax.get_xlim()[1]
  ticks = ax.get_xticks()
  if ticks[-1] > xF:
    ticks = ticks[:-1]
  ticks = np.append(ticks, xF)
  ax.set_xticks(ticks)


def integer_hist(E,N,centrd=False,weights=None,**kwargs):
  """Histogram for integers."""
  ax = plt.gca()
  rnge = (-0.5,N+0.5) if centrd else (0,N+1)
  ax.hist(E,bins=N+1,range=rnge,normed=1,weights=weights,**kwargs)
  ax.set_xlim(rnge)


def not_available_text(ax,txt=None,fs=20):
  if txt is None: txt = '[Not available]'
  else:           txt = '[' + txt + ']'
  ax.text(0.5,0.5,txt,
      fontsize=fs,
      transform=ax.transAxes,
      va='center',ha='center',
      wrap=True)

def plot_err_components(stats):
  """
  Plot components of the error.
  Note: it was chosen to plot(ii, mean_in_time(abs(err_i))),
        and thus the corresponding spread measure is MAD.
        If one chose instead: plot(ii, std_in_time(err_i)),
        then the corresponding measure of spread would have been std.
        This choice was made in part because (wrt. subplot 2)
        the singular values (svals) correspond to rotated MADs,
        and because rms(umisf) seems to convoluted for interpretation.
  """
  fgE = plt.figure(15,figsize=(6,6)).clf()
  set_figpos('1312 mac')

  chrono = stats.HMM.t
  Nx     = stats.xx.shape[1]

  err   = mean( abs(stats.err  .a) ,0)
  sprd  = mean(     stats.mad  .a  ,0)
  umsft = mean( abs(stats.umisf.a) ,0)
  usprd = mean(     stats.svals.a  ,0)

  ax_r = plt.subplot(311)
  ax_r.plot(          arange(Nx),               err,'k',lw=2, label='Error')
  if Nx<10**3:
    ax_r.fill_between(arange(Nx),[0]*len(sprd),sprd,alpha=0.7,label='Spread')
  else:
    ax_r.plot(        arange(Nx),              sprd,alpha=0.7,label='Spread')
  #ax_r.set_yscale('log')
  ax_r.set_title('Element-wise error comparison')
  ax_r.set_xlabel('Dimension index (i)')
  ax_r.set_ylabel('Time-average (_a) magnitude')
  ax_r.set_ylim(bottom=mean(sprd)/10)
  ax_r.set_xlim(right=Nx-1); add_endpoint_xtick(ax_r)
  ax_r.get_xaxis().set_major_locator(MaxNLocator(integer=True))
  plt.subplots_adjust(hspace=0.55) # OR: [0.125,0.6, 0.78, 0.34]
  ax_r.legend()

  ax_s = plt.subplot(312)
  ax_s.set_xlabel('Principal component index')
  ax_s.set_ylabel('Time-average (_a) magnitude')
  ax_s.set_title('Spectral error comparison')
  has_been_computed = np.any(np.isfinite(umsft))
  if has_been_computed:
    L = len(umsft)
    ax_s.plot(        arange(L),      umsft,'k',lw=2, label='Error')
    ax_s.fill_between(arange(L),[0]*L,usprd,alpha=0.7,label='Spread')
    ax_s.set_yscale('log')
    ax_s.set_ylim(bottom=1e-4*usprd.sum())
    ax_s.set_xlim(right=Nx-1); add_endpoint_xtick(ax_s)
    ax_s.get_xaxis().set_major_locator(MaxNLocator(integer=True))
    ax_s.legend()
  else:
    not_available_text(ax_s)

  rmse = stats.rmse.a[chrono.maskObs_BI]
  ax_R = plt.subplot(313)
  ax_R.hist(rmse,bins=30,normed=0)
  ax_R.set_ylabel('Num. of occurence (_a)')
  ax_R.set_xlabel('RMSE')
  ax_R.set_title('Histogram of RMSE values')


def plot_rank_histogram(stats):
  chrono = stats.HMM.t

  has_been_computed = \
      hasattr(stats,'rh') and \
      not all(stats.rh.a[-1]==array(np.nan).astype(int))

  def are_uniform(w):
    """Test inital & final weights, not intermediate (for speed)."""
    (w[0]==1/N).all() and (w[-1]==1/N).all()

  fg = plt.figure(13,figsize=(6,3)).clf()
  set_figpos('3331 mac')
  #
  ax_H = plt.subplot(111)
  ax_H.set_title('(Average of marginal) rank histogram (_a)')
  ax_H.set_ylabel('Freq. of occurence\n (of truth in interval n)')
  ax_H.set_xlabel('ensemble member index (n)')
  ax_H.set_position([0.125,0.15, 0.78, 0.75])
  if has_been_computed:
    w     = stats.w.a [chrono.maskObs_BI]
    ranks = stats.rh.a[chrono.maskObs_BI]
    Nx    = ranks.shape[1]
    N     = w.shape[1]
    if are_uniform(w):
      # Ensemble rank histogram
      integer_hist(ranks.ravel(),N)
    else:
      # Experimental: weighted rank histogram.
      # Weight ranks by inverse of particle weight. Why? Coz, with correct
      # importance weights, the "expected value" histogram is then flat.
      # Potential improvement: interpolate weights between particles.
      w  = w
      K  = len(w)
      w  = np.hstack([w, ones((K,1))/N]) # define weights for rank N+1
      w  = array([ w[arange(K),ranks[arange(K),i]] for i in range(Nx)])
      w  = w.T.ravel()
      w  = np.maximum(w, 1/N/100) # Artificial cap. Reduces variance, but introduces bias.
      w  = 1/w
      integer_hist(ranks.ravel(),N,weights=w)
  else:
    not_available_text(ax_H)
  

def adjustable_box_or_forced():
  "For set_aspect(), adjustable='box-forced' replaced by 'box' since mpl 2.2.0."
  from pkg_resources import parse_version as pv
  return 'box-forced' if pv(mpl.__version__) < pv("2.2.0") else 'box'


def show_figs(fignums=None):
  """Move all fig windows to top"""
  if fignums == None:
    fignums = plt.get_fignums()
  try:
    fignums = list(fignums)
  except:
    fignums = [fignums]
  for f in fignums:
    plt.figure(f)
    fmw = plt.get_current_fig_manager().window
    fmw.attributes('-topmost',1) # Bring to front, but
    fmw.attributes('-topmost',0) # don't keep in front
  
def win_title(fig, string, num=True):
  "Set window title"
  if num:
    N      = fig.number
    string += " [" + str(N) + "]"
  fig.canvas.set_window_title(string)

def set_figpos(loc):
  """
  Place figure on screen, where 'loc' can be either
    NW, E, ...
  or
    4 digits (as str or int) to define grid M,N,i,j.
  """

  #Only works with both:
   #- Patrick's monitor setup (Dell with Mac central-below)
   #- TkAgg backend. (Previously: Qt4Agg)
  if not user_is_patrick or mpl.get_backend() != 'TkAgg':
    return
  fmw = plt.get_current_fig_manager().window

  loc = str(loc)

  # Qt4Agg only:
  #  # Current values 
  #  w_now = fmw.width()
  #  h_now = fmw.height()
  #  x_now = fmw.x()
  #  y_now = fmw.y()
  #  # Constants 
  #  Dell_w = 2560
  #  Dell_h = 1440
  #  Mac_w  = 2560
  #  Mac_h  = 1600
  #  # Why is Mac monitor scaled by 1/2 ?
  #  Mac_w  /= 2
  #  Mac_h  /= 2
  # Append the string 'mac' to place on mac monitor.
  #  if 'mac' in loc:
  #    x0 = Dell_w/4
  #    y0 = Dell_h+44
  #    w0 = Mac_w
  #    h0 = Mac_h-44
  #  else:
  #    x0 = 0
  #    y0 = 0
  #    w0 = Dell_w
  #    h0 = Dell_h

  # TkAgg
  x0 = 0
  y0 = 0
  w0 = 1280
  h0 = 752
  
  # Def place function with offsets
  def place(x,y,w,h):
    #fmw.setGeometry(x0+x,y0+y,w,h) # For Qt4Agg
    geo = str(int(w)) + 'x' + str(int(h)) + \
        '+' + str(int(x)) + '+' + str(int(y))
    fmw.geometry(newGeometry=geo) # For TkAgg

  if not loc[:4].isnumeric():
    if   loc.startswith('NW'): loc = '2211'
    elif loc.startswith('SW'): loc = '2221'
    elif loc.startswith('NE'): loc = '2212'
    elif loc.startswith('SE'): loc = '2222'
    elif loc.startswith('W' ): loc = '1211'
    elif loc.startswith('E' ): loc = '1212'
    elif loc.startswith('S' ): loc = '2121'
    elif loc.startswith('N' ): loc = '2111'

  # Place
  M,N,i,j = [int(x) for x in loc[:4]]
  assert M>=i>0 and N>=j>0
  h0   -= (M-1)*25
  yoff  = 25*(i-1)
  if i>1:
    yoff += 25
  place((j-1)*w0/N, yoff + (i-1)*h0/M, w0/N, h0/M)


# stackoverflow.com/a/7396313
from matplotlib import transforms as mtransforms
def autoscale_based_on(ax, line_handles):
  "Autoscale axis based (only) on line_handles."
  ax.dataLim = mtransforms.Bbox.unit()
  for iL,lh in enumerate(line_handles):
    xy = np.vstack(lh.get_data()).T
    ax.dataLim.update_from_data_xy(xy, ignore=(iL==0))
  ax.autoscale_view()


from matplotlib.widgets import CheckButtons
import textwrap
def toggle_lines(ax=None,autoscl=True,numbering=False,txtwidth=15,txtsize=None,state=None):
  """
  Make checkbuttons to toggle visibility of each line in current plot.
  autoscl  : Rescale axis limits as required by currently visible lines.
  numbering: Add numbering to labels.
  txtwidth : Wrap labels to this length.

  State of checkboxes can be inquired by 
  OnOff = [lh.get_visible() for lh in ax.findobj(lambda x: isinstance(x,mpl.lines.Line2D))[::2]]
  """

  if ax is None: ax = plt.gca()
  if txtsize is None: txtsize = mpl.rcParams['font.size']

  # Get lines and their properties
  lines = {'handle': list(ax.get_lines())}
  for prop in ['label','color','visible']:
    lines[prop] = [plt.getp(x,prop) for x in lines['handle']]
  # Put into pandas for some reason
  lines = pd.DataFrame(lines)
  # Rm those that start with _
  lines = lines[~lines.label.str.startswith('_')]

  # Adjust labels
  if numbering: lines['label'] = [str(i)+': '+lbl for i,lbl in enumerate(lines['label'])]
  if txtwidth:  lines['label'] = [textwrap.fill(lbl,width=txtwidth) for lbl in lines['label']]

  # Set state. BUGGY? sometimes causes MPL complaints after clicking boxes
  if state is not None:
    state = array(state).astype(bool)
    lines.visible = state
    for i,x in enumerate(state):
      lines['handle'][i].set_visible(x)

  # Setup buttons
  # When there's many, the box-sizing is awful, but difficult to fix.
  W       = 0.23 * txtwidth/15 * txtsize/10
  N       = len(lines)
  nBreaks = sum(lbl.count('\n') for lbl in lines['label']) # count linebreaks
  H       = min(1,0.05*(N+nBreaks))
  plt.subplots_adjust(left=W+0.12,right=0.97)
  rax = plt.axes([0.05, 0.5-H/2, W, H])
  check = CheckButtons(rax, lines.label, lines.visible)

  # Adjust button style
  for i in range(N):
    check.rectangles[i].set(lw=0,facecolor=lines.color[i])
    check.labels[i].set(color=lines.color[i])
    if txtsize: check.labels[i].set(size=txtsize)

  # Callback
  def toggle_visible(label):
    ind = lines.label==label
    handle = lines[ind].handle.item()
    vs = not lines[ind].visible.item()
    handle.set_visible( vs )
    lines.loc[ind,'visible'] = vs
    if autoscl:
      autoscale_based_on(ax,lines[lines.visible].handle)
    plt.draw()
  check.on_clicked(toggle_visible)

  # Return focus
  plt.sca(ax)

  # Must return (and be received) so as not to expire.
  return check


@vectorize0
def toggle_viz(h,prompt=False,legend=False):
  """Toggle visibility of the graphics with handle h."""

  # Core functionality: turn on/off
  is_viz = not h.get_visible()
  h.set_visible(is_viz)

  if prompt:
    input("Press <Enter> to continue...")

  # Legend updating. Basic version: works by
  #  - setting line's label to actual_label/'_nolegend_' if is_viz/not
  #  - re-calling legend()
  if legend:
      if is_viz:
        try:
          h.set_label(h.actual_label)
        except AttributeError:
          pass
      else:
        h.actual_label = h.get_label()
        h.set_label('_nolegend_')
      # Legend refresh
      ax = h.axes
      with warnings.catch_warnings():
        warnings.simplefilter("error",category=UserWarning)
        try:
          ax.legend()
        except UserWarning:
          # If all labels are '_nolabel_' then ax.legend() throws warning,
          # and quits before refreshing. => Refresh by creating/rm another legend.
          ax.legend('TMP').remove()

  plt.pause(0.02)
  return is_viz


def freshfig(num=None,figsize=None,*args,**kwargs):
  """
  - If the figure does not exist: create figure it.
    This allows for figure sizing -- even on Macs.
  - Otherwise: clear figure (we avoid closing/opening so as
    to keep (potentially manually set) figure positioning.
  - The rest is the same as:
  >>> fig, ax = suplots()
  """
  fig = plt.figure(num=num,figsize=figsize)
  fig.clf()
  _, ax = plt.subplots(num=fig.number,*args,**kwargs)
  return fig, ax

def savefig_n(f=None):
  """
  Simplify the exporting of a figure, especially when it's part of a series.
  """
  assert savefig_n.index>=0, "Initalize using savefig_n.index = 1 in your script"
  if f is None:
    f = inspect.getfile(inspect.stack()[1][0])   # Get __file__ of caller
    f = save_dir(f)                              # Prep save dir
  f = f + str(savefig_n.index) + '.pdf'          # Compose name
  print("Saving fig to:",f)                      # Print
  plt.savefig(f)                                 # Save
  savefig_n.index += 1                           # Increment index
  plt.pause(0.1)                                 # For safety?
savefig_n.index = -1


from matplotlib.gridspec import GridSpec
def axes_with_marginals(n_joint, n_marg,**kwargs):
  """
  Create a joint axis along with two marginal axes.

  Example:
  >>> ax_s, ax_x, ax_y = axes_with_marginals(4, 1)
  >>> x, y = np.random.randn(2,500)
  >>> ax_s.scatter(x,y)
  >>> ax_x.hist(x)
  >>> ax_y.hist(y,orientation="horizontal")
  """

  N = n_joint + n_marg

  # Method 1
  #fig, ((ax_s, ax_y), (ax_x, _)) = plt.subplots(2,2,num=plt.gcf().number,
      #sharex='col',sharey='row',gridspec_kw={
        #'height_ratios':[n_joint,n_marg],
        #'width_ratios' :[n_joint,n_marg]})
  #_.set_visible(False) # Actually removing would bug the axis ticks etc.
  
  # Method 2
  gs   = GridSpec(N,N,**kwargs)
  fig  = plt.gcf()
  ax_s = fig.add_subplot(gs[n_marg:N     ,0      :n_joint])
  ax_x = fig.add_subplot(gs[0     :n_marg,0      :n_joint],sharex=ax_s)
  ax_y = fig.add_subplot(gs[n_marg:N     ,n_joint:N      ],sharey=ax_s)
  # Cannot delete ticks coz axis are shared
  plt.setp(ax_x.get_xticklabels(), visible=False)
  plt.setp(ax_y.get_yticklabels(), visible=False)

  return ax_s, ax_x, ax_y

from matplotlib.patches import Ellipse
def cov_ellipse(ax, mu, sigma, **kwargs):
    """
    Draw ellipse corresponding to (Gaussian) 1-sigma countour of cov matrix.

    Inspired by stackoverflow.com/q/17952171

    Example:
    >>> ellipse = cov_ellipse(ax, y, R,
    >>>           facecolor='none', edgecolor='y',lw=4,label='$1\\sigma$')
    """

    # Cov --> Width, Height, Theta
    vals, vecs = np.linalg.eigh(sigma)
    x, y       = vecs[:, -1] # x-y components of largest (last) eigenvector
    theta      = np.degrees(np.arctan2(y, x))
    theta      = theta % 180

    h, w       = 2 * np.sqrt(vals.clip(0))

    # Get artist
    e = Ellipse(mu, w, h, theta, **kwargs)

    ax.add_patch(e)
    e.set_clip_box(ax.bbox) # why is this necessary?

    # Return artist
    return e
    



