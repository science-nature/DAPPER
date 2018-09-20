from common import *

#from mpl_toolkits.mplot3d import Axes3D
from mpl_toolkits.mplot3d.art3d import juggle_axes

from mpl_toolkits.axes_grid1 import make_axes_locatable
from mpl_toolkits.axes_grid1.inset_locator import inset_axes

from matplotlib import colors
from matplotlib.ticker import MaxNLocator

tri_corr = True

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

    setup  = stats.setup
    config = stats.config
    m      = setup.f.m
    dt     = setup.t.dt

    # Store
    self.setup = setup
    self.stats = stats
    self.xx    = stats.xx ; xx = stats.xx
    self.yy    = stats.yy ; yy = stats.yy

    # Abbreviate
    mu = stats.mu

    # Set up prompts
    self.is_on     = True
    self.is_paused = False
    print('Initializing liveplotting...')
    print('Press <Enter> to toggle live plot OFF/ON.')
    print('Press <Space> and then <Enter> to pause.')

    #ens_props = {} yields rainbow
    ens_props = {'color': 0.7*RGBs['w'],'alpha':0.3}

    # For periodic functions
    ii,wrap = setup_wrapping(m)


    #####################
    # Correlation plot
    #####################
    if 2 in only and m<1001:
      GS = {'height_ratios':[4, 1],'hspace':0.09,'top':0.95}
      fig_C, (ax_C,ax_AC) = freshfig(2, (5,6), nrows=2, gridspec_kw=GS)
      win_title(fig_C, "Correlations")
      set_figpos('2311')

      if m<=1003:
        # Get cov matrix
        if E is not None:
          C = np.cov(E.T, ddof=1)
        else:
          assert P is not None
          C = P.full if isinstance(P,CovMat) else P
          C = C.copy()
        # Compute corr from cov
        std = sqrt(diag(C))
        C  /= std[:,None]
        C  /= std[None,:]
        # Mask half
        mask = np.zeros_like(C, dtype=np.bool)
        mask   [np.tril_indices_from(mask)] = True
        if tri_corr:
          C  = np.ma.masked_where(mask, C)
        # Make colormap. Log-transform cmap, but not internally in matplotlib,
        # so as to avoid transforming the colorbar too.
        cmap = plt.get_cmap('RdBu')
        trfm = colors.SymLogNorm(linthresh=0.2,linscale=0.2,vmin=-1, vmax=1)
        cmap = cmap(trfm(linspace(-0.6,0.6,cmap.N)))
        cmap = colors.ListedColormap(cmap)
        #
        VM   = 1.0 # abs(np.percentile(C,[1,99])).max()
        im_C = ax_C.imshow(C,cmap=cmap,vmin=-VM,vmax=VM)
        #
        cax = ax_C.figure.colorbar(im_C,ax=ax_C,shrink=0.8)
        plt.box(False)
        ax_C.set_facecolor('w') 
        ax_C.grid(False)
        ax_C.set_title("State correlation matrix:", y=1.07)
        ax_C.xaxis.tick_top()
        
        # ax_AC = inset_axes(ax_C,width="30%",height="60%",loc=3)
        ACF = circulant_ACF(C)
        AAF = circulant_ACF(C,do_abs=True)
        line_AC, = ax_AC.plot(arange(m), ACF, label='Correlation')
        line_AA, = ax_AC.plot(arange(m), AAF, label='Abs. corr.')
        _   = ax_AC.hlines(0,0,m-1,'k','dotted',lw=1)
        # Align ax_AC with ax_C
        bb_AC = ax_AC.get_position()
        bb_C  = ax_C.get_position()
        ax_AC.set_position([bb_C.x0, bb_AC.y0, bb_C.width, bb_AC.height])
        # Tune plot
        ax_AC.set_title("Auto-correlation:")
        ax_AC.set_ylabel("Mean value")
        ax_AC.set_xlabel("Distance (in state indices)")
        ax_AC.set_xticklabels([])
        ax_AC.set_yticks([0,1] + list(ax_AC.get_yticks()[[0,-1]]))
        ax_AC.set_ylim(top=1)
        ax_AC.legend(frameon=True,facecolor='w',
            bbox_to_anchor=(1, 1), loc='upper left', borderaxespad=0.02)

        self.fig_C   = fig_C
        self.ax_C    = ax_C
        self.ax_AC   = ax_AC
        self.im_C    = im_C
        self.line_AC = line_AC
        self.line_AA = line_AA
        self.mask    = mask
      else:
        not_available_text(ax_C)
      

    #####################
    # Spectral error plot
    #####################
      fig_S, ax_S = freshfig(4, (4,3))
      win_title(fig_S, "Spectral view")
      set_figpos('2311')
      ax_S.set_xlabel('Sing. value index')
      ax_S.set_yscale('log')
      ax_S.set_ylim(bottom=1e-5)
      #ax_S.set_ylim([1e-3,1e1])
      try:
        msft = abs(stats.umisf[0])
        sprd =     stats.svals[0]
      except KeyError:
        self.do_spectral_error = False
        not_available_text(ax_S, "Spectral stats not being computed")
      else:
        if np.any(np.isfinite(msft)):
          not_available_text(ax_S, "Spectral stats not finite")
          self.do_spectral_error = False
        else:
          self.do_spectral_error = True

      if self.do_spectral_error:
        self.line_msft, = ax_S.plot(arange(len(msft)),msft,'k',lw=2,label='Error')
        self.line_sprd, = ax_S.plot(arange(len(sprd)),sprd,'b',lw=2,label='Spread',alpha=0.9)
        ax_S.get_xaxis().set_major_locator(MaxNLocator(integer=True))
        ax_S.legend()

      self.ax_S = ax_S


    #####################
    # Diagnostics
    #####################
    if 1 in only:
      GS = {'left':0.125+0.04,'right':0.9+0.04}
      self.fig_pulse, (self.ax_RMS, self.ax_Uni) = freshfig(
          1,(5,3.5),nrows=2,sharex=True,gridspec_kw=GS)
      win_title(self.fig_pulse,"Scalar diagnostics")
      set_figpos('2312')
      self.has_checked_presence = False

      def lin(a,b):
        def f(x):
          y = a + b*x
          return y
        return f

      def divN():
        try:
          N = E.shape[0]
          def f(x): return x/N
          return f
        except AttributeError:
          pass
        

      # --------------
      # RMS
      # --------------
      d_RMS = {
          'rmse' : dict(c='k', label='Error'),
          'rmv'  : dict(c='b', label='Spread', alpha=0.6),
        }

      # --------------
      # Plain
      # --------------
      d_Uni = OrderedDict([
          ('skew'   , dict(c='g', label='Skew')),
          ('kurt'   , dict(c='r', label='Kurt')),
          ('infl'   , dict(c='c', label='(Infl-1)*10',transf=lin(-10,10))),
          ('N_eff'  , dict(c='y', label='N_eff/N'    ,transf=divN(),    step=True)),
          ('iters'  , dict(c='m', label='Iters/2'    ,transf=lin(0,.5), step=True)),
          ('trHK'   , dict(c='k', label='HK')),
          ('resmpl' , dict(c='k', label='Resampl?')),
        ])

      chrono       = setup.t
      chrono.pK    = estimate_good_plot_length(xx,chrono,mult=80)
      chrono.pKObs = int(chrono.pK / chrono.dkObs)

      def raise_field_lvl(dct,fld):
        dct[fld] = dct['plt'][fld]
        del dct['plt'][fld]

      def init_axd(ax,dict_of_dicts):
        new = {}
        for name in dict_of_dicts:

          # Make plt settings a sub-dict
          d = {'plt':dict_of_dicts[name]}
          # Set default lw
          if 'lw' not in d['plt']: d['plt']['lw'] = 2
          # Extract from  plt-dict 'transf' and 'step' fields
          try:             raise_field_lvl(d,'transf')
          except KeyError: d['transf'] = lambda x: x
          try:             raise_field_lvl(d,'step')
          except KeyError: pass

          try: stat = getattr(stats,name) # Check if stat is there.
          # Fails e.g. if assess(0) before creating stat.
          except AttributeError: continue
          try: val0 = stat[0] # Check if stat[0] has been written
          # Fails e.g. if store_u==False and k_tmp==None (init)
          except KeyError:       continue

          if isinstance(stat,np.ndarray):
            if len(stat) != (chrono.KObs+1): raise TypeError(
                "Only len=(KObs+1) ndarrays supported " +
                "[use FAU_series for len=(K+1)]")
            d['data'] = np.full(chrono.pKObs, nan)
            tt_       = chrono.ttObs[arange(chrono.pKObs)]
          else:
            d['data'] = np.full(chrono.pK, nan)
            tt_       = chrono.tt   [arange(chrono.pK)]
          d['data'][0] = d['transf'](val0)
          d['h'],      = ax.plot(tt_,d['data'],**d['plt'])
          new[name]    = d
        return new

      self.d_RMS    = init_axd(self.ax_RMS, d_RMS)
      self.ax_RMS.set_ylabel('RMS')

      self.d_Uni    = init_axd(self.ax_Uni, d_Uni)
      self.ax_Uni.set_ylabel('mean of marginal\n $\sigma$-normalized values',
          fontsize='small', labelpad=0)
      
      self.ax_Uni.set_xlabel('time (t)')



    #####################
    # Weighted histogram
    #####################
    if 4 in only and E is not None and stats._has_w:
      fig_hw, ax_hw = freshfig(4,(6,3), gridspec_kw={'bottom':.15})
      win_title(fig_hw,"Weight histogram")
      set_figpos('2321')
      ax_hw.set_xscale('log')
      ax_hw.set_xlabel('weigth [× N]')
      ax_hw.set_ylabel('count')
      if len(E)<10001:
        hist   = ax_hw.hist(stats.w[0])[2]
        N      = len(E)
        xticks = 1/N * 10**arange(-4,log10(N)+1)
        xtlbls = array(['$10^{'+ str(int(log10(w*N))) + '}$' for w in xticks])
        xtlbls[xticks==1/N] = '1'
        ax_hw.set_xticks(xticks)
        ax_hw.set_xticklabels(xtlbls)
        self.fig_hw = fig_hw
        self.ax_hw  = ax_hw
        self.hist   = hist
      else:
        not_available_text(ax_hw,'Not computed (N > threshold)')



    #####################
    # User-defined state
    #####################
    if 9 in only and hasattr(setup,'liveplotting'):
      self.fig_custom, self.custom = setup.liveplotting(stats,key,E,P)
      win_title(self.fig_custom,"Custom plot")
      set_figpos('2322')
      plot_pause(0.01)

    self.prev_k = 0
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


  def update(self,key,E=None,P=None,**kwargs):
    """Update liveplots"""
    if self.skip_plotting(): return

    k,kObs,f_a_u = key

    stats = self.stats
    mu    = stats.mu
    m     = self.xx.shape[1]

    ii,wrap = setup_wrapping(m)
    
    #####################
    # Correlation plot
    #####################
    if hasattr(self, 'fig_C') and plt.fignum_exists(self.fig_C.number):
      plt.figure(self.fig_C.number)

      if E is not None:
        C = np.cov(E,rowvar=False)
      else:
        assert P is not None
        C = P.full if isinstance(P,CovMat) else P
        C = C.copy()
      std = sqrt(diag(C))
      C  /= std[:,None]
      C  /= std[None,:]
      if tri_corr:
        C = np.ma.masked_where(self.mask, C)
      self.im_C.set_data(C)

      # Auto-corr function
      ACF = circulant_ACF(C)
      AAF = circulant_ACF(C,do_abs=True)
      self.line_AC.set_ydata(ACF)
      self.line_AA.set_ydata(AAF)

      plot_pause(0.01)


    #####################
    # Spectral error plot
    #####################
    if hasattr(self, 'fig_S') and plt.fignum_exists(self.fig_S.number) and self.do_spectral_error:
      plt.figure(self.fig_S.number)
      msft = abs(stats.umisf[k])
      sprd =     stats.svals[k]
      self.line_sprd.set_ydata(sprd)
      self.line_msft.set_ydata(msft)
      update_ylim(msft, self.ax_S)
      plot_pause(0.01)


    #####################
    # Diagnostics
    #####################
    if hasattr(self,'fig_pulse') and plt.fignum_exists(self.fig_pulse.number):
      plt.figure(self.fig_pulse.number)
      chrono = self.setup.t

      # Indices with shift
      kkU    = arange(chrono.pK) + max(0,k-chrono.pK)
      ttU    = chrono.tt[kkU]
      # Indices for Obs-times
      kkA    = kkU[0] <= chrono.kkObs
      kkA   &=           chrono.kkObs <= kkU[-1]

      def update_axd(ax,dict_of_dicts):
        ax.set_xlim(ttU[0], ttU[0] + 1.1*(ttU[-1]-ttU[0]))

        for name, d in dict_of_dicts.items():
          stat = getattr(stats,name)
          if isinstance(stat,np.ndarray):
            tt_              = chrono.ttObs[kkA]
            d['data']        = stat        [kkA]
            if d.get('step',False):
              # Creat "step"-style graph
              d['data']      = d['data'].repeat(2)
              tt_            = tt_      .repeat(2)
              right          = tt_[-1] # use ttU[-1] for continuous extrapolation
              tt_            = np.hstack([ttU[0], tt_[0:-2], right])
            elif stat.dtype == 'bool':
              # Creat "impulse"-style graph
              tt_            = tt_      .repeat(3)
              d['data']      = d['data'].repeat(3)
              tt_     [2::3] = nan
              d['data'][::3] = False
          else:
            tt_              = ttU
            if stat.store_u:
              d['data']      = stat[kkU]
            else: # store .u manually
              tmp = stat[k]
              if self.prev_k not in [k, k-1]:
                # Reset display
                d['data'][:] = nan
              if k >= chrono.pK:
                # Rolling display
                d['data']    = roll_n_sub(d['data'], tmp, -1)
              else:
                # Initial display: append
                d['data'][k] = tmp
          d['data'] = d['transf'](d['data'])
          d['h'].set_data(tt_,d['data'])

      def rm_absent(ax,dict_of_dicts):
        for name in list(dict_of_dicts):
          d = dict_of_dicts[name]
          if not np.any(np.isfinite(d['data'])):
            d['h'].remove()
            del dict_of_dicts[name]
        if dict_of_dicts:
          ax.legend(loc='upper left')

      update_axd(self.ax_RMS,self.d_RMS)
      update_axd(self.ax_Uni,self.d_Uni)

      #if k%(chrono.pK/5) <= 1:
      update_ylim([d['data'] for d in self.d_RMS.values()], self.ax_RMS,
          bottom=0,      cC=0.2,cE=0.9)
      update_ylim([d['data'] for d in self.d_Uni.values()], self.ax_Uni,
          Max=4, Min=-4, cC=0.3,cE=0.9)

      # Check which diagnostics are present
      if (not self.has_checked_presence) and (k>=chrono.kkObs[0]):
        rm_absent(self.ax_RMS,self.d_RMS)
        rm_absent(self.ax_Uni,self.d_Uni)
        self.has_checked_presence = True

      plot_pause(0.01)



    #####################
    # Weight histogram
    #####################
    if kObs and hasattr(self, 'fig_hw') and plt.fignum_exists(self.fig_hw.number):
      plt.figure(self.fig_hw.number)
      ax_hw     = self.ax_hw
      _         = [b.remove() for b in self.hist]
      w         = stats.w[k]
      N         = len(w)
      wmax      = w.max()
      bins      = exp(linspace(log(1e-5/N), log(1), int(N/20)))
      counted   = w>bins[0]
      nC        = np.sum(counted)
      nn,_,pp   = ax_hw.hist(w[counted], bins=bins, color='b')
      self.hist = pp
      #thresh   = '#(w<$10^{'+ str(int(log10(bins[0]*N))) + '}/N$ )'
      ax_hw.set_title('N: {:d}.   N_eff: {:.4g}.   Not shown: {:d}. '.\
          format(N, 1/(w@w), N-nC))
      update_ylim([nn], ax_hw, cC=True)
      plot_pause(0.01)



    #####################
    # User-defined state
    #####################
    if hasattr(self,'fig_custom') and plt.fignum_exists(self.fig_custom.number):
      plt.figure(self.fig_custom.number)
      self.custom(key,E,P)
      plot_pause(0.01)

    # Trackers
    self.prev_k = k




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


def setup_wrapping(m,periodic=True):
  """
  Make periodic indices and a corresponding function
  (that works for ensemble input).
  """

  if periodic:
    ii = np.hstack([-0.5, arange(m), m-0.5])
    def wrap(E):
      midpoint = (E[[0],...] + E[[-1],...])/2
      return ccat(midpoint,E,midpoint)

  else:
    ii = arange(m)
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

# TODO: rename "xtrma"
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



def update_ylim(data,ax,bottom=None,top=None,Min=-1e20,Max=+1e20,cC=0,cE=1):
  """
  Update ylim's intelligently, mainly by computing
  the low/high percentiles of the data.
  - data: iterable of arrays for computing percentiles.
  - bottom/top: override values.
  - Max/Min: bounds.
  - cE: exansion (widenting) rate ∈ [0,1].
      Default: 1, which immediately expands to percentile.
  - cC: compression (narrowing) rate ∈ [0,1].
      Default: 0, which does not allow compression.
  Despite being a little involved,
  the cost of this subroutine is typically not substantial.
  """

  #
  def worth_updating(a,b,curr):
    # Note: should depend on cC and cE
    d = abs(curr[1]-curr[0])
    lower = abs(a-curr[0]) > 0.002*d
    upper = abs(b-curr[1]) > 0.002*d
    return lower and upper
  #
  current = ax.get_ylim()
  # Find "reasonable" limits (by percentiles), looping over data
  maxv = minv = -np.inf # init
  for d in data:
    d = d[np.isfinite(d)]
    if len(d):
      minv, maxv = np.maximum([minv, maxv], \
          array([-1, 1]) * np.percentile(d,[1,99]))
  minv *= -1
  minv, maxv = stretch(minv,maxv,1.02)
  # Pry apart equal values
  if np.isclose(minv,maxv):
    maxv += 0.5
    minv -= 0.5
  # Set rate factor as compress or expand factor. 
  c0 = cC if minv>current[0] else cE
  c1 = cC if maxv<current[1] else cE
  # Adjust
  minv = np.interp(c0, (0,1), (current[0], minv))
  maxv = np.interp(c1, (0,1), (current[1], maxv))
  # Bounds
  maxv = min(Max,maxv)
  minv = max(Min,minv)
  # Overrides
  if top    is not None: maxv = top
  if bottom is not None: minv = bottom
  # Set (if anything's changed)
  #if worth_updating(minv,maxv,current):
    #ax.set_ylim(minv,maxv)
  ax.set_ylim(minv,maxv)


def set_ilim(ax,i,data,zoom=1.0):
  """Set bounds (taken from data) on axis i.""" 
  Min  = data[:,i].min()
  Max  = data[:,i].max()
  lims = round2sigfig([Min, Max])
  lims = inflate_ens(lims,1/zoom)
  if i is 0: ax.set_xlim(lims)
  if i is 1: ax.set_ylim(lims)
  if i is 2: ax.set_zlim(lims)

def set_ilabel(ax,i):
  if i is 0: ax.set_xlabel('x')
  if i is 1: ax.set_ylabel('y')
  if i is 2: ax.set_zlabel('z')



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
  chrono = stats.setup.t

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
  chrono = stats.setup.t
  xx     = stats.xx
  m      = xx.shape[1]
  dims   = equi_spaced_integers(m, min(m, 10))
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
  ax_K.set_xlabel('time (t)')
  ax_K.set_ylabel('mean of marginal\n $\sigma$-normalized values',
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

  m = xx.shape[1]

  if chrono!=None:
    kk,_ = get_plot_inds(xx,chrono,mult=40,**kwargs)
    tt   = chrono.tt[kk]
    ax.set_ylabel('Time (t)')
  else:
    pK   = estimate_good_plot_length(xx,mult=40)
    tt   = arange(pK)
    ax.set_ylabel('Time indices (k)')

  plt.contourf(arange(m),tt,xx[kk],25)
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

  chrono = stats.setup.t
  m      = stats.xx.shape[1]

  err   = mean( abs(stats.err  .a) ,0)
  sprd  = mean(     stats.mad  .a  ,0)
  umsft = mean( abs(stats.umisf.a) ,0)
  usprd = mean(     stats.svals.a  ,0)

  ax_r = plt.subplot(311)
  ax_r.plot(          arange(m),               err,'k',lw=2, label='Error')
  if m<10**3:
    ax_r.fill_between(arange(m),[0]*len(sprd),sprd,alpha=0.7,label='Spread')
  else:
    ax_r.plot(        arange(m),              sprd,alpha=0.7,label='Spread')
  #ax_r.set_yscale('log')
  ax_r.set_title('Element-wise error comparison')
  ax_r.set_xlabel('Dimension index (i)')
  ax_r.set_ylabel('Time-average (_a) magnitude')
  ax_r.set_ylim(bottom=mean(sprd)/10)
  ax_r.set_xlim(right=m-1); add_endpoint_xtick(ax_r)
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
    ax_s.set_xlim(right=m-1); add_endpoint_xtick(ax_s)
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
  chrono = stats.setup.t

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
    m     = ranks.shape[1]
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
      w  = array([ w[arange(K),ranks[arange(K),i]] for i in range(m)])
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
    n      = fig.number
    string += " [" + str(n) + "]"
  fig.canvas.set_window_title(string)

def set_figpos(loc):
  """
  Place figure on screen, where 'loc' can be either
    NW, E, ...
  or
    4 digits (as str or int) to define grid m,n,i,j.
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
  m,n,i,j = [int(x) for x in loc[:4]]
  assert m>=i>0 and n>=j>0
  h0   -= (m-1)*25
  yoff  = 25*(i-1)
  if i>1:
    yoff += 25
  place((j-1)*w0/n, yoff + (i-1)*h0/m, w0/n, h0/m)


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
  for p in ['label','color','visible']:
    lines[p] = [plt.getp(x,p) for x in lines['handle']]
  # Put into pandas for some reason
  lines = pd.DataFrame(lines)
  # Rm those that start with _
  lines = lines[~lines.label.str.startswith('_')]

  # Adjust labels
  if numbering: lines['label'] = [str(i)+': '+lb for i,lb in enumerate(lines['label'])]
  if txtwidth:  lines['label'] = [textwrap.fill(n,width=txtwidth) for n in lines['label']]

  # Set state. BUGGY? sometimes causes MPL complaints after clicking boxes
  if state is not None:
    state = array(state).astype(bool)
    lines.visible = state
    for i,x in enumerate(state):
      lines['handle'][i].set_visible(x)

  # Setup buttons
  # When there's many, the box-sizing is awful, but difficult to fix.
  W = 0.23 * txtwidth/15 * txtsize/10
  N = len(lines)
  H = min(1,0.07*N)
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
def toggle_viz(h,prompt=True,legend=True):
  """
  Toggle visibility of the handle h,
  which can also be a list of handles.
  """

  # Core functionality: turn on/off
  is_viz = not h.get_visible()
  h.set_visible(is_viz)

  if prompt:
    input("Press <Enter> to continue...")

  if legend:
    if is_viz:
      try:
        h.set_label(h.old_label)
      except AttributeError:
        pass
    else:
      h.old_label = h.get_label()
      h.set_label('_nolegend_')

    # Legend update
    ax = h.axes
    with warnings.catch_warnings():
      warnings.simplefilter("error",category=UserWarning)
      try:
        ax.legend()
      except UserWarning:
        # If there's no remaining legend entries, ax.legend() throws warning.
        # And yet it does NOT remove the last entry!
        # Also ax._legend.remove() won't work either coz now _legend does not exist!
        lh = ax.legend('NOT ASSOCIATED TO ANY OBJECT') # Spur mpl back into action
        lh.remove()

  plt.pause(0.04)
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
    



