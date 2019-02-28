from common import *



class LivePlot:
  """
  Live plotting manager. Deals with
   - Pause, skip
   - Which liveploters to call
   - Figure window (title and number)
  Everything else is managed by each specific liveplotter.
  """
  def __init__(self,stats,fignums,key0=(0,None,'u'),E=None,P=None):
    """
    Initialize plots.

    - fignums: List of figures requested for this config,
      referred to by their number. Can also specify using:
         - "default": show all recommended plots (specified here and in HMM)
         - "all"    : show all possible plots
         - boolean  : show all default, or none.
    """

    # HMM-independent
    potential_figures = OrderedDict(
        # name               num  show_by_default init(ializer)
        sliding_diagnostics = (1, 1,              LP_sliding_diagnostics),
        weight_histogram    = (4, 1,              LP_weight_histogram   ),
        )
    # HMM-specific
    for name, (num, dflt, init) in getattr(stats.HMM,'liveplotters',dict()).items():
      assert num>10, "Liveplotters specified in the HMM should have fignum>10."
      potential_figures[name] = (num, dflt, init)

    # Figures requested for this config.
    if isinstance(fignums,str):
      fn = fignums.lower()                               # Yields:
      if   "all" == fn:              fignums = range(99) # All potential_figures
      elif "default" in fn:          fignums = []        # All show_by_default
    elif hasattr(fignums,'__len__'): fignums = fignums   # This list only
    elif fignums:                    fignums = []        # All show_by_default
    else:                            fignums = [None]    # None

    # Call figure inits
    self.any_figs = False
    self.figures = OrderedDict()
    for name, (num, show_by_default, init) in potential_figures.items():
      if num in fignums or (fignums==[] and show_by_default):

        # Startup message
        if not self.any_figs:
          print('Initializing liveplotting...')
          print('Hit <Space> to pause/step.')
          print('Hit <Enter> to resume/skip.')
          self.paused = False
          self.skipping = False
          self.any_figs = True

        updater = init(num,stats,key0,E,P)
        if plt.fignum_exists(num):
          self.figures[name] = (num, updater)
          plt.figure(num).canvas.set_window_title("%s [%d]"%(name,num))
          plot_pause(0.01)


  def update(self,key,E,P):
    """Update liveplots"""

    # Check if there are still open figures
    if self.any_figs:
      open_figns = plt.get_fignums()
      live_figns = set(num for (num, updater) in self.figures.values())
      self.any_figs = bool(live_figns.intersection(open_figns))
    # If no open figures: don't update
    if not self.any_figs:
      return

    # Playback control
    SPACE  = b' '   
    ENTERs = [b'\n', b'\r'] # Linux + Windows
    def pause():
      "Loop until user decision is made."
      ch = read1() 
      while True:
        if ch in ENTERs:
          self.paused = False
        if ch in ENTERs + [SPACE]:
          break
        ch = read1()
        # Pause to enable zoom, pan, etc. of mpl GUI
        plot_pause(0.01) # Don't use time.sleep()!
    #
    if self.paused:
      pause()
    else:
      # Set switches for pause & skipping
      ch = read1()
      if ch==SPACE: # Turn ON pause & turn OFF skipping.
        self.paused = True
        self.skipping = False
        pause()
      elif ch in ENTERs: # Toggle skipping
        self.skipping = not self.skipping 

        
    if not self.skipping:
      # Update figures
      for name, (num, updater) in self.figures.items():
        if plt.fignum_exists(num):
          plt.figure(num)
          updater(key,E,P)
          plot_pause(0.01)



# TODO:
# - iEnKS diagnostics don't work at all when store_u=False
star = "${}^*$"
class LP_sliding_diagnostics:

  def __init__(self,fignum,stats,key0,E,P,Tplot=None):
      GS = {'left':0.125,'right':0.76}
      fig, (ax1, ax2) = freshfig(fignum, (5,3.5), nrows=2, sharex=True, gridspec_kw=GS)

      ax1.set_ylabel('RMS')
      ax2.set_ylabel('Values') 
      ax2.set_xlabel('Time (t)')

      self.T_lag, K_lag, a_lag, self.dt_margin = validate_lag(Tplot, stats.HMM.t)

      def init_ax(ax,style_table):
        plotted_lines = OrderedDict()
        for name in style_table:

            # SKIP -- if stats[name] is not in existence
            # Note: The nan check/deletion comes after the first kObs.
            try: stat = getattr(stats,name)
            except AttributeError: continue
            # try: val0 = stat[key0[0]]
            # except KeyError: continue
            # PS: recall (from series.py) that even if store_u is false, stat[k] is
            # still present if liveplotting=True via the k_tmp functionality.
            
            # Unpack style
            ln = {}
            ln['transf'] = style_table[name][0]
            ln['shape']  = style_table[name][1]
            ln['plt']    = style_table[name][2]

            # Create series
            K_plot, ln['plot_u'] = determine_K_plot(stat,key0,K_lag,a_lag)
            ln['data'] = RollingArray(K_plot)
            ln['tt']   = RollingArray(K_plot)

            # Plot (init)
            ln['handle'], = ax.plot(ln['tt'],ln['data'],**ln['plt'])

            # Plotting only nans yield ugly limits. Revert to defaults.
            ax.set_xlim(0,1)
            ax.set_ylim(0,1)

            plotted_lines[name] = ln
        return plotted_lines

      # Plot
      self.d1 = init_ax(ax1, stats.style1);
      self.d2 = init_ax(ax2, stats.style2);

      # Horizontal line at y=0
      self.baseline0, = ax2.plot(ax2.get_xlim(),[0,0],c=0.5*ones(3),lw=0.7,label='_nolegend_')

      # Store
      self.ax1 = ax1
      self.ax2 = ax2
      self.stats = stats
      self.not_yet_seen_analysis = True

      # Finalize init
      self(key0,E,P)


  # Update plot
  def __call__(self,key,E,P):
      k, kObs, f_a_u = key

      stats  = self.stats
      chrono = stats.HMM.t
      ax1    = self.ax1
      ax2    = self.ax2

      def update_arrays(plotted_lines):
        for name, ln in plotted_lines.items():
          stat = getattr(stats,name)
          t    = chrono.tt[k] # == chrono.ttObs[kObs]
          if isinstance(stat,FAU_series):
            #                 ln['data'] will contain duplicates for f/a times.
            if ln['plot_u']:
              val = stat[key]
              ln['tt']  .insert(k   , t)
              ln['data'].insert(k   , ln['transf'](val))
            elif 'u' not in f_a_u:
              val = stat[key]
              ln['tt']  .insert(kObs, t)
              ln['data'].insert(kObs, ln['transf'](val))
          else:
            if 'a' in f_a_u: # ln['data'] will not contain duplicates, coz only 'a' is input.
              val = stat[kObs]
              ln['tt']  .insert(kObs, t)
              ln['data'].insert(kObs, ln['transf'](val))
            elif 'f' in f_a_u:
              pass

      def update_plot_data(ax,plotted_lines):

          def bend_into(shape, xx, yy):
            # Get arrays. Repeat (to use for intermediate nodes). 
            yy = yy.array.repeat(3)
            xx = xx.array.repeat(3)
            if len(xx)==0:
              pass # shortcircuit any modifications
            elif shape == 'step':
              yy = np.hstack([yy[1:], nan]) # roll leftward
            elif shape == 'dirac':
              nonlocal nDirac
              axW      = np.diff(ax.get_xlim())
              yy[0::3] = False          # set datapoin to 0
              xx[2::3] = nan            # make datapoint disappear
              xx      += nDirac*axW/100 # offset datapoint horizontally
              nDirac  +=1
            return xx, yy

          nDirac = 1
          for name, ln in plotted_lines.items():
            ln['handle'].set_data(*bend_into(ln['shape'], ln['tt'], ln['data']))

      def finalize_init(ax,plotted_lines,mm):
          # Rm lines that only contain NaNs
          for name in list(plotted_lines):
            ln = plotted_lines[name]
            if not np.any(np.isfinite(ln['data'])):
              ln['handle'].remove()
              del plotted_lines[name]
          # Add legends
          if plotted_lines:
            ax.legend(loc='upper left', bbox_to_anchor=(1.01, 1),borderaxespad=0)
            if mm:
              ax.annotate(star+": mean of\nmarginals",
                  xy=(0,-1.5/len(plotted_lines)),
                  xycoords=ax.get_legend().get_frame(),
                  bbox=dict(alpha=0.0), fontsize='small')
          plot_pause(0.01) # coz placement of annotate needs flush sometimes

      # Insert current stats
      update_arrays(self.d1)
      update_arrays(self.d2)

      # Plot
      update_plot_data(ax1, self.d1)
      update_plot_data(ax2, self.d2)

      # Set x-limits (time)
      sliding_xlim(ax1, self.d1['rmse']['tt'], self.dt_margin, self.T_lag)
      self.baseline0.set_xdata(ax1.get_xlim())

      # Set y-limits
      data1 = [ln['data'].array for ln in self.d1.values()]
      data2 = [ln['data'].array for ln in self.d2.values()]
      ax1.set_ylim(0, d_ylim(data1, ax1,                cC=0.2,cE=0.9)[1])
      ax2.set_ylim(  *d_ylim(data2, ax2, Max=4, Min=-4, cC=0.3,cE=0.9))

      # Init legend. Rm nan lines.
      if self.not_yet_seen_analysis and 'a'==f_a_u:
         self.not_yet_seen_analysis = False
         finalize_init(ax1, self.d1, False)
         finalize_init(ax2, self.d2, True)



def sliding_xlim(ax, tt, dt_margin, lag):
  if tt.nFilled==0: return        # Quit
  t1, t2 = tt.span()              # Get suggested span.
  s1, s2 = ax.get_xlim()          # Get previous lims.
  if t1==t2:                      # If zero span (eg tt holds single 'f' and 'a'):
    t1 -= 1                       #   add width
    t2 += 1                       #   add width
  elif np.isnan(t1):              # If user has skipped (too much):
    s2    -= dt_margin            #   Correct for dt_margin.
    span   = s2-s1                #   Compute previous span
    if span < lag:                #   If span<lag:
      span  += (t2-s2)            #     Grow by "dt".
    span   = min(lag, span)       #   Bound
    t1     = t2 - span            #   Set span.
  ax.set_xlim(t1, t2 + dt_margin) # Set xlim to span


class LP_weight_histogram:

  def __init__(self,fignum,stats,key0,E,P):
    if hasattr(stats,'w'):
      fig, ax = freshfig(fignum, (6,3), gridspec_kw={'bottom':.15})
      ax.set_xscale('log')
      ax.set_xlabel('Weigth × N')
      ax.set_ylabel('Count')

      w0 = stats.w[key0]
      if len(w0)<10001:
        hist   = ax.hist(w0)[2]
        N      = len(w0)
        xticks = 1/N * 10**arange(-4,log10(N)+1)
        xtlbls = array(['$10^{'+ str(int(log10(w*N))) + '}$' for w in xticks])
        xtlbls[xticks==1/N] = '1'
        ax.set_xticks(xticks)
        ax.set_xticklabels(xtlbls)
        self.ax    = ax
        self.hist  = hist
        self.stats = stats
      else:
        not_available_text(ax,'Not computed (N > threshold)')

  # Update plot
  def __call__(self,key,E,P):
    if 'a' in key[2]:
      w         = self.stats.w[key]
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


class LP_spectral_errors:

  def __init__(self,fignum,stats,key0,E,P):
    fig, ax = freshfig(fignum, (6,3))
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
      if np.any(np.isinf(self.msft[key0])):
        not_available_text(ax, "Spectral stats not finite")
        do_spectral_error = False
      else:
        do_spectral_error = True

    if do_spectral_error:
      M = len(self.msft[key0])
      self.line_msft, = ax.plot(arange(M),ones(M),'k',lw=2,label='Error')
      self.line_sprd, = ax.plot(arange(M),ones(M),'b',lw=2,label='Spread',alpha=0.9)
      ax.get_xaxis().set_major_locator(MaxNLocator(integer=True))
      ax.legend()

    self.do_spectral_error = do_spectral_error
    self.ax  = ax
    self(key0,E,P)

  # Update plot
  def __call__(self,key,E,P):
    if self.do_spectral_error:
      msft = abs(self.msft[key])
      sprd =     self.sprd[key]
      self.line_sprd.set_ydata(sprd)
      self.line_msft.set_ydata(msft)
      self.ax.set_ylim(*d_ylim(msft))


class LP_correlations:

  # Whether to show half/full (symmetric) corr matrix.
  half = True

  def __init__(self,fignum,stats,key0,E,P):
    GS = {'height_ratios':[4, 1],'hspace':0.09,'top':0.95}
    fig, (ax,ax2) = freshfig(fignum, (5,6), nrows=2, gridspec_kw=GS)

    Nx = len(stats.mu[key0])
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

      self.ax      = ax
      self.ax2     = ax2
      self.im      = im
      self.line_AC = line_AC
      self.line_AA = line_AA
      self.mask    = mask

      self(key0,E,P)
    else:
      not_available_text(ax)

  # Update plot
  def __call__(self,key,E,P):
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


def plot_pause(interval):
  """Similar to plt.pause()"""

  # plt.pause(0) just seems to freeze execution.
  if interval==0:
    return

  try:
    # Implement plt.pause() that doesn't focus window, c.f.
    # github.com/matplotlib/matplotlib/issues/11131, so.com/q/45729092.
    # Only necessary for some platforms (e.g. Windows) and mpl versions.
    # Even then, mere figure creation may steal the focus if script was
    # launched with `$ python example_1.py` rather than ipython's `run`.
    from matplotlib import _pylab_helpers
    def _plot_pause(interval,  focus_figure=True):
        canvas = plt.gcf().canvas
        manager = canvas.manager
        if manager is not None:
            if canvas.figure.stale:
                canvas.draw_idle()
            if focus_figure:
                plt.show(block=False)
            if not is_notebook:
                canvas.start_event_loop(interval)
        else:
            time.sleep(interval)
    _plot_pause(interval, focus_figure=False)

  except:
    # Jupyter notebook support (SO.com/q/34486642)
    # Note: no longer needed with the above _plot_pause()?
    plt.gcf().canvas.draw()
    time.sleep(0.1)



def sliding_marginals(
    obs_inds     = [],
    dims         = [],
    labels       = [],
    Tplot        = None,
    ens_props    = dict(alpha=0.4),
    zoomy        = 1.0,
    pause_f      = 0.0,
    pause_a      = 0.0,
    pause_u      = 0.0,
    ):

  def init(fignum,stats,key0,E,P):
    xx, yy, mu, var, chrono = stats.xx, stats.yy, stats.mu, stats.var, stats.HMM.t

    # Lag settings:
    T_lag, K_lag, a_lag, dt_margin = validate_lag(Tplot, chrono)
    K_plot, plot_u = determine_K_plot(stats.mu,key0,K_lag,a_lag)
    # Extend K_plot forther for adding blanks in resampling (PartFilt):
    has_w = hasattr(stats,'w')
    if has_w: K_plot += a_lag

    # Pre-process dimensions
    DimsX = arange(xx.shape[-1]) if dims==[] else dims              # Chose marginal dims to plot
    iiY   = [i for i,m in enumerate(obs_inds) if m in DimsX]        # Rm inds of obs if not in DimsX
    DimsY = [m for i,m in enumerate(obs_inds) if m in DimsX]        # Rm obs_inds    if not in DimsX
    DimsY = [DimsY.index(m) if m in DimsY else None for m in DimsX] # Get dim (within y) of each x

    Nx = len(DimsX)
    Ny = len(iiY)

    # Set up figure, axes
    fig, axs = freshfig(fignum, (5,7), nrows=Nx, sharex=True)
    if Nx==1: axs = [axs]
    set_figpos('2322')

    # Tune plots
    axs[0].set_title("Marginal time series.")
    for ix, (m,ax) in enumerate(zip(DimsX,axs)):
      ax.set_ylim(*stretch(*span(xx[:,m]), 1/zoomy))
      if labels==[]: ax.set_ylabel("Dim %d"%m)
      else:          ax.set_ylabel(labels[ix])
    axs[-1].set_xlabel('Time (t)')

    # Allocate
    d = Bunch() # data arrays
    h = Bunch() # plot handles
    # Why "if True" ? Just to indent the rest of the line...
    if True         : d.t  = RollingArray((K_plot,))          ; 
    if True         : d.x  = RollingArray((K_plot,Nx))        ; h.x  = []
    if True         : d.y  = RollingArray((K_plot,Ny))        ; h.y  = []
    if E is not None: d.E  = RollingArray((K_plot,len(E),Nx)) ; h.E  = []
    if P is not None: d.mu = RollingArray((K_plot,Nx))        ; h.mu = []
    if P is not None: d.s  = RollingArray((K_plot,2,Nx))      ; h.s  = []

    # Plot (invisible coz everything here is nan, for the moment).
    for ix, (m, iy, ax) in enumerate(zip(DimsX,DimsY,axs)):
      if True     : h.x  +=  ax.plot(d.t, d.x [:,ix]  , 'k')
      if iy!=None : h.y  +=  ax.plot(d.t, d.y [:,iy]  , 'g*', ms=10)
      if 'E'  in d: h.E  += [ax.plot(d.t, d.E [:,:,ix], **ens_props)]
      if 'mu' in d: h.mu +=  ax.plot(d.t, d.mu[:,ix]  , 'b')
      if 's'  in d: h.s  += [ax.plot(d.t, d.s [:,:,ix], 'b--',lw=1)]


    def update(key,E,P):
      nonlocal d 
      k,kObs,f_a_u = key

      EE = duplicate_with_blanks_for_resampled(E, DimsX, key, has_w)

      # Roll data array
      for Ens in EE: # If E is duplicated, so must the others be.
        if 'E'  in d: d.E .insert(k, Ens)
        if 'mu' in d: d.mu.insert(k, mu[key][DimsX])
        if 's'  in d: d.s .insert(k, mu[key][DimsX] + [[1],[-1]]*sqrt(var[key][DimsX]))
        if True     : d.t .insert(k, chrono.tt[k])
        if True     : d.y .insert(k, yy[kObs,iiY] if kObs is not None else nan*ones(Ny))
        if True     : d.x .insert(k, xx[k,DimsX])

      # Update graphs
      for ix, (m, iy, ax) in enumerate(zip(DimsX,DimsY,axs)):
        sliding_xlim(ax, d.t, dt_margin, T_lag)
        if True:       h.x [ix]   .set_data(d.t, d.x [:,ix])
        if iy!=None:   h.y [iy]   .set_data(d.t, d.y [:,iy])
        if 'mu' in d:  h.mu[ix]   .set_data(d.t, d.mu[:,ix])
        if 's'  in d: [h.s [ix][b].set_data(d.t, d.s [:,b,ix]) for b in [0,1]]
        if 'E'  in d: [h.E [ix][n].set_data(d.t, d.E [:,n,ix]) for n in range(len(E))]
        if 'E'  in d: update_alpha(key, stats, h.E[ix])

      if   'f' in f_a_u: plot_pause(pause_f)
      elif 'a' in f_a_u: plot_pause(pause_a)
      elif 'u' in f_a_u: plot_pause(pause_u)

      return # end update()

    # Finalize init
    update(key0,E,P)

    return update # end init()
  return init # end sliding_marginals()


def phase3d(
    obs_inds     = [],
    dims         = [],
    labels       = [],
    Tplot        = None,
    ens_props    = dict(alpha=0.4),
    zoom         = 1.5,
    pause_f      = 0.0,
    pause_a      = 0.0,
    pause_u      = 0.0,
    ):

  M = 3 # Only applicable for 3d plots

  def init(fignum,stats,key0,E,P):

    # Extract data arrays
    k, kObs, f_a_u = key0
    xx, yy, mu, var, chrono = stats.xx, stats.yy, stats.mu, stats.var, stats.HMM.t

    # Lag settings:
    T_lag, K_lag, a_lag, dt_margin = validate_lag(Tplot, chrono)
    K_plot, plot_u = determine_K_plot(stats.mu,key0,K_lag,a_lag)
    # Extend K_plot forther for adding blanks in resampling (PartFilt):
    has_w = hasattr(stats,'w')
    if has_w: K_plot += a_lag
    
    # Dimension settings
    nonlocal dims, labels
    if dims   == []: dims   = [0,1,2]
    if labels == []: labels = ["$x_%s$"%i for i in "123"]
    assert len(dims)==M

    # Set up figure, axes
    fig = plt.figure(fignum, figsize=(5,5))
    ax3 = plt.subplot(111, projection='3d')
    ax3.set_facecolor('w')
    ax3.set_title("Phase space trajectories" +
        ("" if plot_u else "\n(obs times only)"))
    # Tune plot
    for s,i in zip(labels, dims):
      set_ilim(ax3, i, *stretch(*span(xx[:,i]),1/zoom))
      ax3.set_xlabel(labels[0])
      ax3.set_ylabel(labels[1])
      ax3.set_zlabel(labels[2])

    # Allocate
    d = Bunch() # data arrays
    h = Bunch() # plot handles
    s = Bunch() # scatter handles
    if E is not None             : d.E  = RollingArray((K_plot,len(E),M)); h.E = []
    if P is not None             : d.mu = RollingArray((K_plot,M))
    if True                      : d.x  = RollingArray((K_plot,M))
    if list(obs_inds)==list(dims): d.y  = RollingArray((K_plot,M))

    # Plot tails (invisible coz everything here is nan, for the moment).
    if 'E'  in d: h.E  += [ax3.plot(*xn    , **ens_props)[0] for xn in np.transpose(d.E,[1,2,0])]
    if 'mu' in d: h.mu  =  ax3.plot(*d.mu.T, 'b' ,lw=2)  [0]
    if True     : h.x   =  ax3.plot(*d.x .T, 'k' ,lw=3)  [0]
    if 'y'  in d: h.y   =  ax3.plot(*d.y .T, 'g*',ms=14) [0]

    # Scatter
    if 'E'  in d: s.E   =  ax3.scatter(*E.T  [dims],s=3 **2, c=[ hn.get_color() for hn in h.E])
    if 'mu' in d: s.mu  =  ax3.scatter(*nan*ones(M),s=8 **2, c=h.mu.get_color()               )
    if True     : s.x   =  ax3.scatter(*xx[k, dims],s=14**2, c=h.x .get_color(), marker=(5, 1), zorder=99)


    def update(key,E,P):
      k,kObs,f_a_u = key
      show_y = 'y' in d and kObs is not None

      if plot_u:             ind = k
      elif 'u' not in f_a_u: ind = kObs
      else: return

      def update_tail(handle,newdata):
        handle.set_data(newdata[:,0],newdata[:,1])
        handle.set_3d_properties(newdata[:,2])

      EE = duplicate_with_blanks_for_resampled(E, dims, key, has_w)

      for Ens in EE:
        # Roll data array
        if 'E'  in d: d.E .insert(ind, Ens)
        if True     : d.x .insert(ind, xx[k,   dims])
        if 'y'  in d: d.y .insert(ind, yy[kObs,  : ] if show_y else nan*ones(M))
        if 'mu' in d: d.mu.insert(ind, mu[key][dims])

      # Update graph
      s.x._offsets3d = juggle_axes(*d.x[[-1]].T,'z')
      update_tail(h.x, d.x)
      if 'y' in d:
          update_tail(h.y, d.y)
      if 'mu' in d:
          s.mu._offsets3d = juggle_axes(*d.mu[[-1]].T,'z')
          update_tail(h.mu, d.mu)
      else:
          s.E._offsets3d = juggle_axes(*d.E[-1].T,'z')
          for n in range(len(E)):
            update_tail(h.E[n],d.E[:,n,:])
          update_alpha(key, stats, h.E, s.E)

      if   'f' in f_a_u: plot_pause(pause_f)
      elif 'a' in f_a_u: plot_pause(pause_a)
      elif 'u' in f_a_u: plot_pause(pause_u)

      return # end update()

    # Finalize init
    update(key0,E,P)

    return update # end init()
  return init # end phase3d()



def validate_t1t2(chrono,t1,t2): 
  t = chrono
  if t1 is None: t1 = t.tt[0]
  if t2 is None: t2 = t1 + t.Tplot
  assert t1<t2
  t1, t2 = array([t1,t2]).clip(t.tt[0], t.tt[-1])
  return t1, t2

# TODO: add plot_pause for animation
def replay(updater,chrono,t1,t2):
  "Used for non-live plotting (plot_time_series)"
  for k,kObs,t,dt in chrono.ticker:
    if t1 <= t <= t2:
      if kObs is not None:
        updater((k,kObs,'f'),None,None)
        updater((k,kObs,'a'),None,None)
      updater((k,kObs,'u'),None,None)

def plot_time_series(stats,fignum=22,t1=None,t2=None):
  """Plot time series of various statistics."""
  chrono  = stats.HMM.t 
  t1, t2  = validate_t1t2(chrono,t1,t2)
  updater = LP_sliding_diagnostics(fignum,stats,(0,None,'u'),None,None,Tplot=(t2-t1))
  updater.dt_margin = 0
  replay(updater,chrono,t1,t2)

  set_figpos('1313 mac')

def plot_3D_trajectory(stats,fignum=23,t1=None,t2=None,**kwargs):
  """Plot 3D phase-space trajectory."""
  chrono  = stats.HMM.t
  t1, t2  = validate_t1t2(chrono,t1,t2)
  init    = phase3d(Tplot=t2-t1,**kwargs)
  P0      = np.full_like(stats.HMM.X0.C.full, nan) # P0 not None => plots mu
  updater = init(fignum,stats,(0,None,'u'),None,P0)
  replay(updater,chrono,t1,t2)


def update_alpha(key, stats, lines, scatters=None):
  "Adjust color alpha (for particle filters)"

  k,kObs,f_a_u = key
  if kObs is None: return
  if f_a_u=='f': return
  if not hasattr(stats,'w'): return

  # Compute alpha values
  alpha = stats.w[key]
  alpha = (alpha/alpha.max()).clip(0.1,0.4)

  # Set line alpha
  for line, a in zip(lines, alpha):
    line.set_alpha(a)

  # Scatter plot does not have alpha. => Fake it.
  if scatters is not None:
    colors = scatters.get_facecolor()[:,:3]
    if len(colors)==1:
      colors = colors.repeat(len(w),axis=0)
    scatters.set_color(np.hstack([colors, alpha[:,None]]))



def duplicate_with_blanks_for_resampled(E,dims,key,has_w): 
  "Particle filter: insert breaks for resampled particles."
  EE = []

  if has_w and E is not None:
    E = E[:,dims]
    k,kObs,f_a_u = key
    if   f_a_u=='f': pass
    elif f_a_u=='a': _Ea[0] = E[:,0] # Store (1st dim of) ens.
    elif f_a_u=='u' and kObs is not None:
      # Find resampled particles. Insert duplicate ensemble. Write nans (breaks).
      resampled = _Ea[0] != E[:,0]  # Mark as resampled if ens changed.
      EE.append( E.copy() )         # Insert current ensemble (copy to avoid overwriting).
      EE[0][resampled] = nan        # Write breaks

  # Always: append current ensemble
  EE.append(E)                   
  return EE
_Ea = [None] # persistent ens storage


def validate_lag(Tplot, chrono):

  # Defaults
  if Tplot is None:
    Tplot = chrono.Tplot
  
  # Rename
  T_lag = Tplot

  # Validate
  t2 = chrono.tt[-1]
  t1 = max(chrono.tt[0], t2-T_lag)
  T_lag = t2-t1
  
  K_lag = int(T_lag / chrono.dt) + 1 # Lag in indices
  a_lag = K_lag//chrono.dkObs + 1    # Lag in obs indices
  dt_margin = T_lag/20               # Size of leading, empty space

  return T_lag, K_lag, a_lag, dt_margin


def determine_K_plot(stat,key0,K_lag,a_lag):
  """Determine K_plot: the time (in inds) window of plotting,
  i.e. the length of the RollingArray to be used."""

  if isinstance(stat,FAU_series):
    plot_u  = stat.store_u or stat.k_tmp==key0[0]
    K_plot  = 2*a_lag          # f+a series
    if plot_u: K_plot += K_lag # u series.

  else:
    plot_u = False
    K_plot = a_lag

  return K_plot, plot_u



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

  # Some mpl versions don't handle inf limits.
  if not np.isfinite(minv): minv = None
  if not np.isfinite(maxv): maxv = None

  return minv, maxv



from tools.viz import setup_wrapping
def spatial1d(
    obs_inds     = None,
    ens_props    = {'color': 0.7*RGBs['w'],'alpha':0.5},
    periodic     = True,
    conf_mult    = 0,
    pause_f      = 0.5,
    pause_a      = 1.0,
    pause_u      = 0.0,
    ):

  def init(fignum,stats,key0,E,P):

    # Extract data arrays
    xx, yy, mu, Nx = stats.xx, stats.yy, stats.mu, stats.HMM.Nx

    # Make periodic wrapper
    ii, wrap = setup_wrapping(Nx,periodic)

    # Set up figure, axes
    fig, ax = freshfig(fignum, (8,5))
    fig.suptitle("1d amplitude plot")

    # Plot x,y
    nan1 = wrap(nan*ones(Nx))
    line_x, = ax.plot(ii,nan1,'k-',lw=3,label='Truth')
    if obs_inds is not None:
      line_y, = ax.plot(obs_inds, nan*obs_inds,'g*',ms=5,label='Obs')

    # Plot mu/std.dev
    if conf_mult:
      line_mu, = ax.plot(ii,nan1,'b-', lw=2,label='DA mean')
      patch_s  = ax.fill_between(ii,nan1,nan1,
          color='b',alpha=0.4,label=(str(conf_mult) + '$\sigma$'))

    # Plot ensemble/conf
    else:
      if E is not None:
        lines_E  = ax.plot(ii,wrap(E[0] .T),lw=1,**ens_props,label='Ensemble')
        lines_E += ax.plot(ii,wrap(E[1:].T),lw=1,**ens_props)
      else:
        lines_s  = ax.plot(ii,nan1,"b-", lw=1,label=(str(conf_mult) + '$\sigma$'))
        lines_s += ax.plot(ii,nan1,"b-", lw=1)
        line_mu, = ax.plot(ii,nan1,'b-', lw=2,label='DA mean')

    # Tune plot
    ax.set_ylim( *span(xx) )
    ax.set_xlim(stretch(ii[0],ii[-1],1,int=True))
    ax.set_xlabel('State index')
    ax.set_ylabel('Value')
    ax.legend(loc='upper right')

    text_t = ax.text(0.01, 0.01, format_time(None,None,None),
        transform=ax.transAxes,family='monospace',ha='left')

    # Init visibility (must come after legend):
    if obs_inds is not None:
      line_y.set_visible(False)


    def update(key,E,P):
      nonlocal patch_s
      k,kObs,f_a_u = key

      if conf_mult:
        line_mu.set_ydata(wrap(mu[k]))
        patch_s.remove()
        patch_s = ax.fill_between(ii,
          wrap(mu[k] - conf_mult*sqrt(stats.var[k])),
          wrap(mu[k] + conf_mult*sqrt(stats.var[k])),
          color='b',alpha=0.4)

      else:
        if E is not None:

          for i,line in enumerate(lines_E):
            line.set_ydata(wrap(E[i]))

          if hasattr(stats,'w'):
            w    = stats.w[k]
            wmax = w.max()
            for i,line in enumerate(lines_E):
              line.set_alpha((w[i]/wmax).clip(0.1))
            
        else:
          sigma = mu[k] + conf_mult * sqrt(stats.var[k]) * [[1],[-1]]
          lines_s[0].set_ydata(wrap(sigma[0]))
          lines_s[1].set_ydata(wrap(sigma[1]))
          line_mu   .set_ydata(wrap(mu[k]))

      line_x.set_ydata(wrap(xx[k]))

      text_t.set_text(format_time(k,kObs,stats.HMM.t.tt[k]))

      if 'f' in f_a_u:
        if obs_inds is not None:
          line_y.set_ydata(yy[kObs])
          line_y.set_zorder(5)
          line_y.set_visible(True)
        plot_pause(pause_f)

      if 'a' in f_a_u:
        plot_pause(pause_a)

      if 'u' in f_a_u:
        if obs_inds is not None:
          line_y.set_visible(False)
        plot_pause(pause_u)

      return # end update

    # Finalize init
    update(key0,E,P)

    return update # end init()
  return init # end spatial1d()






