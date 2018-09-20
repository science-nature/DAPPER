
from common import *


# Note: the usage of 'dims' in this module is rather unnecessary
# (repleable by [:3], [:], or simply no indexing). 
# It's left for generality, but is probably buggy for systems where m!=3.
m = 3
dims = arange(m)


def LP_setup(
    obs_inds     = None,
    ens_props    = {'alpha':0.3},
    lag          = 40,
    pause_f      = 0.0,
    pause_a      = 0.7,
    weight_alpha = False,
    ):

  def init(stats,key,E,P):

    # Extract data arrays
    k = key[0]
    xx, yy, mu, var = stats.xx, stats.yy, stats.mu, stats.var
    # Alias
    jj = obs_inds
    tt = stats.setup.t.tt

    #####################
    # Set up figure, axes
    #####################
    fig, axs = freshfig(9, (9,5), nrows=m, sharex=True)
    axs[0].set_title("Marginal time series.")
    # Add 3d axes
    fig.subplots_adjust(right=0.45, left=0.1, hspace=0.05)
    bb0 = axs[0].get_position()
    bb2 = axs[2].get_position()
    ax3 = fig.add_axes([0.48, bb2.y0, 0.5, bb0.y1-bb2.y0], projection='3d')
    ax3.set_facecolor('w')
    ax3.set_title("Phase space trajectories")
    # Tune plots
    for s,i in zip("xyz",range(m)):
      set_ilim(ax3   ,i,xx,1.5)
      axs[i].set_ylim(*stretch(*span(xx), 1.01))
      eval("ax3.set_%slabel('%s')"%(s,s), {'ax3':ax3} )
      axs[i].set_ylabel(s)
    axs[2].set_xlabel('t')

    #####################
    # 3d phase space trajectories
    #####################
    if E is None: N = 0
    else:         N = len(E)

    hist_EE = np.full((lag,N,m)    , nan)
    hist_mu = np.full((lag,m)      , nan)
    hist_ss = np.full((lag,2,m)    , nan)
    hist_xx = np.full((lag,m)      , nan)
    hist_yy = np.full((lag,len(jj)), nan)
    hist_tt = np.copy(tt[:lag]) # using nan's would yield bugs.

    PO =  list(jj)==[0,1,2] # PlotObs switch: only if obs operator is id.

    if N : scat_ens =  ax3.scatter(*E.T [dims],s=3 **2, **ens_props)
    else : scat_mu  =  ax3.scatter(*mu[k,dims],s=8 **2, c='b')
    if 1 : scat_x   =  ax3.scatter(*xx[k,dims],s=14**2, c='k',marker=(5, 1))

    if N : tail_ens = [ax3.plot(*hist_EE[:,n,:].T,  **ens_props)[0] for n in range(N)]
    else : tail_mu, =  ax3.plot(*hist_mu       .T, 'b' ,lw=2)
    if 1 : tail_xx, =  ax3.plot(*hist_xx       .T, 'k' ,lw=4)
    if PO: tail_yy, =  ax3.plot(*hist_yy       .T, 'g*',ms=14)

    #####################
    # Marginal time series
    #####################
    line_y  = []
    line_x  = []
    line_mu = []
    lines_s = []
    lines_E = []
    for i, ax in enumerate(axs):
      if N: lines_E += [ax.plot(hist_tt, hist_EE[:,:,i], **ens_props)]
      else:
            line_mu +=  ax.plot(hist_tt, hist_mu  [:,i], 'b')
            lines_s += [ax.plot(hist_tt, hist_ss[:,:,i], 'b--',lw=1)]
      if 1: line_x  +=  ax.plot(hist_tt, hist_xx  [:,i], 'k')
    for i,j in enumerate(jj):
      line_y +=     axs[j].plot(hist_tt, hist_yy  [:,i] ,'g*', ms=10)

    def update(key,E,P):
      nonlocal hist_xx, hist_yy, hist_EE, hist_mu, hist_ss, prev_k, hist_tt, len_roll
      k,kObs,f_a_u = key

      #####################
      # Update rolling data array
      #####################
      a = kObs is not None

      # In case of plot restarting
      if prev_k not in [k, k-1]:
        hist_xx[:] = nan
        hist_yy[:] = nan
        hist_EE[:] = nan
        hist_mu[:] = nan
        hist_ss[:] = nan
        hist_tt[:] = tt[k:k+lag]
        len_roll   = 0
      prev_k = k

      # Abbreviate
      if 1:  _t = tt[k]
      if N:  _E = E  [:,dims]
      else:
            _mu = mu[k,dims]
            _ss = mu[k,dims] + [[1],[-1]]*sqrt(var[k,dims])
      if 1:  _x = xx[k,dims]
      if a:  _y = yy[kObs,:]
      else:  _y = nan*ones(len(jj))

      # Enter current value
      if len_roll<lag: # Grow:
        if 1: hist_tt[len_roll] = _t
        if 1: hist_xx[len_roll] = _x
        if N: hist_EE[len_roll] = _E
        else:
              hist_mu[len_roll] = _mu
              hist_ss[len_roll] = _ss
        if a: hist_yy[len_roll] = _y
        len_roll += 1

      else: # Roll
        if 1: hist_tt = roll_n_sub(hist_tt, _t , -1)
        if 1: hist_xx = roll_n_sub(hist_xx, _x , -1)
        if 1: hist_yy = roll_n_sub(hist_yy, _y , -1)
        if N: hist_EE = roll_n_sub(hist_EE, _E , -1)
        else:
              hist_mu = roll_n_sub(hist_mu, _mu, -1)
              hist_ss = roll_n_sub(hist_ss, _ss, -1)

      #####################
      # Marginal time series
      #####################
      for i, ax in enumerate(axs):
        interval = hist_tt[-1] - hist_tt[0]
        ax.set_xlim(hist_tt[0], hist_tt[-1] + interval/20)
        for n in range(N):      lines_E[i][n] .set_data(hist_tt, hist_EE[:,n,i])
        if not N:
                                line_mu[i]    .set_data(hist_tt, hist_mu[:,i])
                                [lines_s[i][b].set_data(hist_tt, hist_ss[:,b,i]) for b in [0,1]]
        if 1:                   line_x[i]     .set_data(hist_tt, hist_xx[:,i])
      for i,j in enumerate(jj): line_y[i]     .set_data(hist_tt, hist_yy[:,j])

      if 'a' in f_a_u:
        if pause_a: plt.pause(pause_a)

      #####################
      # 3d phase space trajectories
      #####################
      def update_tail(handle,newdata):
        handle.set_data(newdata[:,0],newdata[:,1])
        handle.set_3d_properties(newdata[:,2])

      scat_x._offsets3d = juggle_axes(*tp(xx[k,dims]),'z')
      update_tail(tail_xx, hist_xx)
      if 'tail_yy' in locals():
        update_tail(tail_yy, hist_yy)

      if E is not None:
        scat_ens._offsets3d = juggle_axes(*E[dims].T,'z')

        # Adjust color alpha (for particle filters):
        if weight_alpha:
          colors = scat_ens.get_facecolor()[:,:3]
          if len(colors)==1:
            colors = colors.repeat(N,axis=0)
          alpha = stats.w[k]
          alpha = (alpha/alpha.max()).clip(0.1,0.4)
          scat_ens.set_color(np.hstack([colors, alpha[:,None]]))
          for n in range(N):
            tail_ens[n].set_alpha(alpha[n])
        
        for n in range(N):
          update_tail(tail_ens[n],hist_EE[:,n,:])

      else:
        scat_mu._offsets3d = juggle_axes(*tp(mu[k,dims]),'z')
        update_tail(tail_mu, hist_mu)

      if pause_f: plt.pause(pause_f)

      return # end update()

    # Init
    prev_k   = 0 # Pause/restart monitor
    len_roll = 0 # Growth counter
    update(key,E,P)

    return fig, update # end init()
  return init # end LP_setup()



    
