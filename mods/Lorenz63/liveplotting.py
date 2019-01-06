
from common import *


#TODO
# Note: the usage of 'dims' in this module is rather unnecessary
# (repleable by [:3], [:], or simply no indexing). 
# It's left for generality, but is probably buggy for systems where Nx!=3.
Nx = 3
dims = arange(Nx)


def phase3D(
    obs_inds     = [],
    lag          = 100,
    ens_props    = {'alpha':0.3},
    zoom_3d      = 1.5,
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
    tt = stats.HMM.t.tt

    #####################
    # Set up figure, axes
    #####################
    fig = plt.figure(8, figsize=(5,5))
    # Add 3d axes
    # ax3 = fig.add_axes(projection='3d')
    ax3 = plt.subplot(111, projection='3d')
    ax3.set_facecolor('w')
    ax3.set_title("Phase space trajectories")
    # Tune plots
    for s,i in zip("xyz",range(Nx)):
      set_ilim(ax3,i,  *stretch(*span(xx[:,i]), 1/zoom_3d) )
      eval("ax3.set_%slabel('%s')"%(s,s), {'ax3':ax3} )

    #####################
    # 3d phase space trajectories
    #####################
    if E is None: N = 0
    else:         N = len(E)

    hist_EE = np.full((lag,N,Nx)   , nan)
    hist_mu = np.full((lag,Nx)     , nan)
    hist_ss = np.full((lag,2,Nx)   , nan)
    hist_xx = np.full((lag,Nx)     , nan)
    hist_yy = np.full((lag,len(jj)), nan)
    hist_tt = np.copy(tt[:lag]) # using nan's would yield bugs.

    PO = list(jj)==[0,1,2] # PlotObs switch: only if obs operator is id.

    if N : scat_ens =  ax3.scatter(*E.T [dims],s=3 **2, **ens_props)
    else : scat_mu  =  ax3.scatter(*mu[k,dims],s=8 **2, c='b')
    if 1 : scat_x   =  ax3.scatter(*xx[k,dims],s=14**2, c='k',marker=(5, 1))

    if N : tail_ens = [ax3.plot(*hist_EE[:,n,:].T,  **ens_props)[0] for n in range(N)]
    else : tail_mu, =  ax3.plot(*hist_mu       .T, 'b' ,lw=2)
    if 1 : tail_xx, =  ax3.plot(*hist_xx       .T, 'k' ,lw=4)
    if PO: tail_yy, =  ax3.plot(*hist_yy       .T, 'g*',ms=14)


    def update(key,E,P):
      nonlocal hist_xx, hist_yy, hist_EE, hist_mu, hist_ss, prev_k, hist_tt, len_roll
      k,kObs,f_a_u = key

      #####################
      # Update rolling data array
      #####################
      a = PO and kObs is not None

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



    
