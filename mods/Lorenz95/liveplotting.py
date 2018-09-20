
from common import *

def LP_setup(
    obs_inds     = None,
    ens_props    = {'color': 0.7*RGBs['w'],'alpha':0.5},
    pause_a      = 0.5,
    pause_f      = 0,
    weight_alpha = False,
    conf_mult    = 2,
    conf_patch   = False,
    periodic     = True,
    ):

  def init(stats,key,E,P):

    # Extract data arrays
    xx, yy, mu = stats.xx, stats.yy, stats.mu
    k = key[0]
    tt = stats.setup.t.tt

    K, m = xx.shape
    ii,wrap = setup_wrapping(m,periodic)

    # Set up figure, axes
    fig, ax = freshfig(9, (8,5))
    fig.suptitle("Amplitude plot")
    fig.subplots_adjust(bottom=0.14)

    line_x, = ax.plot(ii,wrap(xx[k]),'k-',lw=3,label='Truth')
    if obs_inds is not None:
      line_y, = ax.plot(obs_inds, 0.0*obs_inds,'g*',ms=5,label='Obs')

    if conf_patch:
      line_mu, = ax.plot(ii,wrap(mu[k]),'b-', lw=2,label='DA mean')
      patch_s  = ax.fill_between(ii,
          wrap(mu[k] - conf_mult*sqrt(stats.var[k])),
          wrap(mu[k] + conf_mult*sqrt(stats.var[k])),
          color='b',alpha=0.4,label=(str(conf_mult) + '$\sigma$'))

    else:
      if E is not None:
        lines_E  = ax.plot(ii,wrap(E[0] .T),lw=1,**ens_props,label='Ensemble')
        lines_E += ax.plot(ii,wrap(E[1:].T),lw=1,**ens_props)
      else:

        sigma    = mu[k] + conf_mult * sqrt(stats.var[k]) * [[1],[-1]]
        lines_s  = ax.plot(ii,wrap(sigma[0]),"b-", lw=1,label=(str(conf_mult) + '$\sigma$'))
        lines_s += ax.plot(ii,wrap(sigma[1]),"b-", lw=1)
        line_mu, = ax.plot(ii,wrap(mu[k]   ),'b-', lw=2,label='DA mean')

    text_f = ax.text(0.014, 0.96,'Forecast',transform=ax.transAxes, color='b')
    text_a = ax.text(0.014, 0.96,'Analysis',transform=ax.transAxes, color='r')

    # Time info
    text_t = ax.text(0.01, 0.01, format_time(None,None,None),
        transform=ax.transAxes,family='monospace',ha='left')

    ax.legend(loc='upper right')

    # Init visibility (must come after legend):
    text_a.set_visible(False)
    text_f.set_visible(False)
    if obs_inds is not None:
      line_y.set_visible(False)

    ax.set_ylim( *span(xx) )
    ax.set_xlim(stretch(ii[0],ii[-1],1,int=True))
    ax.set_xlabel('State index')
    ax.set_ylabel('Value')


    def update(key,E,P):
      nonlocal patch_s

      k,kObs,f_a_u = key

      if 'a' in f_a_u:
        text_f.set_visible(False)
        if pause_a:
          plt.pause(pause_a);
          text_a.set_visible(True)
          plt.pause(pause_a)

      if conf_patch:
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

          if weight_alpha:
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

      text_t.set_text(format_time(k,kObs,tt[k]))

      if pause_f: plt.pause(pause_f)

      if 'f' in f_a_u:
        if obs_inds is not None:
          line_y.set_ydata(yy[kObs])
          line_y.set_visible(True)

      if 'a' in f_a_u:
        if pause_a:
          plt.pause(pause_a)
          text_a.set_visible(False)
        if obs_inds is not None:
          line_y.set_visible(False)
          line_y.set_zorder(5)
        if pause_a:
          plt.pause(pause_a)
          text_f.set_visible(True)
          plt.pause(pause_a)

      return # end update
    return fig, update # end init
  return init # end LP_setup



    

