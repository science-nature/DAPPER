##
from common import *
OD = OrderedDict # alias
DFLT   = OD(ls='-',lw=2,marker='o',ms=6,markeredgewidth=0)
hollow = OD(lw=0,markerfacecolor='none',ms=10,markeredgewidth=1.5)

##
def style(ax):
  for line in ax.get_lines():
    label = line.get_label() 
    if label != 'dontstyle':
      props = compile_style(label)
      line.set(**props)

def compile_style(label):
  """Compile all matching styles for a given label."""
  style_dict = DFLT.copy()       # Init
  for pattern, style in STYLES:  # Loop through style rows
    if re.search(pattern,label): # Match label
      style_dict.update(style)   # Update/overwrite
  return style_dict

STYLES = [
    # Label                         , Style props
    ('(Climatology|OptInterp)'      , OD(c=blend_rgb('k',0.7),ls='-',lw=1.0,ms=0) ),
    ('EnKF'                         , OD(ls=':',ms=0                ) ),
    ('upd_a:Sqrt'                   , OD(c='mlr',marker='s'         ) ), # [.83,.56,.74]
    ('upd_a:PertObs'                , OD(c='mlb',marker='X'         ) ), # [.64,.66,.83]
    ('upd_a:Sqrt.*MDA:1'            , OD(c='mly',marker='^',zorder=0) ), # [.28,.69,.64]
    ('upd_a:PertObs.*MDA:1'         , OD(c='mlg',marker='o',zorder=0) ), # [.89,.31,.13]
    ('upd_a:(Order1|DEnKF)'         , OD(c='mlv',marker='v'         ) ), # 
    ('upd_a:Order1.*MDA:1'          , OD(c='mlv',marker='P',zorder=0) ), #
    ('nIter:10'                     , hollow                          ),
    ]


def edit_entries(ax, conversion_table):
  "Edit legend entries: rename, order, turn on/off."
  # Initialise legend entries with blanks
  invisible_line, = ax.plot(np.NaN, np.NaN, '-', color='none', label='placeholder')
  labels   = np.full((99,99),''            , dtype=object)
  handles  = np.full((99,99),invisible_line, dtype=object)
  lines_with_actual_labels = {}
  # For cropping the completed table
  I = 0
  J = 0
  # Loop through current legend, change entries according to table
  for handle, label in zip(*ax.get_legend_handles_labels()):
    for pattern, sub, (i,j) in conversion_table:
      if all(re.search(x, label) for x in list(pattern)):
        if handles[i,j] is not invisible_line:
          print("\nWarning: some entries from the "+\
              "conversion_table use same coordinates.\n")
        handles[i,j] = handle
        labels [i,j] = sub
        I = max(I,i)
        J = max(J,j)
        lines_with_actual_labels[label] = handle
        break
  ncol = J+1
  # Crop and make into lists
  handles = list( handles[:I+1, :J+1].ravel('F') )
  labels  = list( labels [:I+1, :J+1].ravel('F') )
  return handles, labels, ncol, lines_with_actual_labels

FRMT = [
    # pattern                        , sub                , table_coord
    # ( ['Climatology'               ] , 'System var.'      , (0 , 0) ) ,
    # ( ['OptInterp'                 ] , 'Opt. Interp.'     , (0 , 0) ) ,
    #
    # ( ['Order1' ,'MDA:0','nIter:3' ] , 'Order1 iEnKS 3'   , (0 , 0) ),
    # ( ['Order1' ,'MDA:1','nIter:3' ] , 'Order1 MDA 3'     , (0 , 0) ),
    #
    ( ['EnKF.*PertObs'             ] , 'Stoch. EnKF'      , (0 , 0) ),
    ( ['EnKF.*Sqrt'                ] , 'Determ. EnKF'     , (0 , 0) ),
    ( ['EnKF.*DEnKF'               ] , 'Order1 EnKF'      , (0 , 0) ),
    #
    ( ['PertObs','MDA:1','nIter:3' ] , 'Stoch. MDA 3'     , (0 , 0) ),
    ( ['PertObs','MDA:0','nIter:3' ] , 'EnRML 3'          , (1 , 0) ),
    ( ['Sqrt'   ,'MDA:1','nIter:3' ] , 'Determ. MDA 3'    , (2 , 0) ),
    ( ['Sqrt'   ,'MDA:0','nIter:3' ] , 'iEnKS 3'          , (3 , 0) ),
    #
    ( ['PertObs','MDA:1','nIter:10'] , '10'               , (0 , 1) ),
    ( ['PertObs','MDA:0','nIter:10'] , '10'               , (1 , 1) ),
    ( ['Sqrt'   ,'MDA:1','nIter:10'] , '10'               , (2 , 1) ),
    ( ['Sqrt'   ,'MDA:0','nIter:10'] , '10'               , (3 , 1) ),
    ]



##
R = ResultsTable("data/Stoch_iEnS/bench_L95/P2720L/N_run8"); LAG=2           # dkObs=4
# R = ResultsTable("data/Stoch_iEnS/bench_L95/P2724L/N_run2"); LAG=1            # dkObs=8

# with coloring(): print("Averages over experiment repetition:")
# R.print_mean_field('rmse_a',1,1,cols=slice(0,2))

BaseLineMethods = R.split(['Climatology', 'OptInterp', 'Var3D','ExtKF'])
BaseLineMethods.rm('Var3D')

R.rm("rot:0")
R.rm("upd_a:(Order1|DEnKF)")
# R.rm("upd_a:Order1.*MDA:0")
R.rm("EnKF")


## Plotting
plt.style.use('Stoch_iEnS/paper.mplstyle')

# Make Broken axis using figure with 2 subplots.
fig, axs = plt.subplots(nrows=2,sharex=True,
    gridspec_kw={'height_ratios':[1, 4], 'hspace':0.04})
ax0, ax1 = axs
ax0.set_ylim(1,5)
ax0.set_yticks([2,3,4,5])
ax1.set_ylim(0.25,1)
# Remove border between axes
ax0.spines['bottom'].set_visible(False)
ax1.spines['top']   .set_visible(False)
ax0.xaxis.tick_top()
ax1.xaxis.tick_bottom()
ax0.tick_params(labeltop='off')
# Add wavy/squiggly line to show axis is broken
ax0.set_xlim(12,105)
xs =                    LogSp(*ax0.get_xlim(),401)
ys = 0.8 + 0.1*cos(5*linspace(*ax0.get_xlim(),401))
ax0.plot(xs, ys, 'k-',lw=2, alpha=0.8, clip_on=False,label='dontstyle')

# Plot
field = 'rmse_a'
unique, _, tuning_vals, fieldvals = R.select_optimal(field)
for ax in axs:
  for group_name, tunings, vals in zip(unique, tuning_vals, fieldvals):
    if ax is ax0: # Alter out-of-bounds data to avoid markers halfway appearing along lower edge
      vals = copy(vals); vals[vals<1] = 0.2  # NB slightly changes the angle
    ax.plot(R.xticks, vals, label=group_name)
  for vals, label in zip(BaseLineMethods.mean_field(field)[0], BaseLineMethods.labels):
    ax.plot(R.xticks, vals, label=label)
  style(ax)

# Legend
# if 'checkmarks' not in locals(): checkmarks = []
# checkmarks += [toggle_lines()];
handles, labels, ncol, hh = edit_entries(ax0,FRMT)
leg = ax0.legend(handles, labels, ncol=ncol,
    markerfirst   = 0,             columnspacing = 0.9,
    loc           = 'upper right', handletextpad = 0.2,
    framealpha    = 1.0,           handlelength  = 1.5,
    borderaxespad = 0.1,
    )
ax0.set_zorder(9) # raise legend above ax1

# Axes props
# ax0.set_title(R.meta)
ax1.grid(True,'minor')
ax1.set_xscale('log')
xl = ax1.set_xlabel("Ensemble size ($N$)")
fig.text(0.04, 0.5, "Average RMS error", # "analysis RMSE"
    va='center', rotation='vertical',fontsize=xl.get_fontsize())


# xticks
# ax1.xaxis.set_minor_formatter(mpl.ticker.LogFormatter())
if R.xlabel=='N':
  xt = [20, 30, 40, 60, 80, 100]
elif R.xlabel=='nIter':
  xt = [1, 2, 4, 8, 16, 32]
# xt = R.xticks
ax1.set_xticks(xt); ax1.set_xticklabels(xt)
# Log scaling and associated xtick problems:
# stackoverflow.com/a/14530857/38281
# stackoverflow.com/a/42886695/38281
ax1.xaxis.set_major_formatter(mpl.ticker.ScalarFormatter())
ax1.xaxis.set_minor_formatter(mpl.ticker.NullFormatter())



## Populate plot frame-by-frame for animation effects.
savefig_n.index = 1
def save():
  # savefig_n('data/Stoch_iEnS/figs/prez/'+R.xlabel+'_dkObs_%d'%R.meta['dkObs']+'_i')
  pass

# NB: Don't abbreviate this. I tried a lot (and it was a waste of time):
# - unifying dicts by zipping over all axes (fails coz all axes don't contain all lines)
# - findobj (fails coz complicates getting _legmarker and text)
# - looping over axes (gets messier than this)
objects = {}
objects[ax0] = {line.get_label(): line                          for line in           ax0.get_lines()}
objects[ax1] = {line.get_label(): line                          for line in           ax1.get_lines()}
objects[leg] = {line.get_label(): [line, line._legmarker, text] for line, text in zip(leg.get_lines(), leg.get_texts())}
def tog(label):
  for key in objects:
    toggle_viz(objects[key][label])

  # Text seems sometimes not responsive to set_visible() calls.
  # Hack: plot & rm an arbitrary line to flush the calls.
  ax1.plot(ax1.get_xlim(), ax1.get_ylim())[0].remove()

# Hide all elements
for h in hh: tog(h)

# Animate
save()
tog('iEnKS N:? upd_a:PertObs       Lag:%s MDA:0 nIter:3 '%LAG); save()
tog('iEnKS N:? upd_a:PertObs       Lag:%s MDA:0 nIter:10'%LAG); save()
tog('iEnKS N:? upd_a:PertObs       Lag:%s MDA:1 nIter:3 '%LAG); save()
tog('iEnKS N:? upd_a:PertObs       Lag:%s MDA:1 nIter:10'%LAG); save()
tog('iEnKS N:? upd_a:Sqrt    rot:1 Lag:%s MDA:0 nIter:3 '%LAG); save()
tog('iEnKS N:? upd_a:Sqrt    rot:1 Lag:%s MDA:0 nIter:10'%LAG); save()
tog('iEnKS N:? upd_a:Sqrt    rot:1 Lag:%s MDA:1 nIter:3 '%LAG); save()
tog('iEnKS N:? upd_a:Sqrt    rot:1 Lag:%s MDA:1 nIter:10'%LAG); save()

##


