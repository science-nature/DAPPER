from common import *

sd0 = seed(3)

##############################
# DA Configurations
##############################
from mods.Lorenz95.sak08 import setup
# setup.t.BurnIn = 5
setup.t.T = 4**3.5
setup.t.dkObs = 3

# Get experiment control variable (CtrlVar) from arguments
CtrlVar = sys.argv[1]
# Set range of experimental settings
if CtrlVar == 'N': # ens size
  # xticks = [30]
  xticks = [15, 17, 20, 25, 30, 40, 60, 200]

xticks = array(xticks).repeat(3)

# Parallelization and save-path setup
xticks, save_path, iiRep = distribute(__file__,sys.argv,xticks,CtrlVar)


##############################
# Configs
##############################
cfgs  = List_of_Configs()

# BASELINES
cfgs += Climatology()
cfgs += OptInterp()
cfgs += Var3D()

nIter  = 10
Lag    = 3   # TODO: 8 is optimal ?
infls  = [1.01, 1.02, 1.04, 1.06, 1.10, 1.2, 1.4]
# infls  = [1.01, 1.03, 1.07, 1.25]

for infl in infls:
  cfgs +=   EnKF ('PertObs', N='?', infl=infl,)
  for rot in [False, True]:
    cfgs += EnKF ('Sqrt'   , N='?', infl=infl, rot=rot, )
    cfgs += iEnKS('Sqrt'   , N='?', infl=infl, rot=rot, Lag=Lag, nIter=nIter)
  cfgs +=   iEnKS('EnRML'  , N='?', infl=infl,          Lag=Lag, nIter=nIter)
  cfgs +=   iEnKS('ES-MDA' , N='?', infl=infl,          Lag=Lag, nIter=nIter)


##############################
# Assimilate
##############################
avrgs = np.empty((len(xticks),1,len(cfgs)),dict)
stats = np.empty_like(avrgs)

for iX,(X,iR) in enumerate(zip(xticks,iiRep)):
  with coloring(): print('\n'+"xticks[",iX,'/',len(xticks)-1,"] ",CtrlVar,': ',X,sep="")
  # setattr(setup.t,CtrlVar,X)

  sd    = seed(sd0 + iR)
  xx,yy = simulate(setup)

  for iC,C in enumerate(cfgs):

    if CtrlVar=='N' and hasattr(C,'N'):
      C = C.update_settings(N=X)

    seed(sd)

    stat = C.assimilate(setup,xx,yy)
    avrg = stat.average_in_time()

    stats[iX,0,iC] = stat
    avrgs[iX,0,iC] = avrg
  print_averages(cfgs, avrgs[iX,0],statkeys=['rmse_a','rmv_a','infl'])

#plot_time_series(stats[-1])

np.savez(save_path,
    avrgs      = avrgs,
    xlabel     = CtrlVar,
    xticks     = xticks,
    tuning_tag = 'infl', 
    labels     = cfgs.gen_names(do_tab=True),
    meta       = {'dkObs':setup.t.dkObs})


##############################
# Results load & presentation
##############################
if 'WORKER' in sys.argv: sys.exit(0) # quit if script is running as worker.

R = ResultsTable(save_path)
# R = ResultsTable("data/Stoch_iEnS/bench_L95/P2720L/N_run2")

# with coloring(): print("Averages over experiment repetition:")
# R.print_mean_field('rmse_a',1,1,cols=slice(0,2))

BaseLineMethods = R.split(['Climatology', 'OptInterp', 'Var3D','ExtKF'])

# Plot
fig, ax = plt.subplots()
# R.plot_1d('rmse_a',)
ax, ax_, lhs = R.plot_1d_minz('rmse_a',)
ax.legend()
plt.sca(ax)
# if 'checkmarks' not in locals(): checkmarks = []
# checkmarks += [toggle_lines()];
BaseLineMethods.plot_1d('rmse_a',color='k')

# Adjust plot
ax.set_yscale('log')
ax.set_xscale('log')
ax.grid(True,'minor')
xt = R.xticks
yt = [0.1, 0.2, 0.5, 1, 2, 5]
ax.set_xticks(xt); ax.set_xticklabels(xt)
ax.set_yticks(yt); ax.set_yticklabels(yt)


##

