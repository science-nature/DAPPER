# Benchmark EnKF PertObs where the "denomiator" uses the true R,
# or an estimate made from the observation perturbations.
# Conclusion: true R method is far superior for L95.
#
# To run script, first insert this snippet in da_methods.py:EnKF_analysis():
#     if 'PertObs' in upd_a:
#       if 'Re aug' in upd_a:
#         # here K = AY'/[YY'+DD']
#         D  = center(hnoise.sample(N))
#         Re = D.T @ D
#         C  = Y.T @ Y + Re
#         YC = Y @ tinv(C)
#         KG = A.T @ YC
#       elif 'Re pinv' in upd_a:
#         # here K = A pinv(Y+D) (lin-reg: x --> h(x)+e)
#         D  = center(hnoise.sample(N))
#         YD = Y + D
#         KG = A.T @ tinv(YD.T)
#       elif 'Re sum' in upd_a:
#         # here K = A(Y+D)'/[(Y+D)(Y+D)'] (same as pinv)
#         D  = center(hnoise.sample(N))
#         YD = Y + D
#         C  = YD.T @ YD
#         YC = YD @ tinv(C)
#         KG = A.T @ YC
#       elif 'Re cross0' in upd_a:
#         # here K = AY'/[(Y+D)(Y+D)']
#         D  = center(hnoise.sample(N))
#         YD = Y + D
#         C  = YD.T @ YD
#         YC = Y @ tinv(C)
#         KG = A.T @ YC
#       elif 'Re cross1' in upd_a:
#         # here K = A(Y+D)'/[YY'+DD']
#         D  = center(hnoise.sample(N))
#         Re = D.T @ D
#         C  = Y.T @ Y + Re
#         YD = Y + D
#         YC = YD @ tinv(C)
#         KG = A.T @ YC
#       elif 'Re proj' in upd_a:
#         # here K = AY'/[YY'+ P DD' P], where P is Proj mat onto col(Y))
#         D  = center(hnoise.sample(N))
#         Re = D.T @ D
#         P  = Y.T @ tinv(Y.T)
#         C  = Y.T @ Y + P@Re@P
#         YC = Y @ tinv(C)
#         KG = A.T @ YC
#       else: # i.e. just use the true R
#         # here K = AY'/[YY' + nR]
#         C  = Y.T @ Y + R.full*N1
#         D  = center(hnoise.sample(N))
#         YC = Y @ tinv(C)
#         KG = A.T @ YC
#         HK = Y.T @ YC
#       dE = (KG @ ( y + D - hE ).T).T
#       E  = E + dE
#     elif 'Sqrt' in upd_a:
#       [...]

from common import *

sd0 = seed_init(4)

from mods.Lorenz95.sak08 import setup
setup.t = Chronology(dt=0.05, dkObs=1, T=100, BurnIn=5)

CtrlVar = sys.argv[1] # command-line argument #1

assert CtrlVar == 'N' # Ensemble size.
xticks = [16, 20, 25, 30, 40, 70, 100, 300, 1000]
xticks = array(xticks).repeat(8)

xticks, save_path, rep_inds = distribute(__file__,sys.argv,xticks,CtrlVar)


##############################
# DA Configurations
##############################
cfgs  = List_of_Configs()
cfgs += OptInterp()

for UPD_A in ['PertObs'+s for s in ['',' Re aug',' Re pinv',' Re sum',' Re cross0',' Re cross1',' Re proj']]:
  for INFL in [1.00, 1.02, 1.05, 1.10, 1.20, 1.30, 1.50, 1.70, 3.0]:
    cfgs += EnKF(UPD_A,N='?',infl=INFL)

def adjust_cfg(C,variable,X):
  assert variable=='N'
  if getattr(C,'N',None)=='?': C = C.update_settings(N=X)
  return C


##############################
# Assimilate
##############################
avrgs = np.empty((len(xticks),1,len(cfgs)),dict)

for iX,(X,iR) in enumerate(zip(xticks,rep_inds)):
  with coloring(): print('\n'+"xticks[",iX,'/',len(xticks)-1,"] ",CtrlVar,': ',X,sep="")

  seed(sd0+iR)
  xx,yy = simulate(setup)

  for iC,C in enumerate(cfgs):
    C = adjust_cfg(C,CtrlVar,X)
    seed(sd0+iR)
    avrgs[iX,0,iC] = C.assimilate(setup,xx,yy).average_in_time()

  print_averages(cfgs,avrgs[iX,0])

# Results saved in the format below is supported by DAPPER's ResultsTable, whose main purpose
# is to collect result data from parallelized (or otherwise independent) experiments.
np.savez(save_path,
    avrgs      = avrgs,            # 3D array of dicts, whose fields are the averages.
    xlabel     = CtrlVar,          # The control variable tag (string).
    xticks     = xticks,           # xticks (array).
    tuning_tag = 'infl',           # Tag to search for within labels.
    labels     = cfgs.gen_names()) # List of strings.


##############################
# Results load & presentation
##############################
if 'WORKER' in sys.argv: sys.exit(0) # quit if script is running as worker.

R = ResultsTable(save_path)

##
R = ResultsTable('data/Stoch_iEnS/bench_R_ens/P2720L/N_run2') # Contains the methods I first tested
R          .load('data/Stoch_iEnS/bench_R_ens/P2720L/N_run3') # Added methods: [' Re cross1',' Re proj']

# Print averages of a given field.
# The "subcolumns" show the number of repetitions, crashes and the 1-sigma conf.
with coloring(): print("Averages over experiment repetition:")
R.print_mean_field('rmse_a',1,1,cols=slice(0,2))

# Separate out the baseline methods from the rest
BaseLineMethods = R.split(lambda x: x in ['Climatology', 'OptInterp', 'Var3D','ExtKF'])

# Plot
fig, ax = plt.subplots()
ax, ax_, lhs = R.plot_1d_minz('rmse_a',)

# Baseline plot
plt.sca(ax)
BaseLineMethods.plot_1d('rmse_a',color='k')

# Adjust plot
if R.xlabel=='N':
  ax.loglog()
  ax.legend()
  ax.grid(True,'minor')
  xt = [15,20,30,40,50,70,100, 200, 1000]
  yt = [0.1, 0.2, 0.5, 1, 2, 5]
  ax_.set_xticks(xt); ax.set_xticklabels(xt)
  ax .set_xticks(xt); ax.set_xticklabels(xt)
  ax .set_yticks(yt); ax.set_yticklabels(yt)


##

