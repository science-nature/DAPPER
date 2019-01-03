# Show that, in the linear context, 
# using the wrong covariance normalization in the EnKF
# does not impact the mean estimate (after an intial transitory period)
# because the sqrt-EnKF anyway converges (becomes) a sqrt-KF.

from common import *

sd0 = seed(5)

##############################
# Setup
##############################
# Model is pathololgical: x(k+1) <-- a*x(k), with a>1.
# Hence, x(k) --> infinity, unless started at exactly 0.
# Therefore, you probably want to set the truth trajectory explicitly to 0
# (or let it be a random walk, or something else).

t = Chronology(dt=1,dkObs=5,T=600,BurnIn=500)

Nx = 5;
Ny = Nx;

jj = equi_spaced_integers(Nx,Ny)
Obs = partial_direct_Obs(Nx,jj)
Obs['noise'] = 1.0

X0 = GaussRV(C=1.0,Nx=Nx)

f = linear_model_setup(1.2*eye(Nx))
Dyn['noise'] = 0.0

#other = {'name': os.path.relpath(__file__,'mods/')}
HMM = HiddenMarkovModel(Dyn,Obs,t,X0)
Dyn,Obs,chrono,X0 = HMM.Dyn, HMM.Obs, HMM.t, HMM.X0


##############################
# DA Config
##############################

@DA_Config
def EnKF_wrong(N,**kwargs):
  """
  The Sqrt EnKF, grabbed from da_methods.py,
  BUT using a different normalization for the covariance estimates than 1/(N-1).

  This "change" must also be implemented for model noise,
  which I have not bothered to do.
  Therefore the results are only relevant with NO MODEL ERROR.
  """
  def assimilator(stats,HMM,xx,yy):
    Dyn,Obs,chrono,X0 = HMM.Dyn, HMM.Obs, HMM.t, HMM.X0

    # Normalization
    #NRM= N-1 # Usual
    NRM= 0.2 # Whatever

    # Init
    E = X0.sample(N)
    stats.assess(0,E=E)

    # Loop
    for k,kObs,t,dt in progbar(chrono.ticker):
      E = Dyn(E,t-dt,dt)
      E = add_noise(E, dt, Dyn.noise, kwargs)

      # Analysis update
      if kObs is not None:
        stats.assess(k,kObs,'f',E=E)

        mu = mean(E,0)
        A  = E - mu

        Eo = Obs(E,t)
        xo = mean(Eo,0)
        Y  = Eo-xo
        dy = yy[kObs] - xo

        d,V= eigh(Y @ Obs.noise.C.inv @ Y.T + NRM*eye(N))
        T  = V@diag(d**(-0.5))@V.T * sqrt(NRM)
        Pw = V@diag(d**(-1.0))@V.T
        w  = dy @ Obs.noise.C.inv @ Y.T @ Pw
        E  = mu + w@A + T@A

      stats.assess(k,kObs,E=E)
  return assimilator

cfgs  = List_of_Configs()
cfgs += ExtKF()
cfgs += EnKF('Sqrt',   10,fnoise_treatm='Sqrt-Core')
cfgs += EnKF_wrong(    10,fnoise_treatm='Sqrt-Core')
#cfgs += EnKF('PertObs',10,fnoise_treatm='Stoch',infl=1.3,LP=False)



##############################
# Run experiment
##############################

# Truth/Obs
xx = zeros((chrono.K+1,Dyn.M))
yy = zeros((chrono.KObs+1,Obs.M))
for k,kObs,t,dt in chrono.ticker:
  # DONT USE MODEL. Use  (below), or comment out entirely.
  pass # xx := 0.
  #xx[k] = xx[k-1] + sqrt(dt)*Dyn.noise.sample(1)  # random walk
  if kObs is not None:
    yy[kObs] = Obs(xx[k],t) + Obs.noise.sample(1)

stats = []
avrgs = []

for ic,config in enumerate(cfgs):
  #config.store_u = True
  #config.liveplotting = False
  seed(sd0+2)

  stats += [ config.assimilate(HMM,xx,yy) ]
  avrgs += [ stats[ic].average_in_time() ]
  #print_averages(config, avrgs[-1])
print_averages(cfgs,avrgs)

plot_time_series(stats[-1])

