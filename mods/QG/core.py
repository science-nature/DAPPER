# "QG" Quasi-Geostrophic model.
# 
# Taken from Sakov's enkf-matlab package.
# Model is described in: 
# Sakov, Pavel, and Peter R. Oke.:
#   "A deterministic formulation of the ensemble Kalman filter:
#   an alternative to ensemble square root filters."
#   Tellus A 60.2 (2008): 361-371.
#
# Also see DAPPER/mods/QG/governing_eqn.png,
# and note that:
#  - ψ (psi) is the stream function (i.e. surface elevation)
#  - Doubling time "between 25 and 50"
#

from common import *

#########################
# Model parameters
#########################

# Time settings
dt_internal = 1.25 # CFL ≈ 2.0  # Internal (NO output to DAPPER)
dt          = 4*dt_internal     # DAPPER-facing (outputs to DAPPER)
assert 0==(dt/dt_internal)%1.0, "Must be integer multiple"
# Concerning multiprocessing: the overhead for multiprocessing
# of course decreases as the ratio dt/dt_internal increases.
# But the overhead is minimal already with a ratio of 4.

# These parameters may be interesting to change. 
# In particular, RKH2=2.0e-11 yields a more stable integration,
# and Sakov therefore used it for the ensemble runs (but not the truth).
prms = [
    ["dtout"        , dt         ], # dt registered by DAPPER
    ["dt"           , dt_internal], # dt used internally by Fortran
    ["RKB"          , 0          ], # bottom friction
    ["RKH"          , 0          ], # horizontal friction
    ["RKH2"         , 2.0e-11    ], # biharmonic horizontal friction
    ["F"            , 1600       ], # Froud number
    ["R"            , 1.0e-5     ], # ≈ Rossby number
    ["scheme"       , "'rk4'"    ]  # One of (2ndorder, rk4, dp5)
    ]
# Do NOT change:
prms2 = [
    ["tend"         , 0   ], # Only used by standalone QG
    ["verbose"      , 0   ], # Turn off
    ["rstart"       , 0   ], # Restart switch
    ["restartfname" , "''"], # Read from
    ["outfname"     , "''"]  # Write to
    ]
prms += prms2

# Used list for prms to keep ordering. But it's nice to have a dict too:
prms_dict = {entry[0]: entry[1] for entry in prms}


#########################
# Write parameters to file
#########################
prm_filename = './mods/QG/f90/prms_tmp.txt'
# Create string
prm_txt = """! Parameter file auto generated from python
&parameters"""
for p in prms:
  prm_txt += "\n  " + p[0].ljust(20) + '= ' + str(p[1])
prm_txt += """\n/\n"""
# Write string to file
with open(prm_filename, 'w') as f:
  f.write(prm_txt)


#########################
# Domain management
#########################

# Domain size "hardcoded" in f90/parameters.f90
nx = 129
ny = nx
m  = nx*ny

def square(x):
  psi = x.copy()
  psi = psi.reshape((ny,nx),order='F')
  return psi

def flatten(psi):
  x = psi.ravel(order='F')
  return x

#########################
# Define step functions
#########################

try:
  from mods.QG.f90.py_mod import interface_mod as fortran
except ImportError as err:
  raise Exception("Have you compiled the (fortran) model?" + \
      "\nSee README in folder DAPPER/mods/QG/f90") from err

def step_1(x0, t, dt_):
  """Step a single state vector. Takes care of the copying"""
  assert dt_ == dt
  assert np.isfinite(t)
  assert isinstance(t,float)
  psi = square(x0)
  t   = np.array([t]) # QG is time-indep -- does not matter
  fortran.step(t,psi,prm_filename)
  x   = flatten(psi)
  return x


from tools.utils import multiproc_map
def step(E, t, dt_):
  """Vector and 2D-array (ens) input, with multiproc for ens case."""
  if E.ndim==1:
    return step_1(E,t,dt_)
  if E.ndim==2:
    # Parallelized:
    E = np.array(multiproc_map(step_1, E, t=t, dt_=dt_))
    # Non-parallelized:
    #for n,x in enumerate(E): E[n] = step_1(x,t,dt_)
    return E


#########################
# Free run
#########################
def gen_sample(Len,SpinUp,Spacing):
  sample = zeros((Len,m))
  x = zeros(m)
  n = 0
  for k in progbar(range(Len*Spacing + SpinUp),desc='Simulating'):
    x = step(x,0.0,dt)
    if k>=SpinUp and k%Spacing==0:
      sample[n] = x
      n += 1
  return sample

sample_filename = 'data/samples/QG_samples.npz'
if not os.path.isfile(sample_filename):
  print('Generating a "random" sample with which to start simulations')
  sample = gen_sample(400,500,10)
  np.savez(sample_filename,sample=sample)


#########################
# Plotting
#########################

# Simple animation
#K      = 200
#sample = gen_sample(K,0,10)
#setter = show(sample[0])
#for k in progbar(range(K),desc='Animating'):
  #if k%1==0:
    #setter(sample[k])
    #plt.pause(0.01)

# Although psi is the state variable, q looks cooler.
# q = Del2(psi) - F*psi.
import scipy.ndimage.filters as filters
dx = 1/(nx-1)
def compute_q(psi):
  Lapl = filters.laplace(psi,mode='constant')/dx**2
  return Lapl - prms_dict['F']*psi

# For some reason,
# I cannot get imshow to update inside the DA loop without bugs.
# Its possible to use it after the experiments, though:
#setter = show(xx[0],mu[0],sqrt(var[0]),obs_inds(0))
#for k,kObs,t,dt in progbar(chrono.forecast_range):
#  setter(xx[k],mu[k],var[k],obs_inds(t))
#  plt.suptitle('t = %.1f'%t)
def show_map(stats):
  cmap=plt.cm.viridis
  prep_p = lambda x: square(x)
  prep_q = lambda x: compute_q(square(x))
  prep_m = prep_p
  prep_v = prep_p

  x=stats.xx
  mu=stats.mu.a
  var=stats.var.a

  f, ((ax_p, ax_q), (ax_m, ax_v)) = plt.subplots(2, 2, figsize=(10,10))
  ax_p.set_title('psi')
  ax_q.set_title('q')
  ax_m.set_title('mean estimate of psi')
  ax_v.set_title('std. dev. in psi')
  im_p = ax_p.imshow(prep_p(x[0])        , origin='upper',cmap=cmap)
  im_q = ax_q.imshow(prep_q(x[0])        , origin='upper',cmap=cmap)
  im_m = ax_m.imshow(prep_m(mu[0])       , origin='upper',cmap=cmap)
  im_v = ax_v.imshow(prep_v(sqrt(var[0])), origin='upper',cmap=cmap)
  supt = f.suptitle('step 0/%i'%(x.shape[0]-1))
  #s_yy = ax_p.scatter
  ax_p.set_xlim(0,129)
  ax_p.set_ylim(0,129)
  im_p.set_clim(-35,30)
  im_q.set_clim(-1.2e5,1.2e5)
  im_m.set_clim(-35,30)
  im_v.set_clim(0,1.5)
  plt.pause(0.01)
  print('Type in:\n-"e" to jump to final step,\n-"p" to visualize step-by-step evolution, and "enter" to pause\n-Space to view the first step,\n-"q" to quit.')
  #plt.savefig('QGtestinimap.png')
  
  def setter():
    
    k = getch()
    
    if k == 'p':
      for i in range(1,x.shape[0]-1):
        c = poll_input()
        if c == '\n':
          break
          setter()
        im_p.set_data(prep_p(x[i]))
        im_q.set_data(prep_q(x[i]))
        im_m.set_data(prep_m(mu[i]))
        im_v.set_data(prep_v(sqrt(var[i])))
        supt.set_text('step %i/%i'%(i,x.shape[0]-1))
        #s_yy.set_offsets(array(np.unravel_index(jj,(nx,ny))).T)
        plt.pause(0.01)
      setter()
    
    elif k == 'e' :
      im_p.set_data(prep_p(x[-1]))
      im_q.set_data(prep_q(x[-1]))
      im_m.set_data(prep_m(mu[-1]))
      im_v.set_data(prep_v(sqrt(var[-1])))
      supt.set_text('step %i/%i'%(x.shape[0]-1,x.shape[0]-1))
      plt.pause(0.01)
      setter()
    
    elif k == ' ':
      im_p.set_data(prep_p(x[0]))
      im_q.set_data(prep_q(x[0]))
      im_m.set_data(prep_m(mu[0]))
      im_v.set_data(prep_v(sqrt(var[0])))
      supt.set_text('step %i/%i'%(0,x.shape[0]-1))
      #s_yy.set_offsets(array(np.unravel_index(jj,(nx,ny))).T)
      plt.pause(0.01)
      setter()
    
    elif k == 'q':
      #plt.close()
      return None

    else :
      setter()

  return setter()

def show_error_map(stats):
  error=stats.err.a
  spread=stats.mad.a
  f2, (ax1,ax2) = plt.subplots(1, 2, figsize=(10,6))
  supt = f2.suptitle('Time evolution of RMSE\nstep 0/%i'%(error.shape[0]))
  im_err = ax1.imshow(error[0].reshape((129,129),order='F'),origin='upper',cmap=plt.cm.viridis)
  im_sp = ax2.imshow(spread[0].reshape((129,129),order='F'),origin='upper',cmap=plt.cm.viridis)
  ax1.set_title('RMSE')
  ax2.set_title('Spread')
  plt.pause(0.01)
  print('Type in:\n-"e" to jump to final step,\n-"p" to visualize step-by-step evolution, and "enter" to pause\n-Space to view the first step,\n-"q" to quit.')
  #plt.savefig('QGtestinimap.png')
  def setter():

    k=getch()
    if k == 'p':
      for (i,(r,s)) in enumerate(zip(error[1:],spread[1:])):
        c = poll_input()
        if c == '\n':
          break
          setter()
        im_err.set_data(r.reshape((129,129),order='F'))
        im_sp.set_data(s.reshape((129,129),order='F'))
        supt.set_text('Time evolution of RMSE\nstep %i/%i'%(i,error.shape[0]))
        plt.pause(0.01)
      setter()
    elif k == 'e' :
      im_err.set_data(error[-1].reshape((129,129),order='F'))
      im_sp.set_data(spread[-1].reshape((129,129),order='F'))
      supt.set_text('Time evolution of RMSE\nstep %i/%i'%(error.shape[0],error.shape[0]))
      plt.pause(0.01)
      setter()
    elif k == ' ':
      im_err.set_data(error[0].reshape((129,129),order='F'))
      im_sp.set_data(spread[0].reshape((129,129),order='F'))
      supt.set_text('Time evolution of RMSE\nstep 0/%i'%(error.shape[0]))
      plt.pause(0.01)
      setter()
    elif k == 'q':
      return None
    else:
      setter()

  return setter()

def plot_err(stats,supt=None):
    f,ax = plt.subplots()
    ax.plot(stats.rmse.a,label='RMSE')
    ax.plot(mean(stats.mad.a,axis=1),label='Spread')
    ax.legend()
    ax.set_title('Time through evol. RMSE & Spread')
    ax.set_ylabel('RMSE/Spread')
    if type(supt) == str:
      f.suptitle(supt)

  
#Show tracking of the observations on the map
#Satisfying number of diags for QG constraints are multiples of 2 (divisors of 128)
#Suggested params: diags=16,angle=-1, lighten_factor=3 or diags=8,angle=[-1,1],lighten_factor=3
def show_tracking(m=16641,diags=16,angle=-1,lighten_factor=3):

  jjs = generate_diags(**locals())

  d = int(m**0.5)

  ll=[0]*m

  for i in jjs:
    ll[i]=1
  f , ax= plt.subplots()
  ax.imshow(reshape(ll,(d,d)),cmap=plt.cm.YlGnBu)
  ax.set_title('%i points tracked'%sum(ll))

def generate_diags(m=16641,diags=16,angle=-1,lighten_factor=3):
  if not hasattr(angle,'__iter__'):
    angle = [angle]
  
  jjs = [] 

  d = int(m**0.5)
  for a in angle:
    points = range(m)
    diagos = range(-(d//2)-a*(d//2),d+d//2+1-a*(d//2),(2*(d-1)) // diags)
    tracks = filter(lambda x : x%d - a * (x//d) in diagos, points)
    tracks = filter(lambda x : ( x%d + x%d - a * (x//d) )%lighten_factor == 0,tracks)
    jjs+=tracks

  return array(jjs)



