# "QG" Quasi-Geostrophic model.
# 
# Taken from Sakov's enkf-matlab package.
# Model is described in: 
# Sakov, Pavel, and Peter R. Oke.:
#   "A deterministic formulation of the ensemble Kalman filter:
#   an alternative to ensemble square root filters."
#   Tellus A 60.2 (2008): 361-371.
#
# Also see:
#  - DAPPER/mods/QG/governing_eqn.png
#  - DAPPER/mods/QG/demo.py
# and note that:
#  - ψ (psi) is the stream function (i.e. surface elevation)
#  - Doubling time "between 25 and 50"

from common import *

#########################
# Model parameters
#########################

# Time settings. Can be changed, but this must be done here.
dt_internal = 1.25 # CFL ≈ 2.0  # Internal (NO output to DAPPER)
dt          = 4*dt_internal     # DAPPER-facing (outputs to DAPPER)
assert 0==(dt/dt_internal)%1.0, "Must be integer multiple"

# These parameters may be interesting to change. 
# In particular, RKH2=2.0e-11 yields a more stable integration,
# and Sakov therefore used it for the ensemble runs (but not the truth).
prms = [
    ["dtout"        , dt         ], # dt registered by DAPPER
    ["dt"           , dt_internal], # dt used internally by Fortran
    ["RKB"          , 0          ], # bottom friction
    ["RKH"          , 0          ], # horizontal friction
    ["RKH2"         , 2.0e-12    ], # biharmonic horizontal friction
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

# Used list for prms above to keep ordering. But a dict's nice:
prms_dict = {entry[0]: entry[1] for entry in prms}


#########################
# Write Fortran namelist from prms
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
# Domain size "hardcoded" in f90/parameters.f90.

# "Physical" domain length -- copied from f90/parameters.f90.
# It appears that only square domains generate any dynamics.
NX1 = 2
NY1 = 2
# Resolution level -- copied MREFIN from parameters.f90
res = 7
# Grid lengths.
nx = NX1 * 2 ** (res - 1) + 1 # (axis=1)
ny = NY1 * 2 ** (res - 1) + 1 # (axis=0)
# Actually, the BCs are psi = nabla psi = nabla^2 psi = 0,
# => psi should always be zero on the boundries.
# => it'd be safer to rm boundries from the DA state vector,
#    yielding ndim(state)=(nx-2)*(ny-2), but this is not done here.

# Fortran model (e.g. f90/interface.f90) requires orientation: X[ix,iy].
shape = (nx, ny)
# Passing arrays to/from Fortran requries that flags['F_CONTIGUOUS']==True.
order = 'F'
def    py2f(x):   return x.reshape(            shape      ,order=order)
def    f2py(X):   return X.flatten(                        order=order)
# However, FOR PRINTING/PLOTTING PURPOSES, the y-axis should be vertical
# [imshow(mat) uses the same orientation as print(mat)].
def  square(x):   return x.reshape(            shape[::-1])
def ind2sub(ind): return np.unravel_index(ind, shape[::-1])


#########################
# Define step functions
#########################

try:
  from mods.QG.f90.py_mod import interface_mod as fortran
except ImportError as err:
  raise Exception("Have you compiled the (Fortran) model?" + \
      "\nSee README in folder DAPPER/mods/QG/f90") from err

def step_1(x0, t, dt_):
  """Step a single state vector."""
  assert dt_ == dt                 # dt read through parm file, so don't change it here.
  assert isinstance(t,float)       # Coz Fortran is a typed language.
  assert np.isfinite(t)            # QG is autonomous, but nan/inf wont work.
  psi = py2f(x0.copy())            # Copy coz Fortran will modify in-place.
  fortran.step(t,psi,prm_filename) # Call Fortran model.
  x = f2py(psi)                    # Flattening
  return x


# Model multiprocessing switch
mp = True

def step(E, t, dt_):
  """Vector and 2D-array (ens) input, with multiproc for ens case."""
  if E.ndim==1:
    return step_1(E,t,dt_)
  if E.ndim==2:
    if mp: # PARALLELIZED:
      # Note: the relative overhead for parallelization decreases
      # as the ratio dt/dt_internal increases.
      # But the overhead is already negligible with a ratio of 4.
      E = np.array(multiproc_map(step_1, E, t=t, dt_=dt_))
    else: # NON-PARALLELIZED:
      for n,x in enumerate(E): E[n] = step_1(x,t,dt_)
    return E


#########################
# Free run
#########################

def gen_sample(nSamples,SpinUp,Spacing):
  simulator = make_recursive(step,prog="Simul.")
  K         = SpinUp + nSamples*Spacing
  m         = np.prod(shape) # total state length
  sample    = simulator(zeros(m), K, 0.0, dt)
  return sample[::Spacing]

sample_filename = 'data/samples/QG_samples.npz'
if not os.path.isfile(sample_filename):
  print('Generating a sample from which to initialize experiments')
  sample = gen_sample(400,500,10)
  np.savez(sample_filename,sample=sample)



