# Illustrate how to use DAPPER:
# Basic benchmarking of DA methods.

# Load DAPPER (assumes pwd is <path-to-dapper>)
from common import *

# Load "twin experiment" setup
from mods.Lorenz63.sak12 import setup
setup.t.T = 30 # shorten experiment

# Specify a DA method configuration
config = EnKF('Sqrt', N=10, infl=1.02, rot=True, liveplotting=True)

# Simulate synthetic truth (xx) and noisy obs (yy)
xx,yy = simulate(setup)

# Assimilate yy (knowing the twin setup). Assess estimate (vs xx).
stats = config.assimilate(setup,xx,yy)

# Average stats time series
avrgs = stats.average_in_time()

# Print averages
print_averages(config,avrgs,[],['rmse_a','rmv_a'])

# Plot some diagnostics 
plot_time_series(stats)

# "Explore" objects individually
#print(setup)
#print(config)
#print(stats)
#print(avrgs)

# Excercise1:
# - Try using
#   - The (extended) Kalman filter
#   - 3D-Var
#   - Optimal interpolation
#   - The particle filter
#   Hint: a list of DA configurations can be found in the file
#         from which the setup is imported.

# Excercise2:
# - Repeat the above excercise, now with the model called:
#   - Lorenz95
#   - LA
#   - QG


