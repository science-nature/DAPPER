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

# Excercises:
# - Try using the Lorenz95 model.
#   - Find suitable DA configuration (see suggestions in setup file).
# - Repeat, now with the QG model.


