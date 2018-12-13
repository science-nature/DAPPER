# Illustrate how to use DAPPER
# to benchmark a DA method using a "twin experiment".

# Load DAPPER (assumes pwd is <path-to-dapper>)
from common import *

# Load "twin experiment" setup: a hidden Markov Model (HMM)
from mods.Lorenz63.sak12 import HMM
HMM.t.T = 30 # shorten experiment

# Specify a DA method configuration
config = EnKF('Sqrt', N=10, infl=1.02, rot=True, liveplotting=True)

# Simulate synthetic truth (xx) and noisy obs (yy)
xx,yy = simulate(HMM)

# Assimilate yy (knowing the HMM). Assess estimate (vs xx).
stats = config.assimilate(HMM,xx,yy)

# Average stats time series
avrgs = stats.average_in_time()

# Print averages
print_averages(config,avrgs,[],['rmse_a','rmv_a'])

# Plot some diagnostics 
plot_time_series(stats)

# "Explore" objects individually
# print(HMM)
# print(config)
# print(stats)
# print(avrgs)

# Excercise: Try using
# - The (extended) Kalman filter
# - 3D-Var
# - Optimal interpolation
# - The particle filter
# Hint: suggested DA configs are listed in the HMM file.

# Excercise: Repeat the above excercise, but now with the models:
# - Lorenz95
# - LA
# - QG


