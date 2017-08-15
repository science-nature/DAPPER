#QG test
from common import *
from addons import *
from mods.QG.sak08 import setup,jj,params
from mods.QG.core import *
from draft import showmat

#Francois' config
config=LETKF(loc_rad=15, N=25, taper='Gauss', infl=1.07, rot=True)

obss=jj.copy()
p = setup.h.m
eucl_dist = lambda x,y:sqrt((x%129-y%129)**2+(x//129-y//129)**2)
exp_dec = lambda x,y:exp(-eucl_dist(obss[x],obss[y])/5)
#This a standard distance-decreasing corrmat
R1 = Custom(size=p,f = lambda i,j: 4*exp_dec(i,j))
#now if I get rid of the upper-right and down-left corners, I obtain the correlation where +45 and -45 tracks are not correlated
R2 = Custom(size=p,f = lambda i,j: 4*exp_dec(i,j)*((i<=p//2-1)*(j<=p//2-1) + (i>p//2-1)*(j>p//2-1)) )
#Now a standard diag matrix
R3 = MultiDiag(size=p,Rinfl=4)

#Run
setup.h.noise.C = R2
setup.t.KObs = 200
xx,yy = simulate(setup,reorder=False)
setup.h.noise.C = R1
statsR1 = config.assimilate(setup,xx,yy)

"""
b = Benchmark(Rt=R2,config=config,setup=setup,assimcycles=200)

b += R1
b += R3

b.run()
b.save()
plot_benchmark_analysis(b)
"""