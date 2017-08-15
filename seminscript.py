from addons import *
from tools.viz import plot_benchmark_analysis
from mods.QG.sak08 import setup as setupQG
from mods.Lorenz95.sak08 import setup as setupLorenz

#then run the whole.
if __name__=='__main__':

  config = EnKF('Sqrt',N=50,infl=1.0,rot=True,liveplotting=False)
  p = setupLorenz.h.m

  """

  ################################################################
  ######################## Tunning trust ####################
  ################################################################

  b = Benchmark(setup=setupLorenz,config=config,tunning='Trust',assimcycles=10**4)
  
  b+=MultiDiag(size=p, Rinfl=1.3)

  b.run()
  b.save()

  ################################################################
  ######################## Against MultiDiags ####################
  ################################################################
  
  b = Benchmark(setup=setupLorenz,config=config,tunning='Gain',assimcycles=10**4)
  
  b += MultiDiag(size=p)
  b += MultiDiag(size=p,diags=3,thin=0.9)
  b += MultiDiag(size=p,diags=5,thin=0.85)
  
  b.run()
  b.save()
  del(b)
  
  ################################################################
  ######################## Against diag inflated #################
  ################################################################
  
  b = Benchmark(setup=setupLorenz,config=config,tunning='Gain',assimcycles=10**4)

  b += MultiDiag(size=p,Rinfl=1)
  b += MultiDiag(size=p,Rinfl=1.5)
  b += MultiDiag(size=p,Rinfl=2)
  b += MultiDiag(size=p,Rinfl=3)
  b += MultiDiag(size=p,Rinfl=4)
  b += MultiDiag(size=p,Rinfl=5)

  b.run()
  b.save()
  del(b)
  

  ################################################################
  ######################## Thinning Rt ###########################
  ################################################################
  
  b = Benchmark(Rt=Sparse(size=p),setup=setupLorenz,config=config,tunning='Gain',assimcycles=10**4)
  
  b+=Sparse(size=p,thin=0.8)
  b+=Sparse(size=p,thin=0.6)
  b+=Sparse(size=p,thin=0.4)
  b+=Sparse(size=p,thin=0.3)
  b+=Sparse(size=p,thin=0.2)
  b+=Sparse(size=p,thin=0.1)

  b.run()
  b.save()


  
  ################################################################
  ######################## Thinning SOAR #########################
  ################################################################
  
  b=Benchmark(setup=setupLorenz,config=config,tunning='Gain',assimcycles=10**4)

  b+=SOAR(size=p)
  b+=SOAR(size=p,thin=0.8)
  b+=SOAR(size=p,thin=0.6)
  b+=SOAR(size=p,thin=0.4)
  b+=SOAR(size=p,thin=0.3)
  b+=SOAR(size=p,thin=0.2)
  b+=SOAR(size=p,thin=0.1)

  b.run()
  b.save()
  del(b)
  
  ################################################################
  ######################## Inflate Iden. #########################
  ################################################################

  config = config.update_settings(infl=1.03)
  b=Benchmark(setup=setupLorenz,config=config,tunning=None,assimcycles=10**4)

  b+=MultiDiag(size=p,Rinfl=1)
  b+=MultiDiag(size=p,Rinfl=1.1)
  b+=MultiDiag(size=p,Rinfl=1.2)
  b+=MultiDiag(size=p,Rinfl=1.3)
  b+=MultiDiag(size=p,Rinfl=1.4)
  b+=MultiDiag(size=p,Rinfl=1.5)
  b+=MultiDiag(size=p,Rinfl=1.6)
  b+=MultiDiag(size=p,Rinfl=1.7)
  b+=MultiDiag(size=p,Rinfl=1.8)
  b+=MultiDiag(size=p,Rinfl=1.9)
  b+=MultiDiag(size=p,Rinfl=2.0)
  b+=MultiDiag(size=p,Rinfl=3.0)
  b+=MultiDiag(size=p,Rinfl=4.0)
  b+=MultiDiag(size=p,Rinfl=5.0)
  
  b.run()
  b.save()
  del(b)

  ################################################################
  ######################## Opti = Gain ###########################
  ################################################################

  b=Benchmark(setup=setupLorenz,config=config,tunning='Gain',assimcycles=10**4)

  b+=MultiDiag(size=p,Rinfl=1.0)
  b+=MultiDiag(size=p,Rinfl=1.3)
  b+=MultiDiag(size=p,diags=5)

  b.run()
  b.save()
  del(b)

  
  ################################################################
  ######################## Thinning SOAR (SOAR Rt) ###############
  ################################################################
  
  b=Benchmark(Rt=SOAR(size=p),setup=setupLorenz,config=config,tunning='Trust',assimcycles=10**4)

  b+=SOAR(size=p,thin=0.8)
  b+=SOAR(size=p,thin=0.6)
  b+=SOAR(size=p,thin=0.4)
  b+=SOAR(size=p,thin=0.3)
  b+=SOAR(size=p,thin=0.2)
  b+=SOAR(size=p,thin=0.1)

  b.run()
  b.save()
  del(b)
  
  ################################################################
  ######################## Thinning Sparse (sparse Rt) ###############
  ################################################################

  b = Benchmark(Rt=Sparse(size=p),setup=setupLorenz,config=config,tunning='Gain',assimcycles=10**4)
  
  b+=Sparse(size=p,thin=0.8)
  b+=Sparse(size=p,thin=0.6)
  b+=Sparse(size=p,thin=0.4)
  b+=Sparse(size=p,thin=0.3)
  b+=Sparse(size=p,thin=0.2)
  b+=Sparse(size=p,thin=0.1)

  b.run()
  b.save()

  ################################################################
  ######################## Against MultiDiags ####################
  ################################################################
  
  b = Benchmark(setup=setupLorenz,config=config,tunning='Trust',assimcycles=10**4)
  
  b += MultiDiag(size=p)
  b += MultiDiag(size=p,diags=3,thin=0.9)
  b += MultiDiag(size=p,diags=5,thin=0.85)
  
  b.run()
  b.save()
  del(b)
  
  """
  ################################################################
  ######################## Against diag inflated #################
  ################################################################
  
  b = Benchmark(setup=setupLorenz,config=config,tunning='Gain',assimcycles=10**4)

  b += MultiDiag(size=p,Rinfl=0.5)
  b += MultiDiag(size=p,Rinfl=0.75)
  b += MultiDiag(size=p,Rinfl=1)
  b += MultiDiag(size=p,Rinfl=1.25)
  b += MultiDiag(size=p,Rinfl=1.5)
  b += MultiDiag(size=p,Rinfl=2)

  b.run()
  b.save()
  del(b)
  

  """
  ################################################################
  ######################## Thinning Rt ###########################
  ################################################################
  
  b = Benchmark(Rt=Sparse(size=p),setup=setupLorenz,config=config,tunning='Trust',assimcycles=10**4)
  
  b+=Sparse(size=p,thin=0.8)
  b+=Sparse(size=p,thin=0.6)
  b+=Sparse(size=p,thin=0.4)
  b+=Sparse(size=p,thin=0.3)
  b+=Sparse(size=p,thin=0.2)
  b+=Sparse(size=p,thin=0.1)

  b.run()
  b.save()
  """
