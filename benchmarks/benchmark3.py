##################################################################################
#
#STRONGLY
#changing dtObs^atm
#fully observed
###################################################################################
from common import *
import time
np.random.seed(5)

class bcolors:
#to color the instructions in the console
  HEADER = '\033[95m'
  OKBLUE = '\033[94m'
  OKGREEN = '\033[92m'
  WARNING = '\033[93m'
  FAIL = '\033[91m'
  ENDC = '\033[0m'
  BOLD = '\033[1m'
  UNDERLINE = '\033[4m'

from mods.MAOOAM.maooam16 import setup
f,h,h,chrono,X0 = setup.f, setup.h, setup.h, setup.t, setup.X0
kk_f = chrono.kkObsBI-1


cfg           = DAM(EnKF_N)
cfg.N         = 15
cfg.rot       = False

Tt=time.clock()

for dtObs in [0.4,4.4,8.9,26.7,62.2]:

  #truth
  setup.t= Chronology(0.1,dtObs=dtObs,T=3244.4,BurnIn=0)
  f,h,h,chrono,X0 = setup.f, setup.h, setup.h, setup.t, setup.X0
  xxref=np.loadtxt('./data/truthref2.dat')
  xx = zeros((chrono.K+1,f.m))
  xx=xxref[chrono.kk,:]
  # obs
  yy = zeros((chrono.KObs+1,h.m))
  for k,t in enumerate(chrono.ttObs):
    yy[k] = h.model(xx[chrono.kkObs[k]],t) + h.noise.sample(1)

  print(bcolors.OKBLUE +"truth and observation generated" + bcolors.ENDC)
  print (bcolors.OKGREEN +"Time clock :" + bcolors.ENDC ,time.clock()-Tt)
  np.savetxt('./data/strongly/dtObs/'+str(dtObs)+'/truth.dat',xx)
  np.savetxt('./data/strongly/dtObs/'+str(dtObs)+'/obs.dat',yy)

  print(bcolors.OKBLUE +"method" + bcolors.ENDC,cfg.da_method,bcolors.OKBLUE +"\delta obs"+ bcolors.ENDC,str(dtObs))
  #assimilation
  s = assimilate(setup,cfg,xx,yy,yy)
  
  # if dtObs==4.4 :
  #   #free run because N is constant
  #   fr = zeros((chrono.K+1,f.m))
  #   fr[0] = s.mu[0]
  #   for k,_,t,dt in chrono.forecast_range:
  #     fr[k] = f.model(fr[k-1],t-dt,dt) + sqrt(dt)*f.noise.sample(1)
  #     extr=np.hstack((0,kk_f))
  #     y6=(xx-fr)[extr,:]
  #     np.savetxt("./data/errfr.dat",y6)
  #     np.savetxt('./data/freerun.dat',fr)
    #writing on files
  y1=s.matcovf
  y2=s.matcova
  y3=s.errf
  y4=s.erra
  extr=np.hstack((0,kk_f))
  # y6=(xx-fr)[extr,:]
  l=len(s.Xa[:,1,1])
  np.savetxt('./data/strongly/dtObs/'+str(dtObs)+'/xa1.dat',s.Xa[1,:,:])
  np.savetxt('./data/strongly/dtObs/'+str(dtObs)+'/xa182.dat',s.Xa[floor(l/2),:,:])
  np.savetxt('./data/strongly/dtObs/'+str(dtObs)+'/xa363.dat',s.Xa[l-1,:,:])

  l=len(s.Xf[:,1,1])
  np.savetxt('./data/strongly/dtObs/'+str(dtObs)+'/xf1.dat',s.Xf[1,:,:])
  np.savetxt('./data/strongly/dtObs/'+str(dtObs)+'/xf182.dat',s.Xf[floor(l/2),:,:])
  np.savetxt('./data/strongly/dtObs/'+str(dtObs)+'/xf363.dat',s.Xf[l-1,:,:])

  np.savetxt('./data/strongly/dtObs/'+str(dtObs)+'/obsvar.dat',diag(h.noise.C.C))
  # np.savetxt('./data/strongly/dtObs/'+str(cfg.N)+'/freerun.dat',fr)
  np.savetxt('./data/strongly/dtObs/'+str(dtObs)+'/ensemblemean.dat',s.mu)
  np.savetxt("./data/strongly/dtObs/"+str(dtObs)+"/matcovf.dat",y1)
  np.savetxt("./data/strongly/dtObs/"+str(dtObs)+"/matcova.dat",y2)
  np.savetxt("./data/strongly/dtObs/"+str(dtObs)+"/errf.dat",y3)
  np.savetxt("./data/strongly/dtObs/"+str(dtObs)+"/erra.dat",y4)
  # np.savetxt("./data/strongly/dtObs/"+str(cfg.N)+"/errfr.dat",y6)
  fichier6 = open("./data/strongly/dtObs/"+str(dtObs)+"/info.dat", "w")
  fichier6.write(str(chrono))
  fichier6.write("\n")
  fichier6.write("system's size n :" + str(36))
  fichier6.write("\n")
  fichier6.write("observation's size d :" + str(h.m))
  fichier6.write("\n")
  fichier6.write("method :" + str(cfg.da_method) + "  rot" + str(cfg.rot) )
  fichier6.write("\n")
  fichier6.write("ensemble's size N :" + str(cfg.N))
  fichier6.write("\n")
  fichier6.close()
  print (bcolors.OKGREEN +"Time clock :" + bcolors.ENDC ,time.clock()-Tt)

