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

#weakly
from mods.MAOOAM.maooam16 import setupmix as setup
f,hatm,h,chrono,X0 = setup.f, setup.hatm,setup.h, setup.t, setup.X0
kk_f = chrono.kkObsBI-1


cfg           = DAM(EnKF_N)
cfg.N         = 15
cfg.rot       = False

Tt=time.clock()

for dtObs,cfg.rapportocatm in [[2.2,120],[4.4,60],[8.9,30],[26.7,10],[62.2,4]]:

  #truth
  setup.t= Chronology(0.1,dtObs=dtObs,T=3244.4,BurnIn=0)
  chrono = setup.t
  xxref=np.loadtxt('./data/truthref2.dat')
  xx = zeros((chrono.K+1,f.m))
  xx=xxref[chrono.kk,:]

  yy = zeros((chrono.KObs+1,h.m))
  for k,t in enumerate(chrono.ttObs):
    yy[k] = h.model(xx[chrono.kkObs[k]],t) + h.noise.sample(1)

  yyatm=yy[:,[0,1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19]]

  print(bcolors.OKBLUE +"truth and observation generated" + bcolors.ENDC)
  print (bcolors.OKGREEN +"Time clock :" + bcolors.ENDC ,time.clock()-Tt)
  np.savetxt('./data/strongly/dtobsATM/'+str(dtObs)+'/truth.dat',xx)
  np.savetxt('./data/strongly/dtobsATM/'+str(dtObs)+'/obsatm.dat',yyatm)
  np.savetxt('./data/strongly/dtobsATM/'+str(dtObs)+'/obs.dat',yy)


  print(bcolors.OKBLUE +"method" + bcolors.ENDC,cfg.da_method,bcolors.OKBLUE +"\delta obs"+ bcolors.ENDC,str(dtObs))
  #assimilation
  s = assimilate(setup,cfg,xx,yyatm,yy)
  
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
  np.savetxt('./data/strongly/dtobsATM/'+str(dtObs)+'/xa1.dat',s.Xa[1,:,:])
  np.savetxt('./data/strongly/dtobsATM/'+str(dtObs)+'/xa182.dat',s.Xa[floor(l/2),:,:])
  np.savetxt('./data/strongly/dtobsATM/'+str(dtObs)+'/xa363.dat',s.Xa[l-1,:,:])

  l=len(s.Xf[:,1,1])
  np.savetxt('./data/strongly/dtobsATM/'+str(dtObs)+'/xf1.dat',s.Xf[1,:,:])
  np.savetxt('./data/strongly/dtobsATM/'+str(dtObs)+'/xf182.dat',s.Xf[floor(l/2),:,:])
  np.savetxt('./data/strongly/dtobsATM/'+str(dtObs)+'/xf363.dat',s.Xf[l-1,:,:])

  np.savetxt('./data/strongly/dtobsATM/'+str(dtObs)+'/obsvaratm.dat',diag(hatm.noise.C.C))
  np.savetxt('./data/strongly/dtobsATM/'+str(dtObs)+'/obsvar.dat',diag(h.noise.C.C))
  # np.savetxt('./data/strongly/dtobsATM/'+str(cfg.N)+'/freerun.dat',fr)
  np.savetxt('./data/strongly/dtobsATM/'+str(dtObs)+'/ensemblemean.dat',s.mu)
  np.savetxt("./data/strongly/dtobsATM/"+str(dtObs)+"/matcovf.dat",y1)
  np.savetxt("./data/strongly/dtobsATM/"+str(dtObs)+"/matcova.dat",y2)
  np.savetxt("./data/strongly/dtobsATM/"+str(dtObs)+"/errf.dat",y3)
  np.savetxt("./data/strongly/dtobsATM/"+str(dtObs)+"/erra.dat",y4)
  # np.savetxt("./data/strongly/dtobsATM/"+str(cfg.N)+"/errfr.dat",y6)
  fichier6 = open("./data/strongly/dtobsATM/"+str(dtObs)+"/info.dat", "w")
  fichier6.write(str(chrono))
  fichier6.write("\n")
  fichier6.write("system's size n :" + str(36))
  fichier6.write("\n")
  fichier6.write("observation's size d :" + str(hatm.m) + str(h.m))
  fichier6.write("\n")
  fichier6.write("method :" + str(cfg.da_method) + "  rot" + str(cfg.rot) )
  fichier6.write("\n")
  fichier6.write("ensemble's size N :" + str(cfg.N))
  fichier6.write("\n")
  fichier6.close()
  print (bcolors.OKGREEN +"Time clock :" + bcolors.ENDC ,time.clock()-Tt)

