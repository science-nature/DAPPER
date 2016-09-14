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

from mods.MAOOAM.maooam16 import setupoc as setup
f,hoc,hoc,chrono,X0 = setup.f, setup.hatm,setup.h, setup.t, setup.X0
kk_f = chrono.kkObsBI-1

DAMs = DAM_list()
# DAMs.add(EnKF_N,N=8,rot=False,liveplotting = False) doesn't work inf or nan
DAMs.add(EnKF_N,N=9,rot=False,liveplotting = False)
DAMs.add(EnKF_N,N=10,rot=False,liveplotting = False)
DAMs.add(EnKF_N,N=11,rot=False,liveplotting = False)
DAMs.add(EnKF_N,N=12,rot=False,liveplotting = False)
DAMs.add(EnKF_N,N=13,rot=False,liveplotting = False)
DAMs.add(EnKF_N,N=14,rot=False,liveplotting = False)
DAMs.add(EnKF_N,N=15,rot=False,liveplotting = False)
DAMs.add(EnKF_N,N=20,rot=False,liveplotting = False)

T=time.clock()

xx=np.loadtxt('./data/truthref.dat')
yy=np.loadtxt('./data/obsref.dat')

yyoc=yy[:,[20,21,22,23,24,25,26,27,28,29,30,31,32,33,34,35]]


print(bcolors.OKBLUE +"truth and observation generated" + bcolors.ENDC)
print (bcolors.OKGREEN +"Time clock :" + bcolors.ENDC ,time.clock()-T)
np.savetxt('./data/obsOC/truth.dat',xx)
np.savetxt('./data/obsOC/obs.dat',yy)


for k,cfg in enumerate(DAMs):
	print(bcolors.OKBLUE +"method" + bcolors.ENDC,cfg.da_method,bcolors.OKBLUE +"size N"+ bcolors.ENDC,cfg.N)
	#assimilation
	s = assimilate(setup,cfg,xx,yyoc,yyoc)
	#free run
	# fr = zeros((chrono.K+1,f.m))
	# fr[0] = s.mu[0]
	# for k,_,t,dt in chrono.forecast_range:
	#   fr[k] = f.model(fr[k-1],t-dt,dt) + sqrt(dt)*f.noise.sample(1)

    #writing on files
	y1=s.matcovf
	y2=s.matcova
	y3=s.errf
	y4=s.erra
	# extr=np.hstack((0,kk_f))
	# y6=(xx-fr)[extr,:]
	la=len(s.Xa[:,1,1])
	np.savetxt('./data/obsOC/'+str(cfg.N)+'/xa1.dat',s.Xa[1,:,:])
	np.savetxt('./data/obsOC/'+str(cfg.N)+'/xa182.dat',s.Xa[floor(la/2),:,:])
	np.savetxt('./data/obsOC/'+str(cfg.N)+'/xa363.dat',s.Xa[la-1,:,:])

	lf=len(s.Xf[:,1,1])
	np.savetxt('./data/obsOC/'+str(cfg.N)+'/xf1.dat',s.Xf[1,:,:])
	np.savetxt('./data/obsOC/'+str(cfg.N)+'/xf182.dat',s.Xf[floor(lf/2),:,:])
	np.savetxt('./data/obsOC/'+str(cfg.N)+'/xf363.dat',s.Xf[lf-1,:,:])

	np.savetxt('./data/obsOC/'+str(cfg.N)+'/obsvar.dat',diag(hoc.noise.C.C))
	# np.savetxt('./data/obsOC/'+str(cfg.N)+'/freerun.dat',fr)
	np.savetxt('./data/obsOC/'+str(cfg.N)+'/ensemblemean.dat',s.mu)
	np.savetxt("./data/obsOC/"+str(cfg.N)+"/matcovf.dat",y1)
	np.savetxt("./data/obsOC/"+str(cfg.N)+"/matcova.dat",y2)
	np.savetxt("./data/obsOC/"+str(cfg.N)+"/errf.dat",y3)
	np.savetxt("./data/obsOC/"+str(cfg.N)+"/erra.dat",y4)
	# np.savetxt("./data/obsOC/"+str(cfg.N)+"/errfr.dat",y6)
	fichier6 = open("./data/obsOC/"+str(cfg.N)+"/info.dat", "w")
	fichier6.write(str(chrono))
	fichier6.write("\n")
	fichier6.write("system's size n :" + str(36))
	fichier6.write("\n")
	fichier6.write("observation's size d :" + str(hoc.m))
	fichier6.write("\n")
	fichier6.write("method :" + str(cfg.da_method) + "  rot" + str(cfg.rot) )
	fichier6.write("\n")
	fichier6.write("ensemble's size N :" + str(cfg.N))
	fichier6.write("\n")
	fichier6.close()
	print (bcolors.OKGREEN +"Time clock :" + bcolors.ENDC ,time.clock()-T)

