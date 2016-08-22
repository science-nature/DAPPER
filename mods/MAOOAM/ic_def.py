#ic definition
from mods.MAOOAM.params2 import natm,noc,ndim,t_trans,t_run
from mods.MAOOAM.inprod_analytic import awavenum,owavenum,init_inprod
import os.path



def load_IC() :
	if (ndim==0):
		exit('Number of dimensions is 0!')
	if (os.path.exists('ic.py')):
		#no need to create one
		print ("ic already defined")
	else :
		#init_inprod()
		fichier = open("ic.py", "w")
		
		fichier.write("import numpy as np")
		fichier.write("\n")
		fichier.write("\n")
		
		fichier.write("X0=np.zeros("+str(ndim+1)+")")
		fichier.write("\n")
		fichier.write("\n")
		
		fichier.write("X0[0]=1. #one dimension more to have ndim+1 to support the multiplication with the tensor (j and k are in [|0,ndim|]")
		fichier.write("\n")
		fichier.write("\n")
		
		fichier.write("#psi variables")
		fichier.write("\n")
		for i in range(0,natm):
			fichier.write("X0["+str(i+1)+"]=0.0 # typ="+str(awavenum[i].typ)+", NX0="+str(awavenum[i].Nx)+", Ny= "+str(awavenum[i].Ny))
			fichier.write("\n")

		fichier.write("\n")
		fichier.write("#theta variables")
		fichier.write("\n")
		
		for i in range(0,natm):
			fichier.write("X0["+str(i+1+natm)+"]=0.0 # typ="+str(awavenum[i].typ)+", NX0="+str(awavenum[i].Nx)+", Ny= "+str(awavenum[i].Ny))
			fichier.write("\n")
			
		fichier.write("\n")
		fichier.write("#A variables")
		fichier.write("\n")
		
		for i in range(0,noc):
			fichier.write("X0["+str(i+1+natm+natm)+"]=0.0 # NX0="+str(owavenum[i].Nx)+", Ny= "+str(owavenum[i].Ny))
			fichier.write("\n")

		fichier.write("\n")
		fichier.write("#T variables")
		fichier.write("\n")
		
		for i in range(0,noc):
			fichier.write("X0["+str(i+1+natm+natm+noc)+"]=0.0 # NX0="+str(owavenum[i].Nx)+", Ny= "+str(owavenum[i].Ny))
			fichier.write("\n")
		fichier.close()

if __name__ == "__main__":
	load_IC()
