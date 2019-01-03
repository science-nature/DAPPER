# Sample (1) Nx times from N(0,1), and (2) once from N(0,Id_Nx).
# If prior cov B approaches ones((Nx,Nx)), i.e. perfect degeneracy,
# then the two should become equivalent.

# Explains why it's min(Nx,Ny) that determines weight variation,
# and why serial assimilation of obs is pointless.

############################
# Preamble
############################
from common import *
#sd0 = seed(5)

Ny = 10
N  = 10000
R  = 1

xx = randn(N)
yy = randn(Ny)
y_bar = mean(yy)
ww = sp.stats.norm.pdf(xx,loc=y_bar,scale=sqrt(R/Ny))
ww /= ww.sum()

#B12 = ones((Nx,Nx))
#for i in range(Nx):
  #B12[i,i] = 0.9
#B = B12 @ B12 / Nx

############################
# 
############################

plt.hist(ww,bins=100,normed=1)
