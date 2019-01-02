# Sample (1) M times from N(0,1), and (2) once from N(0,I_m).
# If prior cov B approaches ones((M,M)), i.e. perfect degeneracy,
# then the two should become equivalent.

# Explains why it's min(M,P) that determines weight variation,
# and why serial assimilation of obs is pointless.

############################
# Preamble
############################
from common import *
#sd0 = seed(5)

P = 10
N = 10000
R = 1

xx = randn(N)
yy = randn(P)
y_bar = mean(yy)
ww = sp.stats.norm.pdf(xx,loc=y_bar,scale=sqrt(R/P))
ww /= ww.sum()

#B12 = ones((M,M))
#for i in range(M):
  #B12[i,i] = 0.9
#B = B12 @ B12 / M

############################
# 
############################

plt.hist(ww,bins=100,normed=1)
