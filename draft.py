#draft file for small exec

#Compare effect of inflation factor on different kinds of matrices
"""
d1=data[(data.kind=='MARKOV') | (data.kind=='MARKOV (Rt)')]
f=plt.figure()
ax=f.add_subplot(211)
ax.set_title('MARKOV matrices')
ax.stackplot(d1.infl,d1.RMSE,labels='RMSE',alpha=0.8)
ax.plot(d1.infl,d1.Spread,marker='o',label='Spread',c='black')
ax.set_ylabel('Spread vs RMSE', fontweight='bold')
ax.set_xlabel('Inflation factor',fontweight='bold')
ax.set_ylim(top=1)
ax.legend(loc=2)

d2=data[data.kind=='MultiDiag']
ax2=plt.subplot(212)
ax2.set_title('Identity*1.3 matrices')
ax2.stackplot(d2.infl,d2.RMSE,labels='RMSE',alpha=0.8)
ax2.plot(d2.infl,d2.Spread,marker='o',label='Spread',c='black')
ax2.set_ylabel('Spread vs RMSE',fontweight='bold')
ax2.set_xlabel('Inflation factor',fontweight='bold')
ax2.set_ylim(top=1)
ax2.legend(loc=2)
plt.show()
"""
#Dynamic plot of the QG maps
def show(stats):
  x=stats.xx
  mu=stats.mu.a
  var=stats.var.a

  cmap=plt.cm.viridis
  prep_p = lambda x: square(x)
  prep_q = lambda x: compute_q(square(x))
  prep_m = prep_p
  prep_v = prep_p
          
  f, ((ax_p, ax_q), (ax_m, ax_v)) = plt.subplots(2, 2, figsize=(10,10))
  ax_p.set_title('psi')
  ax_q.set_title('q')
  ax_m.set_title('mean estimate of psi')
  ax_v.set_title('std. dev. in psi')
  supt=plt.suptitle('1/'+str(mu.shape[0]))
  im_p = ax_p.imshow(prep_p(xx[0])        , origin='lower',cmap=cmap)
  im_q = ax_q.imshow(prep_q(xx[0])        , origin='lower',cmap=cmap)
  im_m = ax_m.imshow(prep_m(mu[0])       , origin='lower',cmap=cmap)
  im_v = ax_v.imshow(prep_v(sqrt(var[0])), origin='lower',cmap=cmap)
  ax_p.set_xlim(0,129)
  ax_p.set_ylim(0,129)
  im_p.set_clim(-35,30)
  im_q.set_clim(-1.2e5,1.2e5)
  im_m.set_clim(-35,30)
  im_v.set_clim(0,1.5)
  for i in range(1,mu.shape[0]):
      plt.pause(0.01)
      supt.set_text(str(i+1)+'/'+str(mu.shape[0]))
      im_p.set_data(prep_p(xx[i+1]))
      im_q.set_data(prep_q(xx[i]))
      im_m.set_data(prep_m(mu[i]))
      im_v.set_data(prep_v(sqrt(var[i])))

#Show matrices
from matplotlib import pyplot as plt
def showmat(*args,supt=None,titles=None):
    l = len(args)
    if l%2 == 0 and l//2>1:
      grid = (2,l//2)
      t = t.flatten()
    else:
      grid=(1,l)
    f , t = plt.subplots(grid[0],grid[1],figsize=((grid[1]+1)*3,6))
    if titles == None:
      titles=['']*l
    if not hasattr(t,'__iter__'):
        t=[t]

    for (i,ax) in enumerate(t):
        im = ax.imshow(args[i],cmap=plt.cm.viridis_r,alpha=0.8)
        ax.set_title(titles[i],fontweight='bold')
        ax.axis('off')
    cbaxes = f.add_axes([0.92, 0.1, 0.03, 0.8]) 
    if type(supt) == str:
      f.suptitle(supt,fontweight='bold')
    plt.colorbar(im,cax = cbaxes)
    plt.subplots_adjust(left=0.1, right=0.9, top=0.9, bottom=0.1)
"""
showmat(Sparse(20).full,MARKOV(20).full,SOAR(20).full,MultiDiag(20,diags=3\
  ,decay=2).full,titles=[r'Sparse, $m[i,j]=\frac{1}{1+\frac{d(i,j)}{2}}$' \
  ,r'MARKOV, $m[i,j]=e^{-|i-j|}$',r'SOAR, $m[i,j]=(1+|i-j])e^{-|i-j|}$' \
  ,r'TriDiagonal, $m[i,j]=(|i-j|<3)*2^{|i-j|}$']) \
"""