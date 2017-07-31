from common import *
# Lorenz95.sak08 is fully observated setup: see script -> h=id, forecast noise is zero, 
from mods.Lorenz95.sak08 import setup
#from mods.QG.sak08 import setup
import pickle



"""
The thinning methods could (and should) be built upon the transorm_by method of the CovMat class:
This would imply a big refactorization of this class and this approach: Indeed, the class ObsErrMat is 
redundant as the class CovMat already offers a lot of desired functionalities:
the class obserrmat should only build arrays, used as data for the CovMat class, and the covmat class should store the arguments
used to build the matrix (deltax lr threshold decay whatever they are).Cov mat class should be modified as well so that
<<<<<<< HEAD
the trunc argument could be accessed	not hidden behind self_trunc
=======
the trunc argument could be accessed  not hidden behind self_trunc
>>>>>>> 2bdef7c93fa13f382e83656c206e7918cf54c887

Eg of the thinning method:
#define the mat:
mat=ObsErrMat('MultiDiagonal',diags=5)
C=CovMat(mat)
#The data can then be accessed through C.full (numpy.ndarray)
#modify experiment function so that
C.experiment(xx,yy,setup)
#eventually modify the list_of_matrices class. Benchmark class could remain the same (apparently, eventually add the __repr__ builtin method)

#If thinning is desired:
C.trunc=thinning_value
C.transform_by()
#This would as well require small modification of transform by in order to set the argument function as identity (f: x --> x)
"""


class List_of_matrices(list):

  def __init__(self,*args):
    for mat in args:
      self.append(mat)


  #Probably here add a decorator so that the added matrix is fit to the size requested by the setup
  def __iadd__(self,val):

    if not hasattr(val,'__iter__'):
      val = [val]
    for item in val:
      try:
        self.append(item)
      except RuntimeError:
        print('\nUndefined matrix')
    return self

  def __repr__(self):
    s=''
    for (i,a) in enumerate(self):
      s+=str(i+1)+'.'*4+'\n\t'+a.__str__().replace('\n','\n\t')
      s+='\n-\n'
    return s


class Benchmark(object):
  #Default inflation chosen on paper from sakov 2008

  def __init__(
    self,
    setup=setup,
    config=EnKF('Sqrt', N=50, infl=1.0, rot=True, liveplotting=False),
    Rt=Sparse(size = setup.h.m),
    tunning=None,
    assimcycles=10**4
    ):

    cols=['Setup','R','DA','RMSE','Spread','CPU_Time']

    #set the Rt: this noise will generate the observation
    setup.h.noise.C = Rt
    setup.t.KObs = assimcycles
    self.Rt=Rt
    self._bench=List_of_matrices(Rt)
    self.setup=setup
    self.config=config
    self._output=DataFrame(columns=cols)
    self._is_paused=False
    self.tunning=tunning

    """
    This is so far my best attempt to set default size depending on the current Benchmark in use.
    I of course stopped this investigations to avoid global state.

    Works fine at the first glance however.

    global MARKOV
    MARKOV=functools.partial(MARKOV,size=self.setup.h.m)
    """

  def assess_matrix(self,m,xx,yy):
    #Function aimed at tunning the inflation factor:
    temp=self.config.infl
    #Ugly as possible, should need rework

    r=dict(m.experiment(xx,yy,setup=self.setup,config=self.config))
    if self.tunning is not None:
      word=self.tunning

      def criterion(r,s):
        if word=='Gain':
          return s['RMSE']<=r['RMSE']
        if word=='Trust':
          return s['RMSE']*0.95>s['Spread']

      #t=0
      s=r.copy()
      #Allow a more flexible criterion
      while criterion(r,s):
        #t+=1
        r=s.copy()
        temp+=0.01
        s=dict(m.experiment(xx,yy,setup=self.setup,config=self.config.update_settings(infl=temp)))

      temp=max([temp-0.01,1.0])

    return [('DA',self.config.update_settings(infl=temp))]+list(r.items())


  def run(self):
    #Simulate
    xx,yy=simulate(self.setup)
    #D=diag(diag(self.Rt.matrix)**0.5)

    #Go through the bench
    for (i,m) in enumerate(tqdm.tqdm(self._bench,desc='Benchmark')):
          #self.pause_run()
          if m.rk==m.m:
            r=self.assess_matrix(m,xx,yy)

            #now fill a pandas dataframe with results
            #Create a dictionary of all the useful information for the experiment
            row=dict([('Setup',self.setup.name),('R',m.__class__)]+r)

            if i==0:
              row['kind']+=' (Rt)'
            #Complete the row/dataframe
            for k in set(row.keys())-set(self._output.columns):
              self._output[k]=''
            for l in set(self._output.columns)-set(row.keys()):
              row[l]=''
            self._output.loc[i+1]=row
  
  def save(self):
    
    day=strftime("%m-%d")
    hour=strftime('%H:%M')
    titlecsv='outputcsv-'+hour
    titlepkl='outputpkl-'+hour
    #directory='/Users/remubo/Documents/obs_error_invest/benchmarks/'+day
    directory='data/obs_error_invest/benchmarks/'+day

    if not os.path.exists(directory):
      os.makedirs(directory)

    #Get the arrays column names:
    arrs=[c for c in self.output.columns if type(self.output[c].iloc[0])==numpy.ndarray]
    #Get clean data
    clean=[c for c in self.output.columns if c not in arrs]

    #Build the arrays dictionary:
    dico={name+'-'+str(i):self.output[name].iloc[i] for (name,i) in product(arrs,range(len(self.output)))}

    #Get the rest of the output:
    outcut=self.output[clean]

    #Save the arrays containers
    with open(directory+'/'+titlepkl,'wb') as f:
      pickle.dump(dico,f,pickle.HIGHEST_PROTOCOL)

    #save the csv
    outcut.to_csv(directory+'/'+titlecsv)

    plot_benchmark_analysis(self,skip=True)
    plt.savefig(directory+'/plot-'+hour+'.pdf',bbox_inches='tight')
    plt.close()
    plt.close()
    print('\noutput saved')

  @staticmethod
  def load(day,hour,extra=''):
    #build the paths:
    #directory='/Users/remubo/Documents/obs_error_invest/benchmarks/'+day
    directory='data/obs_error_invest/benchmarks/'+day
    titlecsv='outputcsv-'+hour+extra
    titlepkl='outputpkl-'+hour+extra

    pathpkl=directory+'/'+titlepkl
    pathcsv=directory+'/'+titlecsv

    #Load the csv-saved dataframe:
    df=read_csv(pathcsv,keep_default_na=False)

    #Read the pkled dictionary:
    with open(pathpkl,'rb') as infile:
      dico=pickle.load(infile)

      #Complete the df
      for c in [k.split('-')[0] for k in dico.keys()]:
        df[c]=''
      #Fill in
      for a in dico.items():
        df[a[0].split('-')[0]].iloc[int(a[0].split('-')[1])]=a[1]

    return df

  def pause_run(self):
    key = poll_input()

    if key == '\n' and self._is_paused==False: 
      print('paused')
      self._is_paused = True # pause run
      return None #Avoid the next (unpausing step)
    if self._is_paused:
      # If paused
      ch = getch()
      if ch==' ':
        print('unpaused')
        #remove upper line of tqdm
        CURSOR_UP_ONE = '\x1b[1A'
        ERASE_LINE = '\x1b[2K'
        print(CURSOR_UP_ONE + ERASE_LINE)
        self._is_paused = False

    

  def __repr__(self):
    s=''
    cols=['Setup','Config','Rt','Bench']
    vs=[self.setup.name,self.config,self.Rt,self._bench]
    for (i,j) in zip(cols,vs):
      s+='.'+i+'\n'
      s+='.\t.'+j.__repr__().replace('\n','\n.\t.')+'\n'
    return s

  def __iadd__(self,val):
    self._bench.__iadd__(val)
    return self


  @property
  def bench(self):
    return self._bench

  @property
  def output(self):
    return self._output



#Wrapper for extra setups: