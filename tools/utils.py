# Utilities (non-math)

from common import *

#########################################
# Progressbar
#########################################
def noobar(itrble, desc):
  """Simple progress bar. To be used if tqdm not installed."""
  L  = len(itrble)
  print('{}: {: >2d}'.format(desc,0), end='')
  for k,i in enumerate(itrble):
    yield i
    p = (k+1)/L
    e = '' if k<(L-1) else '\n'
    print('\b\b\b\b {: >2d}%'.format(int(100*p)), end=e)
    sys.stdout.flush()

# Get progbar description by inspecting caller function.
import inspect
def pdesc(desc):
  if desc is not None:
    return desc
  try:
    #stackoverflow.com/q/15608987
    DAC_name  = inspect.stack()[3].frame.f_locals['name_hook']
  except (KeyError, AttributeError):
    #stackoverflow.com/a/900404
    DAC_name  = inspect.stack()[2].function
  return DAC_name 

# Define progbar as tqdm or noobar
try:
  import tqdm
  def progbar(inds, desc=None, leave=1):
    if is_notebook:
      pb = tqdm.tqdm_notebook(inds,desc=pdesc(desc),leave=leave)
    else:
      pb = tqdm.tqdm(inds,desc=pdesc(desc),leave=leave,smoothing=0.3,dynamic_ncols=True)
    # Printing during the progbar loop (may occur with error printing)
    # can cause tqdm to freeze the entire execution. 
    # Seemingly, this is caused by their multiprocessing-safe stuff.
    # Disable this, as per github.com/tqdm/tqdm/issues/461#issuecomment-334343230
    # try: pb.get_lock().locks = []
    # except AttributeError: pass
    return pb
except ImportError as err:
  install_warn(err)
  def progbar(inds, desc=None, leave=1):
    return noobar(inds,desc=pdesc(desc))




#########################################
# Console input / output
#########################################

#stackoverflow.com/q/292095
import select
def poll_input():
  i,o,e = select.select([sys.stdin],[],[],0.0001)
  for s in i: # Only happens if <Enter> has been pressed
    if s == sys.stdin:
      return sys.stdin.readline()
  return None

# Can't get thread solution working (in combination with getch()):
# stackoverflow.com/a/25442391/38281

# stackoverflow.com/a/21659588/38281
# (Wait for) any key:
def _find_getch():
    try:
        import termios
    except ImportError:
        # Non-POSIX. Return msvcrt's (Windows') getch.
        import msvcrt
        return msvcrt.getch
    # POSIX system. Create and return a getch that manipulates the tty.
    import sys, tty
    def _getch():
        fd = sys.stdin.fileno()
        old_settings = termios.tcgetattr(fd)
        try:
            tty.setraw(fd)
            ch = sys.stdin.read(1)
        finally:
            termios.tcsetattr(fd, termios.TCSADRAIN, old_settings)
        return ch
    return _getch
getch = _find_getch()



import inspect
def spell_out(*args):
  """
  Print (args) including variable names.
  Example:
  >>> print(3*2)
  >>> 3*2:
  >>> 6
  """
  frame  = inspect.stack()[1].frame
  lineno = frame.f_lineno

  # This does not work coz they (both) end reading when encountering a func def
  # code  = inspect.getsourcelines(frame)
  # code  = inspect.getsource(frame)

  # Instead, read the source manually...
  f = inspect.getfile(frame)
  if "ipython-input" in f:
    # ( does not work with pure python where f == "<stdin>" )
    line = inspect.getsource(frame)
  else:
    try:
      f = os.path.relpath(f)
      with open(f) as stream:
        line = stream.readlines()[lineno-1]
    except FileNotFoundError:
      line = "(Print command on line number " + str(lineno) + " [of unknown source])"

  # Find opening/closing brackets
  left  = line. find("(")
  right = line.rfind(")")

  # Print header
  with coloring(cFG.MAGENTA):
    print(line[left+1:right] + ":")
  # Print (normal)
  print(*args)


def print_together(*args):
  "Print 1D arrays stacked together."
  print(np.vstack(args).T)


# Local np.set_printoptions. stackoverflow.com/a/2891805/38281
import contextlib
@contextlib.contextmanager
@functools.wraps(np.set_printoptions)
def printoptions(*args, **kwargs):
    original = np.get_printoptions()
    np.set_printoptions(*args, **kwargs)
    yield 
    np.set_printoptions(**original)


import tools.tabulate as tabulate_orig
tabulate_orig.MIN_PADDING = 0
def tabulate(data,headr=(),formatters=(),inds='nice'):
  """
  Pre-processor for tabulate().
  Main task: transpose 'data' (list-of-lists).
  If 'data' is a dict, the 'headr' will be keys.
  'formatter': define formats to apply before relaying to pandas.
               Default: attr.__name__ (when applicable).
  Example:
  >>> print(tabulate(cfgs.distinct_attrs()))
  """

  # Extract dict
  if hasattr(data,'keys'):
    headr = list(data)
    data  = data.values()

  # Default formats
  if not formatters:
    formatters = ({
        'test'  : lambda x: hasattr(x,'__name__'),
        'format': lambda x: x.__name__
        },)
  # Apply formatting (if not applicable, data is just forwarded)
  for f in formatters:
    data = [[f['format'](j) for j in row] if f['test'](row[0]) else row for row in data]

  # Transpose
  data = list(map(list, zip(*data)))

  # Generate nice indices
  if inds=='nice':
    inds = ['[{}]'.format(d) for d in range(len(data))]
  else:
    pass # Should be True or False

  return tabulate_orig.tabulate(data,headr,showindex=inds)


def repr_type_and_name(thing):
  """Print thing's type [and name]"""
  s = "<" + type(thing).__name__ + '>'
  if hasattr(thing,'name'):
    s += ': ' + thing.name
  return s


class MLR_Print:
  """
  Multi-Line, Recursive repr (print) functionality.
  Set class variables to change look:
   - 'indent': indentation per level
   - 'ch': character to use for "spine" (e.g. '|' or ' ')
   - 'ordr_by_linenum': 0: alphabetically, 1: linenumbr, -1: reverse
  """
  indent=3
  ch='.'
  ordr_by_linenum = 0

  # numpy print options
  threshold=10
  precision=None

  # Recursion monitoring.
  _stack=[] # Reference using MLR_Print._stack, ...
  # not self._stack or self.__class__, which reference sub-class "instance".

  # Reference using self.excluded, to access sub-class "instance".
  excluded = [] # Don't include in printing
  excluded.append(re.compile('^_')) # "Private"
  excluded.append('name') # Treated separately

  included = []
  aliases  = {}

  def __repr__(self):
    with printoptions(threshold=self.threshold,precision=self.precision):
      # new line chars
      NL = '\n' + self.ch + ' '*(self.indent-1)

      # Infinite recursion prevention
      is_top_level = False
      if MLR_Print._stack == []:
        is_top_level = True
      if self in MLR_Print._stack:
        return "**Recursion**"
      MLR_Print._stack += [self]

      # Use included or filter-out excluded
      keys = self.included or filter_out(vars(self), *self.excluded)

      # Process attribute repr's
      txts = {}
      for key in keys:
        t = repr(getattr(self,key)) # sub-repr
        if '\n' in t:
          # Activate multi-line printing
          t = t.replace('\n',NL+' '*self.indent) # other lines
          t = NL+' '*self.indent + t             # first line
        t = NL + self.aliases.get(key,key) + ': ' + t # key-name
        txts[key] = t

      def sortr(x):
        if self.ordr_by_linenum:
          return self.ordr_by_linenum*txts[x].count('\n')
        else:
          return x.lower()

      # Assemble string
      s = repr_type_and_name(self)
      for key in sorted(txts, key=sortr):
        s += txts[key]

      # Empty _stack when top-level printing finished
      if is_top_level:
        MLR_Print._stack = []

      return s


class AlignedDict(OrderedDict):
  """Provide aligned-printing for dict."""
  def __init__(self,*args,**kwargs):
    super().__init__(*args,**kwargs)
  def __str__(self):
    L = max([len(s) for s in self.keys()], default=0)
    s = " "
    for key in self.keys():
      s += key.rjust(L)+": "+repr(self[key])+"\n "
    s = s[:-2]
    return s
  def __repr__(self):
    return type(self).__name__ + "(**{\n " + str(self).replace("\n",",\n") + "\n})"
  def _repr_pretty_(self, p, cycle):
    # Is implemented by OrderedDict, so must overwrite.
    if cycle: p.text('{...}')
    else:     p.text(self.__repr__())
  
class Bunch(dict):
  def __init__(self,**kw):
    dict.__init__(self,kw)
    self.__dict__ = self


# From stackoverflow.com/q/22797580 and more
class NamedFunc():
  "Provides custom repr for functions."
  def __init__(self,_func,_repr):
    self._func = _func
    self._repr = _repr
    #functools.update_wrapper(self, _func)
  def __call__(self, *args, **kw):
    return self._func(*args, **kw)
  def __repr__(self):
    argnames = self._func.__code__.co_varnames[
      :self._func.__code__.co_argcount]
    argnames = "("+",".join(argnames)+")"
    return "<NamedFunc>"+argnames+": "+self._repr

class NameFunc():
  "Decorator version"
  def __init__(self,name):
     self.fn_name = name
  def __call__(self,fn):
      return NamedFunc(fn,self.fn_name)



#########################################
# Misc
#########################################

# # Temporarily set attribute values
# @contextlib.contextmanager
# def set_tmp(obj,attr,val):
#     tmp = getattr(obj,attr)
#     setattr(obj,attr,val)
#     yield 
#     setattr(obj,attr,tmp)

import contextlib
@contextlib.contextmanager
def set_tmp(obj, attr, val):
    """
    Temporarily set an attribute.
    code.activestate.com/recipes/577089
    """
    was_there = False
    tmp = None
    if hasattr(obj, attr):
      try:
        if attr in obj.__dict__:
          was_there = True
      except AttributeError:
        if attr in obj.__slots__:
          was_there = True
      if was_there:
        tmp = getattr(obj, attr)
    setattr(obj, attr, val)
    yield #was_there, tmp
    if not was_there: delattr(obj, attr)
    else:             setattr(obj, attr, tmp)



# Better than tic-toc !
import time
class Timer():
  """
  Usage:
  with Timer('<description>'):
    do_stuff()
  """
  def __init__(self, name=None):
      self.name = name

  def __enter__(self):
      self.tstart = time.time()

  def __exit__(self, type, value, traceback):
      #pass # Turn off timer messages
      if self.name:
          print('[%s]' % self.name, end='')
      print('Elapsed: %s' % (time.time() - self.tstart))

def find_1st(xx):
  try:                  return next(x for x in xx if x)
  except StopIteration: return None
def find_1st_ind(xx):
  try:                  return next(k for k in range(len(xx)) if xx[k])
  except StopIteration: return None

# stackoverflow.com/a/2669120
def sorted_human( lst ): 
    """ Sort the given iterable in the way that humans expect.""" 
    convert = lambda text: int(text) if text.isdigit() else text 
    alphanum_key = lambda key: [ convert(c) for c in re.split('([0-9]+)', key) ] 
    return sorted(lst, key = alphanum_key)

def keep_order_unique(arr):
  "Undo the sorting that np.unique() does."
  _, inds = np.unique(arr,return_index=True)
  return arr[np.sort(inds)]

def de_abbreviate(abbrev_d, abbreviations):
  "Expand any abbreviations (list of 2-tuples) that occur in abbrev_d (dict)."
  for a,b in abbreviations:
    if a in abbrev_d:
      abbrev_d[b] = abbrev_d[a]
      del abbrev_d[a]



def filter_out(orig_list,*unwanted,INV=False):
  """
  Returns new list from orig_list with unwanted removed.
  Also supports re.search() by inputting a re.compile('patrn') among unwanted.
  """
  new = []
  for word in orig_list:
    for x in unwanted:
      try:
        # Regex compare
        rm = x.search(word)
      except AttributeError:
        # String compare
        rm = x==word
      if (not INV)==bool(rm):
        break
    else:
      new.append(word)
  return new

def all_but_1_is_None(*args):
  "Check if only 1 of the items in list are Truthy"
  return sum(x is not None for x in args) == 1

# From stackoverflow.com/q/3012421
class lazy_property(object):
    '''
    Lazy evaluation of property.
    Should represent non-mutable data,
    as it replaces itself.
    '''
    def __init__(self,fget):
      self.fget = fget
      self.func_name = fget.__name__

    def __get__(self,obj,cls):
      value = self.fget(obj)
      setattr(obj,self.func_name,value)
      return value



class AssimFailedError(RuntimeError):
    pass

def raise_AFE(msg,time_index=None):
  if time_index is not None:
    msg += "\n(k,kObs,fau) = " + str(time_index) + ". "
  raise AssimFailedError(msg)


def vectorize0(f):
  """
  Vectorize f for its 1st (index 0) argument.

  Compared to np.vectorize:
    - less powerful, but safer to only vectorize 1st argument
    - doesn't evaluate the 1st item twice
    - doesn't always return array

  Example:
  >>> @vectorize0
  >>> def add(x,y):
  >>>   return x+y
  >>> add(1,100)
  101
  >>> x = np.arange(6).reshape((3,-1))
  >>> add(x,100)
  array([[100, 101],
       [102, 103],
       [104, 105]])
  >>> add([20,x],100)
  [120, array([[100, 101],
        [102, 103],
        [104, 105]])]
  """
  @functools.wraps(f)
  def wrapped(x,*args,**kwargs):
    if hasattr(x,'__iter__'):
      out = [wrapped(xi,*args,**kwargs) for xi in x]
      if isinstance(x,np.ndarray):
        out = np.asarray(out)
    else:
      out = f(x,*args,**kwargs)
    return out
  return wrapped


