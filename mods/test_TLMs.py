# Numerical validation of TLM (derivatives).

##
from common import *

eps = 1e-6
def TLM_approx(dxdt,x):
  "Finite-diff TLM"
  E = x + eps*eye(len(x))
  return ( dxdt(E) - dxdt(x) ).T   @   inv(E-x).T

# Transpose explanation:
# Let A be matrix whose cols are realizations. Abbrev: f=dxdt.
# Assume: f(A)-f(x) ≈ F @ (A-x).
# Then  : F ≈ [f(A)-f(x)] @ inv(A-x)         (v1)
#           = [f(A')-f(x')]' @ inv(A'-x')'.  (v2)
# Since DAPPER uses dxdt that is made for A', it's easier to apply v2.
# However, TLM should compute F (not F').

def compare(TLM, dxdt, x):
  F1 = TLM(x)
  F2 = TLM_approx(dxdt,x)
  # rtol=0 => only atol matters.
  return np.allclose(F1, F2, atol=10*eps, rtol=0)


##
from mods.LotkaVolterra.core import dxdt, TLM, Nx
x = 0.5 + 0.1*randn(Nx)
def test_LV(TLM=TLM,dxdt=dxdt,x=x): # capture current values
  assert compare(TLM, dxdt, x)

##
from mods.Lorenz63.core import dxdt, TLM, Nx
x = 10*randn(Nx)
def test_L63(TLM=TLM,dxdt=dxdt,x=x): # capture current values
  assert compare(TLM, dxdt, x)

##
from mods.Lorenz84.core import dxdt, TLM, Nx
x = 10*randn(Nx) # ?
def test_L84(TLM=TLM,dxdt=dxdt,x=x): # capture current values
  assert compare(TLM, dxdt, x)

##
from mods.Lorenz95.core import dxdt, TLM
x = 5 + randn(10)
def test_L95(TLM=TLM,dxdt=dxdt,x=x): # capture current values
  assert compare(TLM, dxdt, x)

##
# TODO
# from mods.LorenzUV.core import model_instance
# LUV = model_instance(nU=10,J=4,F=10)
# x = 5 + randn(LUV.M)
# LUV.dfdt
# def test_LUV(TLM=,dxdt=,x=): # capture current values
  # assert compare(TLM, dxdt, x)

##


