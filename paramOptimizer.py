"""
paramOptimizer.py
Ziga Rojec, EDA group, Faculty of Electrical Engineering, University of Ljubljana, Slovenia. Jan 2018.

This script is used by the additional method, that optimizes numerical parameters of a single topology. This is used when optimize is set to 1 in globalVars.py.

"""

from utils import circuit, slimCircuit
from copy import copy, deepcopy
import random
from utils import returnMaxValueVector, returnMinValueVector

from scoreFunctions import *

#Optimisation modules
from pyopus.optimizer.psade import ParallelSADE
from pyopus.optimizer.hj import HookeJeeves
from pyopus.optimizer.qpmads import QPMADS
from pyopus.optimizer.boxcomplex import BoxComplex
from pyopus.optimizer.base import Reporter, CostCollector, RandomDelay

#Available problems - score functions
problems = {'scoreCirc_CmosVoltageReference_2':scoreCirc_CmosVoltageReference_2,
	    'scoreCirc_ActiveFilter_2':scoreCirc_ActiveFilter_2,
	    'scoreCirc_PassiveBandPass':scoreCirc_PassiveBandPass,
	    'scoreCirc_HighPass':scoreCirc_HighPass,
	    } #set this also in main
PROBLEM = problems[PROBLEMname]

  
#naredi class, ki bo ob inicializaciji vzel topologijo, klic funkcije pa bo vzel ValueVector KUUUL!

class circuitUnderOptimiser:
  """Used for single individual (circuit) parameter optimization.
     Initalizing this class, a topology matrix is needed.
     Calling an object from this class, a ValueVector is needed. 
     When in call, a circuit object is generated.
     We call a cost/score function to evauate an individual in MOEAMODE == False 
  """
  def __init__(self, BigCircuitMatrix):
    self.BigCircuitMatrix = BigCircuitMatrix
    self.fullMatrix = fullRedundancyBigCircuitMatrix(self.BigCircuitMatrix)
    #self.fullMatrix = fullRedundancyBigCircuitMatrix(BigCircuitMatrix)
  def __call__(self, ValueVector):
    """Problem to optimise."""
    circuitHOT = slimCircuit(self.BigCircuitMatrix, ValueVector)
    circuitHOT.fullRedundancyMatrix = self.fullMatrix
    score, results = PROBLEM(circuitHOT, 1, random.randint(1,100000),False) #HERE! Every evaluated file gets a random number between 1 and 100000, since they spawn in the root folder when run in single-thread mode. With a dirty and simple solution comes a stupid hazard, since files two can get the same name with the probability of 1/1000000. You have been warned. 
    return score

probLOW = returnMinValueVector()
probHIGH = returnMaxValueVector()

def optimiseCircuit(topology, values, maxIter):
  """Takes topology and optimizes values using optimisation method.
      Basically a wrapper for circuitUnderOptimiser class.
      See circuitUnderOptimiser class.
      Here, the parameter optimization algorithm is run. 
  """
  print "... optimising values ..."
  problem = deepcopy(circuitUnderOptimiser(topology))
  initValues = copy(values)
  
  opt = ParallelSADE(problem, probLOW, probHIGH, debug=0, maxiter=maxIter)
  #opt = HookeJeeves(problem, probLOW, probHIGH, debug=10, maxiter=10000)
  #opt=BoxComplex(problem, probLOW, probHIGH, maxiter=10000, fstop = 0.1) #settings
  
  cc=CostCollector()
  opt.installPlugin(cc)
  #opt.installPlugin(Reporter(onIterStep=10000))
  #opt.check()
  opt.reset(initValues)
  opt.run()
  cc.finalize()
  print("x=%s f=%e" % (str(opt.x), opt.f))
  return opt.x, opt.f

