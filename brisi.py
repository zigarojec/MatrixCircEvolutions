from reproduction import *
from utils import *

def test():
  vec = createRandomValueVector()
  a = createRandomBigCircuitMatrix(vec)
  makeNetlist(a, 111, 111, a.fullRedundancyMatrix)
  
#In [2]: x = [7, 15, 27, 35, 43, 55, 75]

#In [3]: y = [2.55, 7.8, 18.9, 31.9, 41.3, 66.3, 118]
