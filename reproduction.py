import numpy as np
import random
from copy import copy, deepcopy
import os
from time import time
import pyopus.evaluator.measure as pyo

from globalVars import *
import AdamAndEve as AE
#import runme2
from utils import *
import pickle
import sys
import collections

##Optimisation modules *moved to paramOptimizer.py along with circuitUnderOptimiser class and functions
#from pyopus.optimizer.psade import ParallelSADE
#from pyopus.optimizer.hj import HookeJeeves
#from pyopus.optimizer.qpmads import QPMADS
#from pyopus.optimizer.boxcomplex import BoxComplex
#from pyopus.optimizer.base import Reporter, CostCollector, RandomDelay


#class population:
#	"""TODO."""
#	def __init__(self):
#		self.pool = np.array([])
#		self.NofIndividuals = len(self.pool)
#		self.scores = np.array([])
#		self.matrixDensities = np.array([])
#		self.matrixQuaziIDs = np.array([])
#		self.matrixHashes = np.array([])
#
#	def add_individual(self, individual):
#		self.pool = np.append(self.pool, individual)
#
#	def add_score(self, score):
#		self.scores = np.append(self.scores, score)
#
#	def add_matrixDensity(self, matrixDensity):
#		self.matrixDensities = np.append(self.matrixDensities, matrixDensity)
#
#	def add_matrixQuaziID(self, matrixQuaziID):
#		self.matrixQuaziIDs = np.append(self.matrixQuaziIDs, matrixQuaziID)
#		
#	def add_matrixHash(self, matrixHash):
#		self.matrixHashes = np.append(self.matrixHashes, matrixHash)
		
class population:
	"""TODO."""
	def __init__(self):
		self.pool = np.array([])
		self.NofIndividuals = len(self.pool)
		self.scores = np.array([])
		self.matrixDensities = np.array([])
		self.matrixQuaziIDs = np.array([])
		self.matrixHashes = np.array([])
		self.fronts = []

	def add_individual(self, individual):
		self.pool = np.append(self.pool, individual)

	def add_score(self, score):
		self.scores = np.append(self.scores, score)

	def add_matrixDensity(self, matrixDensity):
		self.matrixDensities = np.append(self.matrixDensities, matrixDensity)

	def add_matrixQuaziID(self, matrixQuaziID):
		self.matrixQuaziIDs = np.append(self.matrixQuaziIDs, matrixQuaziID)
		
	def add_matrixHash(self, matrixHash):
		self.matrixHashes = np.append(self.matrixHashes, matrixHash)

	def __iter__(self):
	    """Allows for iterating over Individuals (Circuits)"""
	    return self.pool.__iter__()		


#Simple continious evolution algorithm (Value changing)
def valueReproduction(ind1, ind2, **crossover):
  """"""
  ind1 = copy(ind1)
  ind2 = copy(ind2)
  
  #print ind1.ValueVector
  #print ind2.ValueVector
  #raw_input()
  
  valueV1 = copy(ind1.ValueVector)
  valueV2 = copy(ind2.ValueVector)  
  
  if crossover['crossover'] == "lineCrossover":
    #Exchange ONE chromoshome between the individuals
    slicePoint = random.randint(0,len(valueV1)-1)
    temp = copy(valueV1[slicePoint])
    valueV1[slicePoint] = copy(valueV2[slicePoint])
    valueV2[slicePoint] = copy(temp)
    
  elif crossover['crossover'] == "sliceCrossover":
    #Exchange several chromosomes between individuals
    slicePoint = random.randint(0,len(valueV1)-1)
    temp = copy(valueV1)
    valueV1 = copy(np.append(valueV2[:slicePoint], temp[slicePoint:]))
    valueV2 = copy(np.append(temp[:slicePoint], valueV2[slicePoint:]))
    
  elif crossover['crossover'] == "rangeCrossover":
    #Exchange ONE chromoshome between the individuals, each chromosome gets new value between the two
    slicePoint = random.randint(0,len(valueV1)-1)
    temp1 = copy(valueV1[slicePoint])
    temp2 = copy(valueV2[slicePoint])
    valueV1[slicePoint] = copy(random.uniform(temp1, temp2))
    valueV2[slicePoint] = copy(abs(valueV1[slicePoint]-temp2))

  elif crossover['crossover'] == "mutation" :
    #Mutate one chromoshome in a ValueVector for both parents. Return mutated parents. 
    slicePoint = random.randint(0,len(valueV1)-1)
    valueV1[slicePoint] = createRandomValueVector()[slicePoint]
    valueV2[slicePoint] = createRandomValueVector()[slicePoint]

  #print valueV1
  #print valueV2
  #raw_input()

  return circuit(ind1.BigCircuitMatrix, valueV1), circuit(ind2.BigCircuitMatrix, valueV2)
  
#Topology evolution tools

def rowwiseCrossover(ind1, ind2):
	"""Takes exactly two individuals of class "circuit" and mates them with simple crossover. 
		More soon."""
	ind1_MX = np.array(ind1.BigCircuitMatrix)
	ind2_MX = np.array(ind2.BigCircuitMatrix)
	
	OutConns_gene_1 = np.zeros([NofOutConns,BigMatrixSize-NofOutConns], dtype=bool)
	OutConns_gene_2 = np.zeros([NofOutConns,BigMatrixSize-NofOutConns], dtype=bool)
	InterConns_gene_1 = np.zeros([BigMatrixSize-NofOutConns, BigMatrixSize-NofOutConns], dtype=bool)
	InterConns_gene_2 = np.zeros([BigMatrixSize-NofOutConns, BigMatrixSize-NofOutConns], dtype=bool)
	
	
	#Extract sexual material - OutConns - Individual_1
	for i in range(len(OutConns_gene_1)):
		OutConns_gene_1[i] = deepcopy(ind1_MX[0:BigMatrixSize-NofOutConns,BigMatrixSize-(i+1)])
	
	#Extract sexual material - OutConns - Individual_2
	for i in range(len(OutConns_gene_2)):
		OutConns_gene_2[i] = deepcopy(ind2_MX[0:BigMatrixSize-NofOutConns,BigMatrixSize-(i+1)])
		
	#Extract sexual material - InterConns - Individual_1	
	for i in range(len(InterConns_gene_1)):
		InterConns_gene_1[i] = deepcopy(ind1_MX[i,:BigMatrixSize-NofOutConns])
	
	#Extract sexual material - InterConns - Individual_2	
	for i in range(len(InterConns_gene_2)):
		InterConns_gene_2[i] = deepcopy(ind2_MX[i,:BigMatrixSize-NofOutConns])	
		
	#--- CROSSOVER ---#
	#DUMB crossover - EVERY 2nd ROW is exchanged
	#Exchange sexual material - OutConns
	for i in range(0,len(OutConns_gene_1),2):
		temp = deepcopy(OutConns_gene_1[i])
		OutConns_gene_1[i] = deepcopy(OutConns_gene_2[i])
		OutConns_gene_2[i] = deepcopy(temp)
	
	#Exchange sexual material
	for i in range(0, BigMatrixSize-NofOutConns, 2):
		temp = deepcopy(InterConns_gene_1[i])
		InterConns_gene_1[i] = deepcopy(InterConns_gene_2[i])
		InterConns_gene_2[i] = deepcopy(temp)
	#--- END of CROSSOVER ---#
	
	#Rebuild a big connections matrix out of genes...
	child1_MX=emptyBigMatrix()
	child2_MX=emptyBigMatrix()
	
	rowsR,columnsR,columnsC,rowsC = sortedNonZeroIndices(OutConns_gene_1)
	for i in range(len(rowsR)):
		child1_MX[columnsR[i], BigMatrixSize-(rowsR[i]+1)] = True
		
	rowsR,columnsR,columnsC,rowsC = sortedNonZeroIndices(OutConns_gene_2)		
	for i in range(len(rowsR)):
		child2_MX[columnsR[i], BigMatrixSize-(rowsR[i]+1)] = True		

	rowsR,columnsR,columnsC,rowsC = sortedNonZeroIndices(InterConns_gene_1)		
	for i in range(len(rowsC)):
		child1_MX[rowsC[i], columnsC[i]] = True
	
	rowsR,columnsR,columnsC,rowsC = sortedNonZeroIndices(InterConns_gene_2)	
	for i in range(len(rowsC)):
		child2_MX[rowsC[i], columnsC[i]] = True
		
	child1 = circuit(child1_MX, ind1.ValueVector)
	child2 = circuit(child2_MX, ind2.ValueVector)
	return child1, child2	
	
def halfhalfCrossover(ind1, ind2):
	"""TODO."""
	ind1_MX = deepcopy(np.array(ind1.BigCircuitMatrix))
	ind2_MX = deepcopy(np.array(ind2.BigCircuitMatrix))
	
	OutConns_gene_1 = np.zeros([NofOutConns,BigMatrixSize-NofOutConns], dtype=bool)
	OutConns_gene_2 = np.zeros([NofOutConns,BigMatrixSize-NofOutConns], dtype=bool)
	InterConns_gene_1 = np.zeros([BigMatrixSize-NofOutConns, BigMatrixSize-NofOutConns], dtype=bool)
	InterConns_gene_2 = np.zeros([BigMatrixSize-NofOutConns, BigMatrixSize-NofOutConns], dtype=bool)
	
	#Extract sexual material - OutConns - Individual_1
	for i in range(len(OutConns_gene_1)):
		OutConns_gene_1[i] = ind1_MX[0:(BigMatrixSize-NofOutConns),(BigMatrixSize-(i+1))]
	
	#Extract sexual material - OutConns - Individual_2
	for i in range(len(OutConns_gene_2)):
		OutConns_gene_2[i] = ind2_MX[0:(BigMatrixSize-NofOutConns),(BigMatrixSize-(i+1))]
		
	#Extract sexual material - InterConns - Individual_1	
	for i in range(len(InterConns_gene_1)):
		InterConns_gene_1[i] = ind1_MX[i,:(BigMatrixSize-NofOutConns)]
	
	#Extract sexual material - InterConns - Individual_2	
	for i in range(len(InterConns_gene_2)):
		InterConns_gene_2[i] = ind2_MX[i,:(BigMatrixSize-NofOutConns)]
	
	#--- CROSSOVER ---#
	#Exchange sexual material - OutConns
	#Vertically sliced half of OuterConns_gene is exchanged
	temp = deepcopy(OutConns_gene_1[:,:int(np.floor((BigMatrixSize-NofOutConns)/2))])
	OutConns_gene_1[:,:int(np.floor((BigMatrixSize-NofOutConns)/2))] = deepcopy(OutConns_gene_2[:,:int(np.floor((BigMatrixSize-NofOutConns)/2))])
	OutConns_gene_2[:,:int(np.floor((BigMatrixSize-NofOutConns)/2))] = deepcopy(temp)
	
	#Exchange sexual material - InterConns
	#Modified crossover - upper triangle of upper triangular matrix is exchanged (half-half)
	temp = np.zeros([BigMatrixSize-NofOutConns,BigMatrixSize-NofOutConns], dtype=bool)
	for i in range(0, int(np.ceil((BigMatrixSize-(NofOutConns+1))/2))):
		temp[i,i:(BigMatrixSize-(NofOutConns+1)-i)] = deepcopy(InterConns_gene_1[i,i:(BigMatrixSize-(NofOutConns+1)-i)])
	for i in range(0, int(np.ceil((BigMatrixSize-(NofOutConns+1))/2))):
		InterConns_gene_1[i,i:(BigMatrixSize-(NofOutConns+1)-i)] = deepcopy(InterConns_gene_2[i,i:(BigMatrixSize-(NofOutConns+1)-i)])
	for i in range(0, int(np.ceil((BigMatrixSize-(NofOutConns+1))/2))):
		InterConns_gene_2[i,i:(BigMatrixSize-(NofOutConns+1)-i)] = deepcopy(temp[i,i:(BigMatrixSize-(NofOutConns+1)-i)])

	#--- END of CROSSOVER ---#	
	#Rebuild a big connections matrix out of genes...
	child1_MX=emptyBigMatrix()
	child2_MX=emptyBigMatrix()
	
	
	#Security - if any empty matrix, return parents!
	if (checkNofOnes(OutConns_gene_1) < 1) | (checkNofOnes(OutConns_gene_2) < 1) | (checkNofOnes(InterConns_gene_1) < 1) | (checkNofOnes(InterConns_gene_2) < 1):
	  print "Error: matrix or part of matrix is empty. Cannot return sorted non zero indices."
	  return ind1, ind2
	else:
	  rowsR,columnsR,columnsC,rowsC = sortedNonZeroIndices(OutConns_gene_1)
	  for i in range(len(rowsR)):
		  child1_MX[columnsR[i], BigMatrixSize-(rowsR[i]+1)] = True
		  
	  rowsR,columnsR,columnsC,rowsC = sortedNonZeroIndices(OutConns_gene_2)		
	  for i in range(len(rowsR)):
		  child2_MX[columnsR[i], BigMatrixSize-(rowsR[i]+1)] = True		

	  rowsR,columnsR,columnsC,rowsC = sortedNonZeroIndices(InterConns_gene_1)		
	  for i in range(len(rowsC)):
		  child1_MX[rowsC[i], columnsC[i]] = True
	  
	  rowsR,columnsR,columnsC,rowsC = sortedNonZeroIndices(InterConns_gene_2)	
	  for i in range(len(rowsC)):
		  child2_MX[rowsC[i], columnsC[i]] = True
		  
	  child1 = circuit(child1_MX, ind1.ValueVector)
	  child2 = circuit(child2_MX, ind2.ValueVector)
	  return child1, child2	

def sliceCrossover(ind1, ind2):
  "To write, yea."

  ind1_MX = deepcopy(np.array(ind1.BigCircuitMatrix))
  ind2_MX = deepcopy(np.array(ind2.BigCircuitMatrix))
  
  OutConns_gene_1 = np.zeros([NofOutConns,BigMatrixSize-NofOutConns], dtype=bool)
  OutConns_gene_2 = np.zeros([NofOutConns,BigMatrixSize-NofOutConns], dtype=bool)
  InterConns_gene_1 = np.zeros([BigMatrixSize-NofOutConns, BigMatrixSize-NofOutConns], dtype=bool)
  InterConns_gene_2 = np.zeros([BigMatrixSize-NofOutConns, BigMatrixSize-NofOutConns], dtype=bool)
  
  #Extract sexual material - OutConns - Individual_1
  for i in range(len(OutConns_gene_1)):
	  OutConns_gene_1[i] = ind1_MX[0:(BigMatrixSize-NofOutConns),(BigMatrixSize-(i+1))]
  
  #Extract sexual material - OutConns - Individual_2
  for i in range(len(OutConns_gene_2)):
	  OutConns_gene_2[i] = ind2_MX[0:(BigMatrixSize-NofOutConns),(BigMatrixSize-(i+1))]
	  
  #Extract sexual material - InterConns - Individual_1	
  for i in range(len(InterConns_gene_1)):
	  InterConns_gene_1[i] = ind1_MX[i,:(BigMatrixSize-NofOutConns)]
  
  #Extract sexual material - InterConns - Individual_2	
  for i in range(len(InterConns_gene_2)):
	  InterConns_gene_2[i] = ind2_MX[i,:(BigMatrixSize-NofOutConns)]
	  
  #--- CROSSOVER ---#
  #Exchange sexual material - OutConns
  choice = random.randint(0,1)
  if choice == 1:
    slicePoint = random.randint(1,len(OutConns_gene_1)-1)
    #print "Vertical slice point is: ", slicePoint
    
    #Horizontally sliced part of OuterConns_gene is exchanged
    temp = deepcopy(OutConns_gene_1[0:slicePoint])
    OutConns_gene_1[0:slicePoint] = deepcopy(OutConns_gene_2[0:slicePoint])
    OutConns_gene_2[0:slicePoint] = deepcopy(temp)
  
  #Exchange sexual material - InterConns
  
  if choice == 0:
    slicePoint = random.randint(1,len(InterConns_gene_1)-1)
    #print "Horizontal slice point is: ", slicePoint
    ##Horizontally sliced part of InterConns_gene is exchanged
    temp = deepcopy(InterConns_gene_1[0:slicePoint])
    InterConns_gene_1[0:slicePoint] = deepcopy(InterConns_gene_2[0:slicePoint])
    InterConns_gene_2[0:slicePoint] = deepcopy(temp)

    #--- END of CROSSOVER ---#	  
  
  #Rebuild a big connections matrix out of genes...
  child1_MX=emptyBigMatrix()
  child2_MX=emptyBigMatrix()
  
  rowsR,columnsR,columnsC,rowsC = sortedNonZeroIndices(OutConns_gene_1)
  for i in range(len(rowsR)):
	  child1_MX[columnsR[i], BigMatrixSize-(rowsR[i]+1)] = True
	  
  rowsR,columnsR,columnsC,rowsC = sortedNonZeroIndices(OutConns_gene_2)		
  for i in range(len(rowsR)):
	  child2_MX[columnsR[i], BigMatrixSize-(rowsR[i]+1)] = True		

  rowsR,columnsR,columnsC,rowsC = sortedNonZeroIndices(InterConns_gene_1)		
  for i in range(len(rowsC)):
	  child1_MX[rowsC[i], columnsC[i]] = True
  
  rowsR,columnsR,columnsC,rowsC = sortedNonZeroIndices(InterConns_gene_2)	
  for i in range(len(rowsC)):
	  child2_MX[rowsC[i], columnsC[i]] = True
  
  #ind1_ValueVector, ind2_ValueVector =  valueReproduction(ind1.ValueVector, ind2.ValueVector, crossover="sliceCrossover")
  
  child1 = circuit(child1_MX, copy(ind1.ValueVector))
  child2 = circuit(child2_MX, copy(ind2.ValueVector))
  return child1, child2	  

def lineExchangeCrossover(ind1, ind2):
  """Exchange a random line in a matrix"""
  ind1_MX = deepcopy(np.array(ind1.BigCircuitMatrix))
  ind2_MX = deepcopy(np.array(ind2.BigCircuitMatrix))
  
  OutConns_gene_1 = np.zeros([NofOutConns,BigMatrixSize-NofOutConns], dtype=bool)
  OutConns_gene_2 = np.zeros([NofOutConns,BigMatrixSize-NofOutConns], dtype=bool)
  InterConns_gene_1 = np.zeros([BigMatrixSize-NofOutConns, BigMatrixSize-NofOutConns], dtype=bool)
  InterConns_gene_2 = np.zeros([BigMatrixSize-NofOutConns, BigMatrixSize-NofOutConns], dtype=bool)
  
  #Extract sexual material - OutConns - Individual_1
  for i in range(len(OutConns_gene_1)):
	  OutConns_gene_1[i] = ind1_MX[0:(BigMatrixSize-NofOutConns),(BigMatrixSize-(i+1))]
  
  #Extract sexual material - OutConns - Individual_2
  for i in range(len(OutConns_gene_2)):
	  OutConns_gene_2[i] = ind2_MX[0:(BigMatrixSize-NofOutConns),(BigMatrixSize-(i+1))]
	  
  #Extract sexual material - InterConns - Individual_1	
  for i in range(len(InterConns_gene_1)):
	  InterConns_gene_1[i] = ind1_MX[i,:(BigMatrixSize-NofOutConns)]
  
  #Extract sexual material - InterConns - Individual_2	
  for i in range(len(InterConns_gene_2)):
	  InterConns_gene_2[i] = ind2_MX[i,:(BigMatrixSize-NofOutConns)]
	  
  #--- CROSSOVER ---#
  #Exchange sexual material - OutConns
  innerNodeMatingProb = 0.5
  if random.uniform(0,1) < innerNodeMatingProb:
    lineIndex = random.randint(1,len(OutConns_gene_1)-1)   
    #Horizontally sliced part of OuterConns_gene is exchanged
    temp = deepcopy(OutConns_gene_1[lineIndex])
    OutConns_gene_1[lineIndex] = deepcopy(OutConns_gene_2[lineIndex])
    OutConns_gene_2[lineIndex] = deepcopy(temp)
  #Exchange sexual material - InterConns
  else:
    lineIndex = random.randint(1,len(InterConns_gene_1)-1)
    #print "lineIndex is(InterConns): ", lineIndex
    ##Horizontally sliced part of InterConns_gene is exchanged
    temp = deepcopy(InterConns_gene_1[lineIndex])
    InterConns_gene_1[lineIndex] = deepcopy(InterConns_gene_2[lineIndex])
    InterConns_gene_2[lineIndex] = deepcopy(temp)

    #--- END of CROSSOVER ---#	  
  
  #Rebuild a big connections matrix out of genes...
  child1_MX=emptyBigMatrix()
  child2_MX=emptyBigMatrix()
  
  rowsR,columnsR,columnsC,rowsC = sortedNonZeroIndices(OutConns_gene_1)
  for i in range(len(rowsR)):
	  child1_MX[columnsR[i], BigMatrixSize-(rowsR[i]+1)] = True
	  
  rowsR,columnsR,columnsC,rowsC = sortedNonZeroIndices(OutConns_gene_2)		
  for i in range(len(rowsR)):
	  child2_MX[columnsR[i], BigMatrixSize-(rowsR[i]+1)] = True		

  rowsR,columnsR,columnsC,rowsC = sortedNonZeroIndices(InterConns_gene_1)		
  for i in range(len(rowsC)):
	  child1_MX[rowsC[i], columnsC[i]] = True
  
  rowsR,columnsR,columnsC,rowsC = sortedNonZeroIndices(InterConns_gene_2)	
  for i in range(len(rowsC)):
	  child2_MX[rowsC[i], columnsC[i]] = True
  
  #ind1_ValueVector, ind2_ValueVector =  valueReproduction(ind1.ValueVector, ind2.ValueVector, crossover="rangeCrossover")
  child1 = circuit(child1_MX, copy(ind1.ValueVector))
  child2 = circuit(child2_MX, copy(ind2.ValueVector))
  return child1, child2	    
  
def pinViseCrossover(ind1, ind2):
  """Exchange a random line in a matrix. In Inner connection matrix it exchanges a horisontal AND vertical line."""
  ind1_MX = deepcopy(np.array(ind1.BigCircuitMatrix))
  ind2_MX = deepcopy(np.array(ind2.BigCircuitMatrix))
  
  OutConns_gene_1 = np.zeros([NofOutConns,BigMatrixSize-NofOutConns], dtype=bool)
  OutConns_gene_2 = np.zeros([NofOutConns,BigMatrixSize-NofOutConns], dtype=bool)
  InterConns_gene_1 = np.zeros([BigMatrixSize-NofOutConns, BigMatrixSize-NofOutConns], dtype=bool)
  InterConns_gene_2 = np.zeros([BigMatrixSize-NofOutConns, BigMatrixSize-NofOutConns], dtype=bool)
  
  #Extract sexual material - OutConns - Individual_1
  for i in range(len(OutConns_gene_1)):
	  OutConns_gene_1[i] = ind1_MX[0:(BigMatrixSize-NofOutConns),(BigMatrixSize-(i+1))]
  
  #Extract sexual material - OutConns - Individual_2
  for i in range(len(OutConns_gene_2)):
	  OutConns_gene_2[i] = ind2_MX[0:(BigMatrixSize-NofOutConns),(BigMatrixSize-(i+1))]
	  
  #Extract sexual material - InterConns - Individual_1	
  for i in range(len(InterConns_gene_1)):
	  InterConns_gene_1[i] = ind1_MX[i,:(BigMatrixSize-NofOutConns)]
  
  #Extract sexual material - InterConns - Individual_2	
  for i in range(len(InterConns_gene_2)):
	  InterConns_gene_2[i] = ind2_MX[i,:(BigMatrixSize-NofOutConns)]
	  
  #--- CROSSOVER ---#
  #Exchange sexual material - OutConns
  innerNodeMatingProb = 0.5
  if random.uniform(0,1) < innerNodeMatingProb:	
    lineIndex = random.randint(1,len(OutConns_gene_1)-1)   
    #Horizontally sliced part of OuterConns_gene is exchanged
    temp = deepcopy(OutConns_gene_1[lineIndex])
    OutConns_gene_1[lineIndex] = deepcopy(OutConns_gene_2[lineIndex])
    OutConns_gene_2[lineIndex] = deepcopy(temp)
    
  #Exchange sexual material - InterConns
  else:
    lineIndex = random.randint(1,len(InterConns_gene_1)-1)
    #print "lineIndex is(InterConns): ", lineIndex
    ##Horizontally sliced part of InterConns_gene is exchanged
    tempH = deepcopy(InterConns_gene_1[lineIndex])
    InterConns_gene_1[lineIndex] = deepcopy(InterConns_gene_2[lineIndex])
    InterConns_gene_2[lineIndex] = deepcopy(tempH)
    
    ##Vertically sliced part of InterConns_gene is exchanged TODO
    tempV = deepcopy(InterConns_gene_1[:,lineIndex])
    InterConns_gene_1[:,lineIndex] = deepcopy(InterConns_gene_2[:,lineIndex])
    InterConns_gene_2[:,lineIndex] = deepcopy(tempV)

    #--- END of CROSSOVER ---#	  
  
  #Rebuild a big connections matrix out of genes...
  child1_MX=emptyBigMatrix()
  child2_MX=emptyBigMatrix()
  
  rowsR,columnsR,columnsC,rowsC = sortedNonZeroIndices(OutConns_gene_1)
  for i in range(len(rowsR)):
	  child1_MX[columnsR[i], BigMatrixSize-(rowsR[i]+1)] = True
	  
  rowsR,columnsR,columnsC,rowsC = sortedNonZeroIndices(OutConns_gene_2)		
  for i in range(len(rowsR)):
	  child2_MX[columnsR[i], BigMatrixSize-(rowsR[i]+1)] = True		

  rowsR,columnsR,columnsC,rowsC = sortedNonZeroIndices(InterConns_gene_1)		
  for i in range(len(rowsC)):
	  child1_MX[rowsC[i], columnsC[i]] = True
  
  rowsR,columnsR,columnsC,rowsC = sortedNonZeroIndices(InterConns_gene_2)	
  for i in range(len(rowsC)):
	  child2_MX[rowsC[i], columnsC[i]] = True
  
  #ind1_ValueVector, ind2_ValueVector =  valueReproduction(ind1.ValueVector, ind2.ValueVector, crossover="rangeCrossover")
  child1 = circuit(child1_MX, copy(ind1.ValueVector))
  child2 = circuit(child2_MX, copy(ind2.ValueVector))
  return child1, child2	    


def pinViseCrossover2(ind1, ind2):
  """Exchange a random line in a matrix. In Inner connection matrix it exchanges a horizontal AND vertical line."""
  ind1_MX = deepcopy(np.array(ind1.BigCircuitMatrix))
  ind2_MX = deepcopy(np.array(ind2.BigCircuitMatrix))

	  
  #--- CROSSOVER ---#
  for i in range(0, random.randint(1, 5)):
    #Exchange sexual material - BigCircuitMatrix
    lineIndex = random.randint(1,len(ind1_MX)-1)
    #print "lineIndex is(InterConns): ", lineIndex
    ##Horizontally sliced part of whole matrix is exchanged
    tempH = deepcopy(ind1_MX[lineIndex])
    ind1_MX[lineIndex] = deepcopy(ind2_MX[lineIndex])
    ind2_MX[lineIndex] = deepcopy(tempH)
    
    ##Vertically sliced part of whole matrix is exchanged
    tempV = deepcopy(ind1_MX[:,lineIndex])
    ind1_MX[:,lineIndex] = deepcopy(ind2_MX[:,lineIndex])
    ind2_MX[:,lineIndex] = deepcopy(tempV)
    #--- END of CROSSOVER ---#	  
  #ind1_ValueVector, ind2_ValueVector =  valueReproduction(ind1.ValueVector, ind2.ValueVector, crossover="rangeCrossover")
  rowsR,columnsR,columnsC,rowsC = sortedNonZeroIndices(ind1_MX)
  if (columnsC == None) | (rowsC == None):
    return ind1, ind2
  rowsR,columnsR,columnsC,rowsC = sortedNonZeroIndices(ind2_MX)
  if (columnsC == None) | (rowsC == None):
    return ind1, ind2  
  
  child1 = circuit(ind1_MX, copy(ind1.ValueVector))
  child2 = circuit(ind1_MX, copy(ind2.ValueVector))
  return child1, child2	    


	
def mutationADDnode(ind1):
	"""TODO."""
	ind1_MX = deepcopy(np.array(ind1.BigCircuitMatrix))
	OutConns_gene_1 = np.zeros([NofOutConns,BigMatrixSize-NofOutConns], dtype=bool)
	InterConns_gene_1 = np.zeros([BigMatrixSize-NofOutConns, BigMatrixSize-NofOutConns], dtype=bool)

	#Extract sexual material - OutConns - Individual_1
	for i in range(len(OutConns_gene_1)):
		OutConns_gene_1[i] = ind1_MX[0:BigMatrixSize-NofOutConns,BigMatrixSize-(i+1)]	

	#Extract sexual material - InterConns - Individual_1	
	for i in range(len(InterConns_gene_1)):
		InterConns_gene_1[i] = ind1_MX[i,:BigMatrixSize-NofOutConns]

	#--- addMUTATION ---#
	#Pick a random gene. Put a random connection. 
	
	choice = random.randint(0,1)	#make a choice between mutation Inter (0) or Outer (1) conns
	done = 0
	SafeCount = 0
	if choice == 0:					#Mutate interconns
		while done == 0:				#Be sure to do a change
			i = random.randint(0,len(InterConns_gene_1)-1)
			j = random.randint(i,len(InterConns_gene_1[i])-1)
			SafeCount = SafeCount + 1
			#print "mutate interconns"
			#print "i=", i, " j=", j
			if 	InterConns_gene_1[i][j] != True:
				InterConns_gene_1[i][j] = True
				done = 1
			if SafeCount > 500 :
				done = 1
	if choice == 1:					#Mutate outerconns
		while done == 0:				#Be sure to do a change
			i = random.randint(0,len(OutConns_gene_1)-1)
			j = random.randint(0,len(OutConns_gene_1[i])-1)
			SafeCount = SafeCount + 1
			if 	OutConns_gene_1[i][j] != True:
				OutConns_gene_1[i][j] = True
				done = 1
			if SafeCount > 500 :
				done = 1
	#--- END of addMUTATION ---#
	
	#Rebuild a big connections matrix out of genes...
	child1_MX=emptyBigMatrix()

	rowsR,columnsR,columnsC,rowsC = sortedNonZeroIndices(OutConns_gene_1)
	for i in range(len(rowsR)):
		child1_MX[columnsR[i], BigMatrixSize-(rowsR[i]+1)] = True
		
	rowsR,columnsR,columnsC,rowsC = sortedNonZeroIndices(InterConns_gene_1)		
	for i in range(len(rowsC)):
		child1_MX[rowsC[i], columnsC[i]] = True

	child1 = circuit(child1_MX, ind1.ValueVector)
	return child1
	
def mutationREMOVEnode(ind1):
	"""TODO."""
	ind1_MX = deepcopy(np.array(ind1.BigCircuitMatrix))
	OutConns_gene_1 = np.zeros([NofOutConns,BigMatrixSize-NofOutConns], dtype=bool)
	InterConns_gene_1 = np.zeros([BigMatrixSize-NofOutConns, BigMatrixSize-NofOutConns], dtype=bool)

	#Extract sexual material - OutConns - Individual_1
	for i in range(len(OutConns_gene_1)):
		OutConns_gene_1[i] = ind1_MX[0:BigMatrixSize-NofOutConns,BigMatrixSize-(i+1)]	

	#Extract sexual material - InterConns - Individual_1	
	for i in range(len(InterConns_gene_1)):
		InterConns_gene_1[i] = ind1_MX[i,:BigMatrixSize-NofOutConns]

	#--- addMUTATION ---#
	#Pick a random gene. Delete a random connection. 
	#Mutate ONLY INTERCONNS. Every outer conn must have at least one connection.
	done = 0
	SafeCount = 0
	while done == 0:				#Be sure to do a change
			rowsR,columnsR,columnsC,rowsC = sortedNonZeroIndices(InterConns_gene_1)
			i = random.randint(0, len(rowsR)-1)
			j = columnsR[i]
			SafeCount = SafeCount + 1
			if rowsR[i] != columnsR[i]:
				InterConns_gene_1[rowsR[i]][columnsR[i]] = False
				done = 1
			if SafeCount > 500 :
				done = 1
	#--- END of addMUTATION ---#
	
	#Rebuild a big connections matrix out of genes...
	child1_MX=emptyBigMatrix()	
	
	rowsR,columnsR,columnsC,rowsC = sortedNonZeroIndices(OutConns_gene_1)
	if (columnsC == None) | (rowsC == None):
	  return ind1
	for i in range(len(rowsR)):
		child1_MX[columnsR[i], BigMatrixSize-(rowsR[i]+1)] = True
		
	rowsR,columnsR,columnsC,rowsC = sortedNonZeroIndices(InterConns_gene_1)
	if (columnsC == None) | (rowsC == None):
	  return ind1
	for i in range(len(rowsC)):
		child1_MX[rowsC[i], columnsC[i]] = True

	child1 = circuit(child1_MX, ind1.ValueVector)
	return child1	
	

def makeNetlist(circuit, generationNum, individualNum, FullBigCircuitMatrix):
	"""Make netlist from circuit object"""
	BigCircuitMatrix = circuit.BigCircuitMatrix
	ValueVector = circuit.ValueVector
	HOTCIRC = "HOT_CIRCUIT"		#Name of the current netlist
 
	nodesBank = np.zeros(shape=(BigMatrixSize),dtype=int)	#empty 1D container of node numbers for each element pins
	OutConnsBank = {}	#will contain where nodes GND, Vp, Vn and so on are meant to be
	
	#REPAIR BigCircuitMatrix -- ADD FULL REDUNDANCE
	FullBigCircuitMatrix = circuit.fullRedundancyMatrix
	
	#Locations of non zero elements
	rowsR,columnsR,columnsC,rowsC = sortedNonZeroIndices(FullBigCircuitMatrix)
	#matrixDensity = float(len(rowsR))/float((BigMatrixSize*BigMatrixSize/2))	#(ones/(all/2))
	#print "Matrix density :", matrixDensity
	#matrixQuaziID = sum(rowsR)+sum(columnsR)-BigMatrixSize*(BigMatrixSize-1)
	#print "Matrix Quazi ID: ", matrixQuaziID
	
	#if sum(i > (BigMatrixSize - NofOutConns) for i in rowsR) > NofOutConns:
	#	print "Gresniki=", sum(i > (BigMatrixSize - NofOutConns) for i in rowsR)
	#	raw_input("Pritisni!")
	
	#make nodes bank with just looking "up" in BigCircuitMatrix
	for i in range(0,BigMatrixSize):
		up = np.where(columnsC == i)
		right = np.where(rowsR == i)
		#look up
		for j in up:
			# Commented part was written with element-not-used function in mind. 
			#if (len(up[0]) > 1) | (len(right[0]) > 1):
			#	nodesBank[i] = min(rowsC[j])+1 	#increment the node number to keep 0 for Spice GND
			#else:
			#	nodesBank[i] = 0	#Here 0 is used as a marker to put nothing there (elemnt pin not used).
			
			#Let's assume every given element must be used. (A node might get unused - dangling)
			#Yes, but though, both SpiceOpus and HSpice hangs in some occasions when there is a dangling node. Need to fix this. 
			nodesBank[i] = min(rowsC[j])+1 	#increment the node number to keep 0 for Spice GND
		if i > (BigMatrixSize-1-NofOutConns):
			OutConnsBank[i] = nodesBank[i]
	
	#print "OutConnsBank: ",OutConnsBank
	
	#Open empty file for netlist
	ID = "g_" + str(generationNum) + "_i_" + str(individualNum)
	netlistName = ID + "_subckt.cir"
	circ = open(ID + "_subckt.cir", "w")
	circ.write("*%s \n" %netlistName)

	circ.write(".SUBCKT %s %d " %(HOTCIRC, OutConnsBank[BigMatrixSize-NofOutConns]))	#GND
	circ.write("%d " %(OutConnsBank[BigMatrixSize-NofOutConns+1]))				#Vsupp
	#circ.write("%d " %(OutConnsBank[BigMatrixSize-NofOutConns+2]))				#Iref
	circ.write("%d " %(OutConnsBank[BigMatrixSize-NofOutConns+2]))				#Vout
	circ.write("vsupp ")									#vsupp and vsupn are global nodes - so no need here 
	circ.write("vsupn ")
	circ.write("\n")
	#circ.write("vdd \n")
	#circ.write("%d " %(OutConnsBank[BigMatrixSize-NofOutConns+3]))
	#circ.write("%d " %(OutConnsBank[BigMatrixSize-NofOutConns+4]))
	#circ.write("%d " %(OutConnsBank[BigMatrixSize-NofOutConns+5]))
	#circ.write("%d \n" %(OutConnsBank[BigMatrixSize-NofOutConns+6]))

	#---------------------------	
	#---------------------------	

	value_index = 0		#goes with ValueVector which contains element values
	multipl_index = 0	#goes with Multipl, which contains multiplication factors for mos (only integer)
	#go thru all lines except connections (GND, Vp, Vn, Vin1, Vin2, Vout1, Vout2 - they shall not connect to each other!)
	i=0
	while i < BigMatrixSize:
		#Here, every used element is written into a netlist and given nodes and values/models
		#if-statements must strictly follow the sequence of elements in BigCircuitMatrix.
		if nodesBank[i] != 0:
			if i < np.sum(ALLPINS[0:1]): #NOTE: see globalVars.py for ALLPINS
				#put resistor------------------------------------------------------------------------
				if (nodesBank[i]==nodesBank[i+1]) | (sum(nodesBank == nodesBank[i]) < 2):
				  circ.write("*") #exclude element from netlist if element pins connected together or if dangeling node detected
				circ.write("R_%s %s " %(i, nodesBank[i]))	#put name of element and + node
				circ.write("%s %s \n" %(nodesBank[i+1], ValueVector[value_index]))	#put - node and value
				value_index = value_index + 1									#increment value_index
				i = i+2
				#-------------------------------------------------------------------------------------
			elif i < np.sum(ALLPINS[0:2]):
				#put capacitor------------------------------------------------------------------------
				if (nodesBank[i]==nodesBank[i+1]) | (sum(nodesBank == nodesBank[i]) < 2):
				  circ.write("*") #exclude element from netlist
				circ.write("C_%s %s " %(i, nodesBank[i]))	#put name of element and + node
				circ.write("%s %s \n" %(nodesBank[i+1], ValueVector[value_index]))	#put - node and value
				value_index = value_index + 1									#increment value_index
				i = i+2
				#-------------------------------------------------------------------------------------
			elif i < np.sum(ALLPINS[0:3]):
				#put inductor------------------------------------------------------------------------
				if (nodesBank[i]==nodesBank[i+1]) | (sum(nodesBank == nodesBank[i]) < 2):
				  circ.write("*") #exclude element from netlist
				circ.write("L_%s %s " %(i, nodesBank[i]))	#put name of element and + node
				circ.write("%s %s \n" %(nodesBank[i+1], ValueVector[value_index]))	#put - node and value
				value_index = value_index + 1									#increment value_index
				i = i+2
				#-------------------------------------------------------------------------------------
			elif i < np.sum(ALLPINS[0:4]):
				#put zener------------------------------------------------------------------------
				if (nodesBank[i]==nodesBank[i+1]) | (sum(nodesBank == nodesBank[i]) < 2):
				  circ.write("*") #exclude element from netlist
				circ.write("dz_%s %s " %(i, nodesBank[i]))	#put name of element and + node
				circ.write("%s %s \n" %(nodesBank[i+1], MODEL_ZENER))	#put - node and model name
				i = i+2
				#-------------------------------------------------------------------------------------
			elif i < np.sum(ALLPINS[0:5]):
				#put NPN------------------------------------------------------------------------
				if (len(set([nodesBank[i], nodesBank[i+1], nodesBank[i+2]]))==1) | (sum(nodesBank == nodesBank[i]) < 2):
				  circ.write("*") #exclude element from netlist
				circ.write("q_%s %s " %(i, nodesBank[i]))	#put name of element and COLLECTOR node
				circ.write("%s " %(nodesBank[i+1]))			#put BASE node
				circ.write("%s %s \n" %(nodesBank[i+2], MODEL_NPN))	#put EMITTER node and model name
				i = i+3
				#----------------------------------------------------------------------------	
			elif i < np.sum(ALLPINS[0:6]):
				#put 3 parallel PNPs (subckt)------------------------------------------------
				if (len(set([nodesBank[i], nodesBank[i+1], nodesBank[i+2]]))==1) | (sum(nodesBank == nodesBank[i]) < 2):
				  circ.write("*") #exclude element from netlist
				circ.write("x_%s %s " %(i, nodesBank[i]))	#put name of element and COLLECTOR node
				circ.write("%s " %(nodesBank[i+1]))			#put BASE node
				circ.write("%s %s \n" %(nodesBank[i+2], "par3pnp"))	#put EMITTER node and subckt name
				i = i+3
				#----------------------------------------------------------------------------	
			elif i < np.sum(ALLPINS[0:7]):
				#put PNP------------------------------------------------------------------------
				if (len(set([nodesBank[i], nodesBank[i+1], nodesBank[i+2]]))==1) | (sum(nodesBank == nodesBank[i]) < 2):
				  circ.write("*") #exclude element from netlist
				circ.write("q_%s %s " %(i, nodesBank[i]))	#put name of element and COLLECTOR node
				circ.write("%s " %(nodesBank[i+1]))			#put BASE node
				circ.write("%s %s \n" %(nodesBank[i+2], MODEL_PNP))	#put EMITTER node and model name
				i = i+3
				#----------------------------------------------------------------------------	
			
			elif i < np.sum(ALLPINS[0:8]):
				#put NMOS------------------------------------------------------------------------
				#if (len(set([nodesBank[i], nodesBank[i+1], nodesBank[i+2]]))==1) | (sum(nodesBank == nodesBank[i]) < 2):
				#  circ.write("*") #exclude element from netlist
				circ.write("xmn_%s %s " %(i, nodesBank[i]))	#put drain
				circ.write("%s " %(nodesBank[i+1]))		#put gate
				circ.write("%s %s %s " %(nodesBank[i+2], "0", MODEL_NMOS))	#put source, bulk(to vss), and model name
				circ.write("w=%s l=%s m=%s \n" %(ValueVector[value_index], ValueVector[value_index+1], Multipl[multipl_index]))
				value_index = value_index + 2
				i = i+3
				multipl_index = multipl_index + 1
			
			elif i < np.sum(ALLPINS[0:9]):
				#put PMOS------------------------------------------------------------------------
				#if (len(set([nodesBank[i], nodesBank[i+1], nodesBank[i+2]]))==1) | (sum(nodesBank == nodesBank[i]) < 2):
				#  circ.write("*") #exclude element from netlist
				circ.write("xmp_%s %s " %(i, nodesBank[i]))	#put drain
				circ.write("%s " %(nodesBank[i+1]))		#put gate
				circ.write("%s %s %s " %(nodesBank[i+2], "vdd", MODEL_PMOS))	#put source, bulk(to vdd), and model name
				circ.write("w=%s l=%s m=%s \n" %(ValueVector[value_index], ValueVector[value_index+1], Multipl[multipl_index]))
				value_index = value_index + 2
				i = i+3
				multipl_index = multipl_index + 1
			
			elif i < np.sum(ALLPINS[0:10]):
				#put OPAMP1------------------------------------------------------------------------
				if (len(set([nodesBank[i], nodesBank[i+1], nodesBank[i+2]]))==1) | (sum(nodesBank == nodesBank[i]) < 2):
				  circ.write("*") #exclude element from netlist
				circ.write("x_%s %s " %(i, nodesBank[i]))	#put name of element and NON-INVERTING node 
				circ.write("%s " %(nodesBank[i+1]))		#put INVERTING node
				circ.write("vsupp ")				#put POSITIVE POWER SUPPLY node NOTE: Change according to topdc.cir circuit
				circ.write("vsupn ")				#put NEGATIVE POWER SUPPLY node NOTE: Change according to topdc.cir circuit
				circ.write("%s %s \n" %(nodesBank[i+2], MODEL_OPAMP1))	#put OUTPUT node and model name
				i = i+3
			
			#Writing various subcircuits to netlist:
			elif i < np.sum(ALLPINS[0:11]):#UNDER CONSTRUCTION from here on
				#put PMosCurrSrc1Stg------------------------------------------------------------------
				if (len(set([nodesBank[i], nodesBank[i+1], nodesBank[i+2]]))==1) | (sum(nodesBank == nodesBank[i]) < 2):
				  circ.write("*") #exclude element from netlist
				circ.write("x_%s %s " %(i, nodesBank[i]))	#put name of element and 1st node 
				circ.write("%s " %(nodesBank[i+1]))		#put 2nd node
				circ.write("%s " %(nodesBank[i+2]))		#put 3rd node
				circ.write("PMosCurrSrc1stg ")			#put SUBCIRCUIT name
				circ.write("w=%s l=%s \n" %(ValueVector[value_index], ValueVector[value_index+1]))
				value_index = value_index + 2
				i = i+3
			  
			elif i < np.sum(ALLPINS[0:12]):
				#put NofCascPMosCurrSrc1stg------------------------------------------------------------------
				if (len(set([nodesBank[i], nodesBank[i+1], nodesBank[i+2]]))==1) | (sum(nodesBank == nodesBank[i]) < 2):
				  circ.write("*") #exclude element from netlist
				circ.write("x_%s %s " %(i, nodesBank[i]))	#put name of element and 1st node 
				circ.write("%s " %(nodesBank[i+1]))		#put 2nd node
				circ.write("%s " %(nodesBank[i+2]))				#put 3rd node
				circ.write("CascPMosCurrSrc1stg ")		#put SUBCIRCUIT name
				circ.write("w=%s l=%s \n" %(ValueVector[value_index], ValueVector[value_index+1]))
				value_index = value_index + 2
				i = i+3
			  
			elif i < np.sum(ALLPINS[0:13]):
				#put NofCascPMosCurrSrc1stg------------------------------------------------------------------
				if (len(set([nodesBank[i], nodesBank[i+1], nodesBank[i+2]]))==1) | (sum(nodesBank == nodesBank[i]) < 2):
				  circ.write("*") #exclude element from netlist
				circ.write("x_%s %s " %(i, nodesBank[i]))	#put name of element and 1st node 
				circ.write("%s " %(nodesBank[i+1]))		#put 2nd node
				circ.write("%s " %(nodesBank[i+2]))				#put 3rd node
				circ.write("NMosAmp1ResOnSrc ")		#put SUBCIRCUIT name
				circ.write("w=%s l=%s r=%s \n" %(ValueVector[value_index], ValueVector[value_index+1], ValueVector[value_index+2]))
				value_index = value_index + 3
				i = i+3
			#From here on, put elements with 4 pins:
			elif i < np.sum(ALLPINS[0:14]):
				#put BJTNPNCurrSink------------------------------------------------------------------
				if (len(set([nodesBank[i], nodesBank[i+1], nodesBank[i+2], nodesBank[i+3]]))==1) | (sum(nodesBank == nodesBank[i]) < 2):
				  circ.write("*") #exclude element from netlist
				circ.write("x_%s %s " %(i, nodesBank[i]))	#put name of element and 1st node 
				circ.write("%s " %(nodesBank[i+1]))		#put 2nd node
				circ.write("%s " %(nodesBank[i+2]))		#put 3rd node
				circ.write("%s " %(nodesBank[i+3]))		#put 4rd node
				circ.write("BJTNPNCurrSink \n")			#put SUBCIRCUIT name
				i = i+4	  
			  
			elif i < np.sum(ALLPINS[0:15]):
				#put BJTNPNCurrSink------------------------------------------------------------------
				if (len(set([nodesBank[i], nodesBank[i+1], nodesBank[i+2], nodesBank[i+3]]))==1) | (sum(nodesBank == nodesBank[i]) < 2):
				  circ.write("*") #exclude element from netlist
				circ.write("x_%s %s " %(i, nodesBank[i]))	#put name of element and 1st node 
				circ.write("%s " %(nodesBank[i+1]))		#put 2nd node
				circ.write("%s " %(nodesBank[i+2]))		#put 3rd node
				circ.write("%s " %(nodesBank[i+3]))		#put 4rd node
				circ.write("BJTPNPCurrSrc \n")			#put SUBCIRCUIT name
				i = i+4	   
			  
			elif i < np.sum(ALLPINS[0:16]):
				#put NMosCurrMirr------------------------------------------------------------------
				if (len(set([nodesBank[i], nodesBank[i+1], nodesBank[i+2], nodesBank[i+3]]))==1) | (sum(nodesBank == nodesBank[i]) < 2):
				  circ.write("*") #exclude element from netlist
				circ.write("x_%s %s " %(i, nodesBank[i]))	#put name of element and 1st node 
				circ.write("%s " %(nodesBank[i+1]))		#put 2nd node
				circ.write("%s " %(nodesBank[i+2]))		#put 3rd node
				circ.write("%s " %(nodesBank[i+3]))		#put 4rd node
				circ.write("NMosCurrMirr ")		#put SUBCIRCUIT name
				circ.write("w=%s l=%s \n" %(ValueVector[value_index], ValueVector[value_index+1]))
				value_index = value_index + 2
				i = i+4
			  
			elif i < np.sum(ALLPINS[0:17]):
				#put CascNMosCurrMirr------------------------------------------------------------------
				if (len(set([nodesBank[i], nodesBank[i+1], nodesBank[i+2], nodesBank[i+3]]))==1) | (sum(nodesBank == nodesBank[i]) < 2):
				  circ.write("*") #exclude element from netlist
				circ.write("x_%s %s " %(i, nodesBank[i]))	#put name of element and 1st node 
				circ.write("%s " %(nodesBank[i+1]))		#put 2nd node
				circ.write("%s " %(nodesBank[i+2]))		#put 3rd node
				circ.write("%s " %(nodesBank[i+3]))		#put 4rd node
				circ.write("CascNMosCurrMirr ")		#put SUBCIRCUIT name
				circ.write("w=%s l=%s \n" %(ValueVector[value_index], ValueVector[value_index+1]))
				value_index = value_index + 2
				i = i+4	  
			  
			elif i < np.sum(ALLPINS[0:18]):
				#put PMosCurrSrc2Stg------------------------------------------------------------------
				if (len(set([nodesBank[i], nodesBank[i+1], nodesBank[i+2], nodesBank[i+3]]))==1) | (sum(nodesBank == nodesBank[i]) < 2):
				  circ.write("*") #exclude element from netlist
				circ.write("x_%s %s " %(i, nodesBank[i]))	#put name of element and 1st node 
				circ.write("%s " %(nodesBank[i+1]))		#put 2nd node
				circ.write("%s " %(nodesBank[i+2]))		#put 3rd node
				circ.write("%s " %(nodesBank[i+3]))		#put 4rd node
				circ.write("PMosCurrSrc2Stg ")		#put SUBCIRCUIT name
				circ.write("w1=%s l1=%s w2=%s l2=%s \n" %(ValueVector[value_index], ValueVector[value_index+1], ValueVector[value_index+2], ValueVector[value_index+3]))
				value_index = value_index + 4
				i = i+4	 	  
			  
			elif i < np.sum(ALLPINS[0:19]):
				#put CascPMosCurrSrc2Stg------------------------------------------------------------------
				if (len(set([nodesBank[i], nodesBank[i+1], nodesBank[i+2], nodesBank[i+3]]))==1) | (sum(nodesBank == nodesBank[i]) < 2):
				  circ.write("*") #exclude element from netlist
				circ.write("x_%s %s " %(i, nodesBank[i]))	#put name of element and 1st node 
				circ.write("%s " %(nodesBank[i+1]))		#put 2nd node
				circ.write("%s " %(nodesBank[i+2]))		#put 3rd node
				circ.write("%s " %(nodesBank[i+3]))		#put 4rd node
				circ.write("CascPMosCurrSrc2Stg ")		#put SUBCIRCUIT name
				circ.write("w1=%s l1=%s w2=%s l2=%s \n" %(ValueVector[value_index], ValueVector[value_index+1], ValueVector[value_index+2], ValueVector[value_index+3]))
				value_index = value_index + 4
				i = i+4	 	  
			  
			elif i < np.sum(ALLPINS[0:20]):
				#put PMosCascode------------------------------------------------------------------
				if (len(set([nodesBank[i], nodesBank[i+1], nodesBank[i+2], nodesBank[i+3]]))==1) | (sum(nodesBank == nodesBank[i]) < 2):
				  circ.write("*") #exclude element from netlist
				circ.write("x_%s %s " %(i, nodesBank[i]))	#put name of element and 1st node 
				circ.write("%s " %(nodesBank[i+1]))		#put 2nd node
				circ.write("%s " %(nodesBank[i+2]))		#put 3rd node
				circ.write("%s " %(nodesBank[i+3]))		#put 4rd node
				circ.write("PMosCascode ")			#put SUBCIRCUIT name
				circ.write("w1=%s l1=%s w2=%s l2=%s \n" %(ValueVector[value_index], ValueVector[value_index+1], ValueVector[value_index+2], ValueVector[value_index+3]))
				value_index = value_index + 4
				i = i+4	 	  
			  
			elif i < np.sum(ALLPINS[0:21]):
				#put NMosCascode------------------------------------------------------------------
				if (len(set([nodesBank[i], nodesBank[i+1], nodesBank[i+2], nodesBank[i+3]]))==1) | (sum(nodesBank == nodesBank[i]) < 2):
				  circ.write("*") #exclude element from netlist
				circ.write("x_%s %s " %(i, nodesBank[i]))	#put name of element and 1st node 
				circ.write("%s " %(nodesBank[i+1]))		#put 2nd node
				circ.write("%s " %(nodesBank[i+2]))		#put 3rd node
				circ.write("%s " %(nodesBank[i+3]))		#put 4rd node
				circ.write("NMosCascode ")		#put SUBCIRCUIT name
				circ.write("w1=%s l1=%s w2=%s l2=%s \n" %(ValueVector[value_index], ValueVector[value_index+1], ValueVector[value_index+2], ValueVector[value_index+3]))
				value_index = value_index + 4
				i = i+4	 
			  
			else:
				i = i+1
		else:
			i = i+1
	#convergence aid
	#add 1 GOhm between every node and GND to help HSpice converge if dangling nodes in circuit. Added 13.5.2016. Looks it works fine.
	circ.write("\n*Convergence-aid resistors:\n")
	for i in set(nodesBank):  
	  circ.write("R_%s %s " %(np.max(nodesBank)+1+i, i))	#put name of element and + node
	  circ.write("%s %s \n" %("0", "1e9"))	#put - node and value	
	  
	circ.write(".ends\n")
	circ.close()
	#print "**********************************************************************"
	#print "* Generation %s, individual %s's netlist created. *" %(generationNum, individualNum)
	#print "**********************************************************************"


def makeNetlist2(circuit, generationNum, individualNum, FullBigCircuitMatrix):
	"""Make netlist from circuit object
	  Experimental - to test nodesBank creation without need of fullRedundancyBigCircuitMatrix.
	  Does not work. Do not use.
	"""
	BigCircuitMatrix = circuit.BigCircuitMatrix
	ValueVector = circuit.ValueVector
	HOTCIRC = "HOT_CIRCUIT"		#Name of the current netlist
 
	nodesBank = np.zeros(shape=(BigMatrixSize),dtype=int)	#empty 1D container of node numbers for each element pins
	OutConnsBank = {}	#will contain where nodes GND, Vp, Vn and so on are meant to be
	
	#REPAIR BigCircuitMatrix -- ADD FULL REDUNDANCE
	FullBigCircuitMatrix = circuit.fullRedundancyMatrix
	
	#Locations of non zero elements
	rowsR,columnsR,columnsC,rowsC = sortedNonZeroIndices(FullBigCircuitMatrix)
	#matrixDensity = float(len(rowsR))/float((BigMatrixSize*BigMatrixSize/2))	#(ones/(all/2))
	#print "Matrix density :", matrixDensity
	#matrixQuaziID = sum(rowsR)+sum(columnsR)-BigMatrixSize*(BigMatrixSize-1)
	#print "Matrix Quazi ID: ", matrixQuaziID
	
	#if sum(i > (BigMatrixSize - NofOutConns) for i in rowsR) > NofOutConns:
	#	print "Gresniki=", sum(i > (BigMatrixSize - NofOutConns) for i in rowsR)
	#	raw_input("Pritisni!")
	
	#make nodes bank with just looking "up" in BigCircuitMatrix
	for i in range(0,BigMatrixSize):
		up = np.where(columnsC == i)
		right = np.where(rowsR == i)
		#look up
		for j in up:
			# Commented part was written with element-not-used function in mind. 
			#if (len(up[0]) > 1) | (len(right[0]) > 1):
			#	nodesBank[i] = min(rowsC[j])+1 	#increment the node number to keep 0 for Spice GND
			#else:
			#	nodesBank[i] = 0	#Here 0 is used as a marker to put nothing there (elemnt pin not used).
			
			#Let's assume every given element must be used. (A node might get unused - dangling)
			#Yes, but though, both SpiceOpus and HSpice hangs in some occasions when there is a dangling node. Need to fix this. 
			nodesBank[i] = min(rowsC[j])+1 	#increment the node number to keep 0 for Spice GND
		if i > (BigMatrixSize-1-NofOutConns):
			OutConnsBank[i] = nodesBank[i]	
	print nodesBank
	#print "OutConnsBank: ",OutConnsBank
	
	#another solution to gain nodesBank without need for fullRedundancyBigCircuitMatrix  #
	#WARNING Does not work as predicted                                                  #
	rowsR,columnsR,columnsC,rowsC = sortedNonZeroIndices(BigCircuitMatrix)               #
	nodesBank2 = np.zeros(shape=(BigMatrixSize),dtype=int)                               #
	for i in range(0,BigMatrixSize):                                                     #
	  up = np.where(columnsC == i)[0]                                                    #
	  print up                                                                           #
	  if len(up)==1:                                                                     #
	    nodesBank2[i] = rowsC[np.min(up)]                                                #
	  elif len(up)>1:                                                                    #
	    while len(up)>1:                                                                 #
	      #hodi gor, levo,                                                               #
	      up = np.where(columnsC == rowsC[np.min(up)])[0]                                #
	      print "\t",up                                                                  #
	    nodesBank2[i] = rowsC[np.min(up)]                                                #
	nodesBank2 = nodesBank2 + 1                                                          #
	print nodesBank2                                                                     #
	
	#test
	if not all(nodesBank == nodesBank2):
	  print nodesBank == nodesBank2
	  raw_input()
	
	
	#Open empty file for netlist
	ID = "g_" + str(generationNum) + "_i_" + str(individualNum)
	netlistName = ID + "_subckt.cir"
	circ = open(ID + "_subckt.cir", "w")
	circ.write("*%s \n" %netlistName)

	circ.write(".SUBCKT %s %d " %(HOTCIRC, OutConnsBank[BigMatrixSize-NofOutConns]))	#GND
	circ.write("%d " %(OutConnsBank[BigMatrixSize-NofOutConns+1]))				#Vsupp
	#circ.write("%d " %(OutConnsBank[BigMatrixSize-NofOutConns+2]))				#Iref
	circ.write("%d " %(OutConnsBank[BigMatrixSize-NofOutConns+2]))				#Vout
	circ.write("vsupp ")									#vsupp and vsupn are global nodes - so no need here 
	circ.write("vsupn ")
	circ.write("\n")
	#circ.write("vdd \n")
	#circ.write("%d " %(OutConnsBank[BigMatrixSize-NofOutConns+3]))
	#circ.write("%d " %(OutConnsBank[BigMatrixSize-NofOutConns+4]))
	#circ.write("%d " %(OutConnsBank[BigMatrixSize-NofOutConns+5]))
	#circ.write("%d \n" %(OutConnsBank[BigMatrixSize-NofOutConns+6]))

	#---------------------------	
	#---------------------------	

	value_index = 0		#goes with ValueVector which contains element values
	multipl_index = 0	#goes with Multipl, which contains multiplication factors for mos (only integer)
	#go thru all lines except connections (GND, Vp, Vn, Vin1, Vin2, Vout1, Vout2 - they shall not connect to each other!)
	i=0
	while i < BigMatrixSize:
		#Here, every used element is written into a netlist and given nodes and values/models
		#if-statements must strictly follow the sequence of elements in BigCircuitMatrix.
		if nodesBank[i] != 0:
			if i < np.sum(ALLPINS[0:1]): #NOTE: see globalVars.py for ALLPINS
				#put resistor------------------------------------------------------------------------
				if (nodesBank[i]==nodesBank[i+1]) | (sum(nodesBank == nodesBank[i]) < 2):
				  circ.write("*") #exclude element from netlist if element pins connected together or if dangeling node detected
				circ.write("R_%s %s " %(i, nodesBank[i]))	#put name of element and + node
				circ.write("%s %s \n" %(nodesBank[i+1], ValueVector[value_index]))	#put - node and value
				value_index = value_index + 1									#increment value_index
				i = i+2
				#-------------------------------------------------------------------------------------
			elif i < np.sum(ALLPINS[0:2]):
				#put capacitor------------------------------------------------------------------------
				if (nodesBank[i]==nodesBank[i+1]) | (sum(nodesBank == nodesBank[i]) < 2):
				  circ.write("*") #exclude element from netlist
				circ.write("C_%s %s " %(i, nodesBank[i]))	#put name of element and + node
				circ.write("%s %s \n" %(nodesBank[i+1], ValueVector[value_index]))	#put - node and value
				value_index = value_index + 1									#increment value_index
				i = i+2
				#-------------------------------------------------------------------------------------
			elif i < np.sum(ALLPINS[0:3]):
				#put inductor------------------------------------------------------------------------
				if (nodesBank[i]==nodesBank[i+1]) | (sum(nodesBank == nodesBank[i]) < 2):
				  circ.write("*") #exclude element from netlist
				circ.write("L_%s %s " %(i, nodesBank[i]))	#put name of element and + node
				circ.write("%s %s \n" %(nodesBank[i+1], ValueVector[value_index]))	#put - node and value
				value_index = value_index + 1									#increment value_index
				i = i+2
				#-------------------------------------------------------------------------------------
			elif i < np.sum(ALLPINS[0:4]):
				#put zener------------------------------------------------------------------------
				if (nodesBank[i]==nodesBank[i+1]) | (sum(nodesBank == nodesBank[i]) < 2):
				  circ.write("*") #exclude element from netlist
				circ.write("dz_%s %s " %(i, nodesBank[i]))	#put name of element and + node
				circ.write("%s %s \n" %(nodesBank[i+1], MODEL_ZENER))	#put - node and model name
				i = i+2
				#-------------------------------------------------------------------------------------
			elif i < np.sum(ALLPINS[0:5]):
				#put NPN------------------------------------------------------------------------
				if (len(set([nodesBank[i], nodesBank[i+1], nodesBank[i+2]]))==1) | (sum(nodesBank == nodesBank[i]) < 2):
				  circ.write("*") #exclude element from netlist
				circ.write("q_%s %s " %(i, nodesBank[i]))	#put name of element and COLLECTOR node
				circ.write("%s " %(nodesBank[i+1]))			#put BASE node
				circ.write("%s %s \n" %(nodesBank[i+2], MODEL_NPN))	#put EMITTER node and model name
				i = i+3
				#----------------------------------------------------------------------------	
			elif i < np.sum(ALLPINS[0:6]):
				#put 3 parallel PNPs (subckt)------------------------------------------------
				if (len(set([nodesBank[i], nodesBank[i+1], nodesBank[i+2]]))==1) | (sum(nodesBank == nodesBank[i]) < 2):
				  circ.write("*") #exclude element from netlist
				circ.write("x_%s %s " %(i, nodesBank[i]))	#put name of element and COLLECTOR node
				circ.write("%s " %(nodesBank[i+1]))			#put BASE node
				circ.write("%s %s \n" %(nodesBank[i+2], "par3pnp"))	#put EMITTER node and subckt name
				i = i+3
				#----------------------------------------------------------------------------	
			elif i < np.sum(ALLPINS[0:7]):
				#put PNP------------------------------------------------------------------------
				if (len(set([nodesBank[i], nodesBank[i+1], nodesBank[i+2]]))==1) | (sum(nodesBank == nodesBank[i]) < 2):
				  circ.write("*") #exclude element from netlist
				circ.write("q_%s %s " %(i, nodesBank[i]))	#put name of element and COLLECTOR node
				circ.write("%s " %(nodesBank[i+1]))			#put BASE node
				circ.write("%s %s \n" %(nodesBank[i+2], MODEL_PNP))	#put EMITTER node and model name
				i = i+3
				#----------------------------------------------------------------------------	
			
			elif i < np.sum(ALLPINS[0:8]):
				#put NMOS------------------------------------------------------------------------
				#if (len(set([nodesBank[i], nodesBank[i+1], nodesBank[i+2]]))==1) | (sum(nodesBank == nodesBank[i]) < 2):
				#  circ.write("*") #exclude element from netlist
				circ.write("xmn_%s %s " %(i, nodesBank[i]))	#put drain
				circ.write("%s " %(nodesBank[i+1]))		#put gate
				circ.write("%s %s %s " %(nodesBank[i+2], "0", MODEL_NMOS))	#put source, bulk(to vss), and model name
				circ.write("w=%s l=%s m=%s \n" %(ValueVector[value_index], ValueVector[value_index+1], Multipl[multipl_index]))
				value_index = value_index + 2
				i = i+3
				multipl_index = multipl_index + 1
			
			elif i < np.sum(ALLPINS[0:9]):
				#put PMOS------------------------------------------------------------------------
				#if (len(set([nodesBank[i], nodesBank[i+1], nodesBank[i+2]]))==1) | (sum(nodesBank == nodesBank[i]) < 2):
				#  circ.write("*") #exclude element from netlist
				circ.write("xmp_%s %s " %(i, nodesBank[i]))	#put drain
				circ.write("%s " %(nodesBank[i+1]))		#put gate
				circ.write("%s %s %s " %(nodesBank[i+2], "vdd", MODEL_PMOS))	#put source, bulk(to vdd), and model name
				circ.write("w=%s l=%s m=%s \n" %(ValueVector[value_index], ValueVector[value_index+1], Multipl[multipl_index]))
				value_index = value_index + 2
				i = i+3
				multipl_index = multipl_index + 1
			
			elif i < np.sum(ALLPINS[0:10]):
				#put OPAMP1------------------------------------------------------------------------
				if (len(set([nodesBank[i], nodesBank[i+1], nodesBank[i+2]]))==1) | (sum(nodesBank == nodesBank[i]) < 2):
				  circ.write("*") #exclude element from netlist
				circ.write("x_%s %s " %(i, nodesBank[i]))	#put name of element and NON-INVERTING node 
				circ.write("%s " %(nodesBank[i+1]))		#put INVERTING node
				circ.write("vsupp ")				#put POSITIVE POWER SUPPLY node NOTE: Change according to topdc.cir circuit
				circ.write("vsupn ")				#put NEGATIVE POWER SUPPLY node NOTE: Change according to topdc.cir circuit
				circ.write("%s %s \n" %(nodesBank[i+2], MODEL_OPAMP1))	#put OUTPUT node and model name
				i = i+3
			
			#Writing various subcircuits to netlist:
			elif i < np.sum(ALLPINS[0:11]):#UNDER CONSTRUCTION from here on
				#put PMosCurrSrc1Stg------------------------------------------------------------------
				if (len(set([nodesBank[i], nodesBank[i+1], nodesBank[i+2]]))==1) | (sum(nodesBank == nodesBank[i]) < 2):
				  circ.write("*") #exclude element from netlist
				circ.write("x_%s %s " %(i, nodesBank[i]))	#put name of element and 1st node 
				circ.write("%s " %(nodesBank[i+1]))		#put 2nd node
				circ.write("%s " %(nodesBank[i+2]))		#put 3rd node
				circ.write("PMosCurrSrc1stg ")			#put SUBCIRCUIT name
				circ.write("w=%s l=%s \n" %(ValueVector[value_index], ValueVector[value_index+1]))
				value_index = value_index + 2
				i = i+3
			  
			elif i < np.sum(ALLPINS[0:12]):
				#put NofCascPMosCurrSrc1stg------------------------------------------------------------------
				if (len(set([nodesBank[i], nodesBank[i+1], nodesBank[i+2]]))==1) | (sum(nodesBank == nodesBank[i]) < 2):
				  circ.write("*") #exclude element from netlist
				circ.write("x_%s %s " %(i, nodesBank[i]))	#put name of element and 1st node 
				circ.write("%s " %(nodesBank[i+1]))		#put 2nd node
				circ.write("%s " %(nodesBank[i+2]))				#put 3rd node
				circ.write("CascPMosCurrSrc1stg ")		#put SUBCIRCUIT name
				circ.write("w=%s l=%s \n" %(ValueVector[value_index], ValueVector[value_index+1]))
				value_index = value_index + 2
				i = i+3
			  
			elif i < np.sum(ALLPINS[0:13]):
				#put NofCascPMosCurrSrc1stg------------------------------------------------------------------
				if (len(set([nodesBank[i], nodesBank[i+1], nodesBank[i+2]]))==1) | (sum(nodesBank == nodesBank[i]) < 2):
				  circ.write("*") #exclude element from netlist
				circ.write("x_%s %s " %(i, nodesBank[i]))	#put name of element and 1st node 
				circ.write("%s " %(nodesBank[i+1]))		#put 2nd node
				circ.write("%s " %(nodesBank[i+2]))				#put 3rd node
				circ.write("NMosAmp1ResOnSrc ")		#put SUBCIRCUIT name
				circ.write("w=%s l=%s r=%s \n" %(ValueVector[value_index], ValueVector[value_index+1], ValueVector[value_index+2]))
				value_index = value_index + 3
				i = i+3
			#From here on, put elements with 4 pins:
			elif i < np.sum(ALLPINS[0:14]):
				#put BJTNPNCurrSink------------------------------------------------------------------
				if (len(set([nodesBank[i], nodesBank[i+1], nodesBank[i+2], nodesBank[i+3]]))==1) | (sum(nodesBank == nodesBank[i]) < 2):
				  circ.write("*") #exclude element from netlist
				circ.write("x_%s %s " %(i, nodesBank[i]))	#put name of element and 1st node 
				circ.write("%s " %(nodesBank[i+1]))		#put 2nd node
				circ.write("%s " %(nodesBank[i+2]))		#put 3rd node
				circ.write("%s " %(nodesBank[i+3]))		#put 4rd node
				circ.write("BJTNPNCurrSink \n")			#put SUBCIRCUIT name
				i = i+4	  
			  
			elif i < np.sum(ALLPINS[0:15]):
				#put BJTNPNCurrSink------------------------------------------------------------------
				if (len(set([nodesBank[i], nodesBank[i+1], nodesBank[i+2], nodesBank[i+3]]))==1) | (sum(nodesBank == nodesBank[i]) < 2):
				  circ.write("*") #exclude element from netlist
				circ.write("x_%s %s " %(i, nodesBank[i]))	#put name of element and 1st node 
				circ.write("%s " %(nodesBank[i+1]))		#put 2nd node
				circ.write("%s " %(nodesBank[i+2]))		#put 3rd node
				circ.write("%s " %(nodesBank[i+3]))		#put 4rd node
				circ.write("BJTPNPCurrSrc \n")			#put SUBCIRCUIT name
				i = i+4	   
			  
			elif i < np.sum(ALLPINS[0:16]):
				#put NMosCurrMirr------------------------------------------------------------------
				if (len(set([nodesBank[i], nodesBank[i+1], nodesBank[i+2], nodesBank[i+3]]))==1) | (sum(nodesBank == nodesBank[i]) < 2):
				  circ.write("*") #exclude element from netlist
				circ.write("x_%s %s " %(i, nodesBank[i]))	#put name of element and 1st node 
				circ.write("%s " %(nodesBank[i+1]))		#put 2nd node
				circ.write("%s " %(nodesBank[i+2]))		#put 3rd node
				circ.write("%s " %(nodesBank[i+3]))		#put 4rd node
				circ.write("NMosCurrMirr ")		#put SUBCIRCUIT name
				circ.write("w=%s l=%s \n" %(ValueVector[value_index], ValueVector[value_index+1]))
				value_index = value_index + 2
				i = i+4
			  
			elif i < np.sum(ALLPINS[0:17]):
				#put CascNMosCurrMirr------------------------------------------------------------------
				if (len(set([nodesBank[i], nodesBank[i+1], nodesBank[i+2], nodesBank[i+3]]))==1) | (sum(nodesBank == nodesBank[i]) < 2):
				  circ.write("*") #exclude element from netlist
				circ.write("x_%s %s " %(i, nodesBank[i]))	#put name of element and 1st node 
				circ.write("%s " %(nodesBank[i+1]))		#put 2nd node
				circ.write("%s " %(nodesBank[i+2]))		#put 3rd node
				circ.write("%s " %(nodesBank[i+3]))		#put 4rd node
				circ.write("CascNMosCurrMirr ")		#put SUBCIRCUIT name
				circ.write("w=%s l=%s \n" %(ValueVector[value_index], ValueVector[value_index+1]))
				value_index = value_index + 2
				i = i+4	  
			  
			elif i < np.sum(ALLPINS[0:18]):
				#put PMosCurrSrc2Stg------------------------------------------------------------------
				if (len(set([nodesBank[i], nodesBank[i+1], nodesBank[i+2], nodesBank[i+3]]))==1) | (sum(nodesBank == nodesBank[i]) < 2):
				  circ.write("*") #exclude element from netlist
				circ.write("x_%s %s " %(i, nodesBank[i]))	#put name of element and 1st node 
				circ.write("%s " %(nodesBank[i+1]))		#put 2nd node
				circ.write("%s " %(nodesBank[i+2]))		#put 3rd node
				circ.write("%s " %(nodesBank[i+3]))		#put 4rd node
				circ.write("PMosCurrSrc2Stg ")		#put SUBCIRCUIT name
				circ.write("w1=%s l1=%s w2=%s l2=%s \n" %(ValueVector[value_index], ValueVector[value_index+1], ValueVector[value_index+2], ValueVector[value_index+3]))
				value_index = value_index + 4
				i = i+4	 	  
			  
			elif i < np.sum(ALLPINS[0:19]):
				#put CascPMosCurrSrc2Stg------------------------------------------------------------------
				if (len(set([nodesBank[i], nodesBank[i+1], nodesBank[i+2], nodesBank[i+3]]))==1) | (sum(nodesBank == nodesBank[i]) < 2):
				  circ.write("*") #exclude element from netlist
				circ.write("x_%s %s " %(i, nodesBank[i]))	#put name of element and 1st node 
				circ.write("%s " %(nodesBank[i+1]))		#put 2nd node
				circ.write("%s " %(nodesBank[i+2]))		#put 3rd node
				circ.write("%s " %(nodesBank[i+3]))		#put 4rd node
				circ.write("CascPMosCurrSrc2Stg ")		#put SUBCIRCUIT name
				circ.write("w1=%s l1=%s w2=%s l2=%s \n" %(ValueVector[value_index], ValueVector[value_index+1], ValueVector[value_index+2], ValueVector[value_index+3]))
				value_index = value_index + 4
				i = i+4	 	  
			  
			elif i < np.sum(ALLPINS[0:20]):
				#put PMosCascode------------------------------------------------------------------
				if (len(set([nodesBank[i], nodesBank[i+1], nodesBank[i+2], nodesBank[i+3]]))==1) | (sum(nodesBank == nodesBank[i]) < 2):
				  circ.write("*") #exclude element from netlist
				circ.write("x_%s %s " %(i, nodesBank[i]))	#put name of element and 1st node 
				circ.write("%s " %(nodesBank[i+1]))		#put 2nd node
				circ.write("%s " %(nodesBank[i+2]))		#put 3rd node
				circ.write("%s " %(nodesBank[i+3]))		#put 4rd node
				circ.write("PMosCascode ")			#put SUBCIRCUIT name
				circ.write("w1=%s l1=%s w2=%s l2=%s \n" %(ValueVector[value_index], ValueVector[value_index+1], ValueVector[value_index+2], ValueVector[value_index+3]))
				value_index = value_index + 4
				i = i+4	 	  
			  
			elif i < np.sum(ALLPINS[0:21]):
				#put NMosCascode------------------------------------------------------------------
				if (len(set([nodesBank[i], nodesBank[i+1], nodesBank[i+2], nodesBank[i+3]]))==1) | (sum(nodesBank == nodesBank[i]) < 2):
				  circ.write("*") #exclude element from netlist
				circ.write("x_%s %s " %(i, nodesBank[i]))	#put name of element and 1st node 
				circ.write("%s " %(nodesBank[i+1]))		#put 2nd node
				circ.write("%s " %(nodesBank[i+2]))		#put 3rd node
				circ.write("%s " %(nodesBank[i+3]))		#put 4rd node
				circ.write("NMosCascode ")		#put SUBCIRCUIT name
				circ.write("w1=%s l1=%s w2=%s l2=%s \n" %(ValueVector[value_index], ValueVector[value_index+1], ValueVector[value_index+2], ValueVector[value_index+3]))
				value_index = value_index + 4
				i = i+4	 
			  
			else:
				i = i+1
		else:
			i = i+1
	#convergence aid
	#add 1 GOhm between every node and GND to help HSpice converge if dangling nodes in circuit. Added 13.5.2016. Looks it works fine.
	circ.write("\n*Convergence-aid resistors:\n")
	for i in set(nodesBank):  
	  circ.write("R_%s %s " %(np.max(nodesBank)+1+i, i))	#put name of element and + node
	  circ.write("%s %s \n" %("0", "1e9"))	#put - node and value	
	  
	circ.write(".ends\n")
	circ.close()
	#print "**********************************************************************"
	#print "* Generation %s, individual %s's netlist created. *" %(generationNum, individualNum)
	#print "**********************************************************************"
	

#BW = 1900#1395.200
#CUTOFF = 2000#1855.83
#PARTPASS = 1.0 #(gain)
#PARTSTOP = 0.0 #(gain)



def geneticOperation(parent1, parent2, generationNum):
  """Takes two parnets of class 'circuit', choses mating or mutation based on given probabilities. 
  Returns mini population with both parents and two newborn babies."""
  family = population()

  #add parents to generation - copy all attributes and results
  family.add_individual(deepcopy(parent1))
  family.add_individual(deepcopy(parent2))

  if random.uniform(0,1) < matingProb:
    if random.uniform(0,1) < topologyGenOperProb:
      #sliceCrossover
      #child1, child2 = sliceCrossover(parent1, parent2)
      #family.add_individual(deepcopy(child1))
      #family.add_individual(deepcopy(child2))
      #lineExchangeCrossover
      #child1, child2 = lineExchangeCrossover(parent1, parent2)
      child1, child2 = pinViseCrossover2(parent1, parent2)
      family.add_individual(deepcopy(child1))
      family.add_individual(deepcopy(child2))
    else:
      ##rangeValueCrossover - (interpolation)
      #child1, child2 = valueReproduction(parent1, parent2,crossover="rangeCrossover")
      #family.add_individual(deepcopy(child1))
      #family.add_individual(deepcopy(child2))
      #sliceValueCrossover
      child1, child2 = valueReproduction(parent1, parent2,crossover="sliceCrossover")
      family.add_individual(deepcopy(child1))
      family.add_individual(deepcopy(child2))
  else:
    if random.uniform(0,1) < topologyGenOperProb:
      choice = random.randint(0,2)
      if choice == 0:
	mut1 = mutationADDnode(parent1)
	mut2 = mutationADDnode(parent2)
	family.add_individual(deepcopy(mut1))
	family.add_individual(deepcopy(mut2))
      elif choice == 1:
	mut1 = mutationREMOVEnode(parent1)
	mut2 = mutationREMOVEnode(parent2)
	family.add_individual(deepcopy(mut1))
	family.add_individual(deepcopy(mut2))
      elif choice == 2:
	mut1 = mutationREMOVEnode(parent1)
	mut1 = mutationADDnode(mut1)
	mut2 = mutationREMOVEnode(parent2)
	mut2 = mutationADDnode(mut2)		    
	family.add_individual(deepcopy(mut1))
	family.add_individual(deepcopy(mut2))
    else:
      #ValueVector mutation
      mut1, mut2 = valueReproduction(parent1, parent2, crossover="mutation")
      family.add_individual(deepcopy(mut1))
      family.add_individual(deepcopy(mut2))
  
  return family


def geneticOperation3(parent1, parent2, generationNum):
  """Takes two parnets of class 'circuit', choses mating or mutation based on given probabilities. 
  Returns only array with two children.
  While every child differs from every parent, this procedure is repeating. This is how we do not get duplicates."""
  debug = 0
  duplicantsFound = True
  while duplicantsFound:
    children = []
    matingRND = random.uniform(0,1)
    topologyGenOperRND = random.uniform(0,1)
    mutationRND = random.randint(0,3)
    if debug: print "Generated random numbers for tournamnet:", matingRND, topologyGenOperRND, mutationRND
  
    if matingRND < matingProb:
      if topologyGenOperRND < topologyGenOperProb:
	child1, child2 = pinViseCrossover2(parent1, parent2)
	children.append(copy(child1))#deep
	children.append(copy(child2))#deep
	if debug: print "mtrx pinVise"
      else:
	#sliceValueCrossover
	child1, child2 = valueReproduction(parent1, parent2,crossover="sliceCrossover")
	children.append(copy(child1))#deep
	children.append(copy(child2))#deep
	if debug: print "valu sliceValueXO"
    else:
      if topologyGenOperRND < topologyGenOperProb:
	choice = mutationRND
	if choice == 0:
	  mut1 = mutationADDnode(parent1)
	  mut2 = mutationADDnode(parent2)
	  children.append(copy(mut1))#deep
	  children.append(copy(mut2))#deep
	  if debug: print "mtrx ADD node"
	elif choice == 1:
	  mut1 = mutationREMOVEnode(parent1)
	  mut2 = mutationREMOVEnode(parent2)
	  children.append(copy(mut1))#deep
	  children.append(copy(mut2))#deep
	  if debug: print "mtrx REMOVE node"
	elif choice == 2:
	  mut1 = mutationREMOVEnode(parent1)
	  mut1 = mutationADDnode(mut1)
	  mut2 = mutationREMOVEnode(parent2)
	  mut2 = mutationADDnode(mut2)
	  if debug: print "mtrx MOVE node"
	  children.append(copy(mut1))#deep
	  children.append(copy(mut2))#deep
	elif choice == 3:
	  children.append(copy(createRandomBigCircuitMatrix(copy(createRandomValueVector()))))#deep
	  children.append(copy(createRandomBigCircuitMatrix(copy(createRandomValueVector()))))#deep
      else:
	#ValueVector mutation
	mut1, mut2 = valueReproduction(parent1, parent2, crossover="mutation")
	children.append(copy(mut1))#deep
	children.append(copy(mut2))#deep
	if debug: print "valu mutate"
	
    #duplicantsFound = False
    #print "reduMtrxHash", [parent1.reduMtrxHash, parent2.reduMtrxHash, children[0].reduMtrxHash, children[1].reduMtrxHash]
    #print "fullHash", [parent1.fullHash, parent2.fullHash, children[0].fullHash, children[1].fullHash]
    #print set([parent1.fullHash, parent2.fullHash, children[0].fullHash, children[1].fullHash])
    if len(set([parent1.fullHash, parent2.fullHash, children[0].fullHash, children[1].fullHash])) == 4:
      duplicantsFound = False
  
  return children

  
  
def tournament(generation, NofElite, tournamentSize, POP_SIZE, sortedPool_Indices):
  """Takes whole generation and returns mating pool and sorted indices by score of each individual.
  Procedure copies NofElite of best individuals from given generation.
  Tournament size is adjustable by parameter.
  Each tournament, one of best among random chosen proceeds to mating pool. """
  gen_Pool = copy(generation.pool)
  gen_Scores = copy(generation.scores)
  gen_matrixDensities = copy(generation.matrixDensities)
  gen_matrixQuaziIDs = copy(generation.matrixQuaziIDs)
  
  tour_Pool = []
  tour_Scores = []
  #tour_matrixDensities = []
  #tour_matrixQuaziIDs = []
  
  #copy elite straightforward
  tour_Pool = np.append(tour_Pool, generation.pool[sortedPool_Indices[:NofElite]])
  tour_Scores = np.append(tour_Scores, generation.scores[sortedPool_Indices[:NofElite]])
  #tour_matrixDensities = np.append(tour_matrixDensities, generation.matrixDensities[sortedPool_Indices[:NofElite]])
  #tour_matrixQuaziIDs = np.append(tour_matrixQuaziIDs, generation.matrixQuaziIDs[sortedPool_Indices[:NofElite]])
  
  for i in range(0,POP_SIZE/2-NofElite):
    randomIndices = []
    scoresTemp = []
    for j in range(0, tournamentSize):
     randomIndices.append(random.randint(NofElite,POP_SIZE-1))
    for j in range(0, tournamentSize):
     scoresTemp.append(gen_Scores[randomIndices[j]])
    tourIndices = np.argsort(scoresTemp, kind='mergesort')
    
    # here, best one from tournament is chosen and put to tour_Pool
    tour_Pool = np.append(tour_Pool, gen_Pool[randomIndices[tourIndices[0]]])
    tour_Scores = np.append(tour_Scores, gen_Scores[randomIndices[tourIndices[0]]])
    #tour_matrixDensities = np.append(tour_matrixDensities, gen_matrixDensities[randomIndices[tourIndices[0]]])
    #tour_matrixQuaziIDs = np.append(tour_matrixQuaziIDs, gen_matrixQuaziIDs[randomIndices[tourIndices[0]]]) 

  matingPool = population()
  matingPool.pool = np.append(matingPool.pool, tour_Pool)
  matingPool.scores = np.append(matingPool.scores, tour_Scores)
  #matingPool.matrixDensities = np.append(matingPool.matrixDensities, tour_matrixDensities)
  #matingPool.matrixQuaziIDs = np.append(matingPool.matrixQuaziIDs, tour_matrixQuaziIDs)
  
  mP_SortedIndices = np.argsort(matingPool.scores)
  
  return matingPool, mP_SortedIndices

  
  
  
def dispatchRandomCircuitObjectGeneration(i):
  """This function was created in order to dispatch a time consuming fullRedundancyBigCircuitMatrix creation.
  This function is particularly for random circuits."""
  a = createRandomBigCircuitMatrix(copy(createRandomValueVector()))
  #a.fullRedundancyMatrix = fullRedundancyBigCircuitMatrix(a.BigCircuitMatrix)
  a.fullRedundancyMatrix = a.fullRedundancyMatrix
  return a

def dispatchCircuitObjectGeneration(somecircuit, i):
  """This function was created in order to dispatch a time consuming fullRedundancyBigCircuitMatrix creation."""
  #somecircuit.fullRedundancyMatrix = fullRedundancyBigCircuitMatrix(somecircuit.BigCircuitMatrix) #TEST 14.6.2017
  somecircuit.fullRedundancyMatrix = somecircuit.fullRedundancyMatrix
  return somecircuit



  