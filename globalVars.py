"""
globalVars.py
Ziga Rojec, EDA group, Faculty of Electrical Engineering, University of Ljubljana, Slovenia. Jan 2018.



In this script, the main parameters for the evolutionary algorithm are set, such as population size, random seed, various probabilities and so on. 

Here, the problem is also chosen (PROBLEMname). This is the cost function (high-level circuit definition) name from the script scoreFunctions.py. 
"""

import numpy as np
from buildingBlocksBank import buildBlocks, paramBounds

#Define the cost function--------------------------
PROBLEMname = 'scoreCirc_CmosVoltageReference_2'
#---------------------------------------------------
#Problem names available:
# scoreCirc_ActiveFilter_2
# scoreCirc_CmosVoltageReference_2
# scoreCirc_PassiveBandPass
# scoreCirc_HighPass

seedN = 0		#random seed
continuee = 	1	#to continue the algorithm from the last run
optimise = 	1	#turn on/off the global parameter optimizer PSADE
insertAdam = 	0	#insert a working circuit to further evolve

			#(if continuee == 1) 
			#	copy backdata.pkl into main folder manually! 

#Debug
debug = 0

#Evolutionary algorythm parameters:-----------------
multipl = 	4
POP_SIZE = 	100*multipl
NofElite =  	2*multipl
NofRANDOMS = 	20*multipl
tournamentSize = 3

matingProb = 0.4
topologyGenOperProb = 0.6 #<-- 1 for MIDEM article

#---------------------------------------------------
endingGenNum = 4000 	#top number of generations
minimalScore = 0.01 	#stopping criteria based on Score (cost function)
#---------------------------------------------------
DONE = 0		#global variable that is set to 1 when an individual fulfills criteria
#---------------------------------------------------

NofOutConns = 3		#number of outerConnections - such as Vin, Vout, GND. NOTE: When changing this, one has also zo update the makeNetlist method from reproduction.py . 

#TODO the code below to be moved to buildingBlocksBank.py maybe some time.

#Global values NOTE: strictly follow the sequence in buildBlocks if changing the software!

#Calculating various global variables connected to matrix representation of a circuit
ALLPINS = np.array([], dtype=int)	#number of pins for each device type in an array
for element in buildBlocks:
  ALLPINS = np.append(ALLPINS,  element['NofPins']*element['Quantity'])

BigMatrixSize = 0
for element in buildBlocks:
  BigMatrixSize = BigMatrixSize + element['NofPins']*element['Quantity']
BigMatrixSize = BigMatrixSize + NofOutConns

NofPARAMS = 0
for element in buildBlocks:
  NofPARAMS = NofPARAMS + len(element['ParamTypes'])*element['Quantity']
	  
Nof2poles = 0
Nof3poles = 0
Nof4poles = 0
for element in buildBlocks:
  if element['NofPins'] == 2:
    Nof2poles += element['Quantity']
  if element['NofPins'] == 3:
    Nof3poles += element['Quantity']
  if element['NofPins'] == 4:
    Nof4poles += element['Quantity']

#MODELS (Obsolete. Already set in buildingBlocksBank.)
MODEL_ZENER = 	"zd4v7"
MODEL_NPN = 	"bc238b"
SUBCKT_par3pnp = "par3pnp"
MODEL_PNP = 	"BC308B"
MODEL_OPAMP1 = 	"LM348T"#"LM741N"
MODEL_NMOS = 	"submodn"
MODEL_PMOS = 	"submodp"

#Mutiplication parameters for nmos and pmos 
#TODO: Move to buildingBlocksBank.
#Multipl = np.array([2, 2, 2, 2, 2, 2, 2, 2, 2, 2])
Multipl = np.array([1, 1, 1, 1, 1, 1, 1, 1, 1, 1])

