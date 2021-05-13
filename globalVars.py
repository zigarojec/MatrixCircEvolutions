"""
globalVars.py
Ziga Rojec, EDA group, Faculty of Electrical Engineering, University of Ljubljana, Slovenia. Jan 2018.



In this script, the main parameters for the evolutionary algorithm are set, such as population size, random seed, various probabilities and so on. 

Here, the problem is also chosen (PROBLEMname). This is the cost function (high-level circuit definition) name from the script scoreFunctions.py. 
"""

import numpy as np
from buildingBlocksBank import buildBlocks, paramBounds

#Define the cost function--------------------------
PROBLEMname = 'scoreCirc_commonEmitterAmp_resilenceMode'
#---------------------------------------------------
#Problem names available:
# scoreCirc_ActiveFilter_2
# scoreCirc_CmosVoltageReference_2
# scoreCirc_PassiveBandPass
# scoreCirc_HighPass
# scoreCirc_commonEmitterAmp_resilenceMode

seedN = 0		#random seed
continuee = 	0	#to continue the algorithm from the last run
optimise = 	0	#turn on/off the global parameter optimizer PSADE
insertAdam = 	0	#insert a working circuit to further evolve
robustMode = 1      #iterate through all models listed for device when evaluating

			#(if continuee == 1) 
			#	copy backdata.pkl into main folder manually! 

#Debug
debug = 0

#Evolutionary algorythm parameters:-----------------
multipl = 	1
POP_SIZE = 	100*multipl
NofElite =  	2*multipl
NofRANDOMS = 	20*multipl
tournamentSize = 3

matingProb = 0.4
topologyGenOperProb = 1     # 0...only parameter optimization
                            # 1...only topology optimization

#---------------------------------------------------
endingGenNum = 4000 	#top number of generations
minimalScore = 0.01 	#stopping criteria based on Score (cost function)
#---------------------------------------------------
DONE = 0		#global variable that is set to 1 when an individual fulfills criteria
#---------------------------------------------------



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

