#Main script
#Evolutionary algorythm for automatic topology synthesis
#2015-2016 Ziga Rojec - EDA, FE, UL 

import numpy as np
import random
from time import time, strftime
from copy import copy, deepcopy
import collections, pickle, shutil
import os, sys
from os import mkdir
import resource
from threading import Thread

#Paralellism modules
from pyopus.parallel import base
from pyopus.parallel.cooperative import cOS
from pyopus.parallel.mpi import MPI

#Optimisation modules
from pyopus.optimizer.psade import ParallelSADE
from pyopus.optimizer.base import Reporter, CostCollector, RandomDelay

#My modules
from globalVars import *
import AdamAndEve as AE #old adam style...
import adam_dict #The new Adam
from reproduction import *
from scoreFunctions import *
from paramOptimizer import *
from buildingBlocksBank import *

#Other settings
np.set_printoptions(threshold=sys.maxsize, linewidth=1000)	#print whole matrix when finished
#Settings
random.seed(seedN)	#Fixed seed
np.random.seed(seedN)

problems = {'scoreCirc_CmosVoltageReference_2':scoreCirc_CmosVoltageReference_2,
	    'scoreCirc_ActiveFilter_2':scoreCirc_ActiveFilter_2,
	    'scoreCirc_PassiveBandPass':scoreCirc_PassiveBandPass,
	    'scoreCirc_HighPass':scoreCirc_HighPass,
	    'scoreCirc_commonEmitterAmp_resilenceMode':scoreCirc_commonEmitterAmp_resilenceMode,
	    } #set this also in paramOptimizer

PROBLEM = problems[PROBLEMname]
MOEAMODE = 0 # DO NOT CHANGE

def check_input():
    print("Starting listener thread. Type STOP if you want to end the algorithm gently.")
    while True:
        _in = input()
        print("received input: " + _in)
        if _in.lower() == "stop":
            os.mknod("./STOP")  # Creates stop file. 
            break



if __name__=='__main__':
  #---------------------------------------------------
  
  t0 = time()
  print()
  print("+++	 MATRIX (R)EVOLUTIONS STARTED	+++\n")

  ###MINI EVOLUTION CYCLE	
  startdate = strftime("%Y_%m_%d")
  starttime = strftime("%H-%M")
  print("Starting date: ", startdate, ", starting time: ", starttime)
  
  # GET YOUR WORKING DIRECTORY!
  working_directory_path = os.getcwd()
  os.chdir("../_MAIN_data")

  datadirname = "data_" + startdate + "_" + starttime
  os.mkdir(datadirname)
  os.mkdir(datadirname + "/" + "diversityPlots")
  shutil.copy2(working_directory_path + '/globalVars.py', datadirname + '/globalVars.py')	#copy input scripts for logging reasons
  shutil.copy2(working_directory_path + '/buildingBlocksBank.py', datadirname + '/buildingBlocksBank.py')	
  shutil.copy2(working_directory_path + '/scoreFunctions.py', datadirname + '/scoreFunctions.py')  
  output = open("data.pkl", "wb")
  os.chdir(working_directory_path)  


  # Set up MPI for parallel computing
  cOS.setVM(MPI(mirrorMap={
      #TODO set models in home folder for MPI. Fix that...
    #'models_for_start.inc':'.', 
    #'topdc_robust_commonemitter.cir':'.', 
  }))
  
  generationNum = 0
  generations = []
  
  hotGen = population() #hot generation (current)
  oldGen = population() #previous generation (old)
  
  
  currentBestScore = np.Inf
  bestScoresList = []
  averageScoresList = []
  
  
  if os.path.exists("./STOP"):
      input("Warning: A STOP file exists in this directory. Evolution will only last for 1 gen. You can remove the STOP file now. Press key when ready.")
  
  if continuee:
    with open("backdata.pkl","rb") as pkl_file: # Changed to rb to load pickle correctly
      data = pickle.load(pkl_file)
    pkl_file.close()
    
    print("Ressurecting old population...")
    generation = data[0]
    generationNum = 0#data[1]
    bestScoresList = data[2]
    result = data[3]
    bestI = data[4]
    #datadirname = data[5]
    BigMatrixSize = data[6]
    POP_SIZE = data[7]
    averageScoresList = data[8]
    #generations.append(generation)
    hotGen = deepcopy(generation)
    
    input("...old population resurrected.")
  else:
    print("Creating initial population...")
    
    NEWindividuals = cOS.dispatch(jobList=((dispatchRandomCircuitObjectGeneration, [i]) for i in range(0,POP_SIZE)), remote=True)
    for i in range(0,POP_SIZE):
        NEWindividuals[i].generationNum = generationNum
        NEWindividuals[i].individualNum = i
        hotGen.add_individual(NEWindividuals[i])
    
    if insertAdam:
      hotGen.pool[0] = AE.adam 	#Insert an already designed circuit to optimise.
      hotGen.pool[0].fullRedundancyMatrix = fullRedundancyBigCircuitMatrix(AE.adam.BigCircuitMatrix)
      print("Adam-circuit inserted into population.")
      # Ubder construction. Read from adam_dict netlist!
    input("...initial population created.")
    
  
  #---CREATE INITIAL POPULATION---#
  #---EVALUATE & SORT INITIAL POPULATION---# 
  print("EVALUATING GENERATION %d" %generationNum)
  stw0 = time()
  
  listener_thread = Thread(target = check_input)  # Starting the STOP signal listener. 
  listener_thread.start()
  
  #Parallel!! :)
  results=cOS.dispatch(jobList=((PROBLEM, [hotGen.pool[i], generationNum, i, False]) for i in range(0,len(hotGen.pool))), remote=True)
  results = np.array(results)
  
  stw1 = time()
  print("Evaluation of initial population lasted for %f s" %(stw1-stw0))
  #input("Tisni!")

  hotGen.scores = np.transpose(results[:,0])
  #hotGen.matrixDensities = np.transpose(results[:,1])
  #hotGen.matrixQuaziIDs = np.transpose(results[:,2])
  
  #Sort population
  sortedPool_Indices = np.argsort(hotGen.scores, kind='mergesort')
  currentBestScore = hotGen.scores[sortedPool_Indices[0]]
  
  #najboljsi v generaciji je...
  print(":::GENERATION %04.d - BEST ONE::: %f ::YEAH!::" %(generationNum,currentBestScore))
  
  printer(results[sortedPool_Indices[0]], stw0, generationNum, problem=PROBLEMname)
  bestScoresList.append(hotGen.scores[sortedPool_Indices[0]])
  averageScoresList.append(np.average(hotGen.scores))
  
  DONE = 0
  while DONE == 0:
    generationNum = generationNum + 1	#risen the generations counter
    oldGen = deepcopy(hotGen)
    hotGen = population() #hot generation (current) initialization
    
    #generations.append(population(generationNum))		#append new generation instance to a list
    
    #set up mating pool - tournament
    matingPool, mP_SortedIndices = tournament(oldGen, NofElite, tournamentSize, POP_SIZE, sortedPool_Indices)
    
    tempPool = population()
    for i in range(0, len(mP_SortedIndices)//2):
      parent1 = matingPool.pool[i]
      parent2 = matingPool.pool[random.randint(len(mP_SortedIndices)//2,len(mP_SortedIndices)-1)]
      family = geneticOperation(parent1, parent2, generationNum)
      for j in range(0, len(family.pool)):
        tempPool.add_individual(family.pool[j])
    
    #dispach creating of full redundance in circuit matrices
    individuals = cOS.dispatch(jobList=((dispatchCircuitObjectGeneration, [tempPool.pool[i], i]) for i in range(0,len(tempPool.pool))), remote=True)
    for i in range(0,len(tempPool.pool)):
        # WATCH Here add gennum ind num to individual object
        individuals[i].generationNum = generationNum
        individuals[i].individualNum = i
        tempPool.pool[i] = individuals[i]
   
    #remove duplicated individuals created during evolution----------------
    #Warning. When optimizing ValueVector, duplicates are not among circuits with same topology!
    tempPool_hashes = []	##every tempPool individual gets a hash
    
    for i in range(0, len(tempPool.pool)):
      #tempPool_hashes.append(hash(tempPool.pool[i].BigCircuitMatrix.tostring()) + hash(tempPool.pool[i].ValueVector.tostring())) #TEST: getting hash of FULL MATRIX
      tempPool_hashes.append(hash(tempPool.pool[i].fullRedundancyMatrix.tostring()))# + hash(tempPool.pool[i].ValueVector.tostring())) #TEST: getting hash of FULL MATRIX

    dups = collections.defaultdict(list)
    for i, item in enumerate(tempPool_hashes):
      dups[item].append(i)
    
    uniqateI = []
    for i in set(tempPool_hashes):
      uniqateI.append(dups[i][0])
    
    temptempPool = deepcopy(tempPool.pool)
    tempPool.pool = []
    
    for i in uniqateI:
      tempPool.add_individual(temptempPool[i])
    #----------------------------------------------------------------------
    
    #Evaluate the big temporary pool, sort results and append to current generation
    results=cOS.dispatch(jobList=((PROBLEM, [tempPool.pool[i], generationNum, i, False]) for i in range(0,len(tempPool.pool))), remote=True)
    results = np.array(results)
	
    tempPool.scores = np.append(tempPool.scores,np.transpose(results[:,0]))
    #tempPool.matrixDensities = np.append(tempPool.matrixDensities,np.transpose(results[:,1]))
    #tempPool.matrixQuaziIDs = np.append(tempPool.matrixQuaziIDs , np.transpose(results[:,2]))
    #Sort population
    sortedTempPool_Indices = np.argsort(tempPool.scores, kind='mergesort')

    #-------------#      
    #--OPTIMISE SOME OF BEST AND TAKE THEM WITH YOU--#
    
    bSL_npA = np.array(bestScoresList)
    deltaScore = abs(np.average(bSL_npA[-10:])-tempPool.scores[sortedTempPool_Indices[0]])
    
    #if (optimise == True) & (currentBestScore < 300) & ((currentBestScore < 100) | (deltaScore < 1.0)) & (not generationNum%20):
    if (optimise == True) & (((generationNum > 30) & (deltaScore < 1.0) & (not generationNum%20)) | (generationNum < 2)):
      bestToOptimise = [0, random.randint(1,NofElite), random.randint(NofElite+1, 10*NofElite)]
      if generationNum < 2:
        bestToOptimise = [0] #in first generation optimise just adam!
     
      for i in bestToOptimise:
        topology = copy(tempPool.pool[sortedTempPool_Indices[i]].BigCircuitMatrix)
        values = copy(tempPool.pool[sortedTempPool_Indices[i]].ValueVector)

        print("Circuit ", i, "from gen", generationNum, "...")
        maxiter = 3000 if generationNum < 100 else 8000
        x, f = optimiseCircuit(topology, values, maxiter)
        tempPool.pool = np.append(deepcopy(circuit(topology, x)), tempPool.pool)
        tempPool.pool[0] = dispatchCircuitObjectGeneration(tempPool.pool[0],0)
        
        tempPool.scores = np.append(copy(np.float64(f)),tempPool.scores)
        #tempPool.matrixDensities = np.append(copy(tempPool.matrixDensities[sortedTempPool_Indices[i]]), tempPool.matrixDensities)
        #tempPool.matrixQuaziIDs = np.append(copy(tempPool.matrixQuaziIDs[sortedTempPool_Indices[i]]), tempPool.matrixQuaziIDs)
        
        #sort them again
        sortedTempPool_Indices = np.argsort(tempPool.scores, kind='mergesort')
        sortedTempPool_Indices = list(sortedTempPool_Indices)
    
    #-------------#  


    #adjust NofRANDOMS
    if len(tempPool.pool) < (POP_SIZE - NofRANDOMS):
      NofRANDOMS = POP_SIZE - len(tempPool.pool)    

    #adding random new individuals"
    indN = len(hotGen.pool)
    
    NEWindividuals = cOS.dispatch(jobList=((dispatchRandomCircuitObjectGeneration, [i]) for i in range(0,NofRANDOMS)), remote=True)
    for i in range(0,NofRANDOMS):
      hotGen.pool = np.append(hotGen.pool, deepcopy(NEWindividuals[i]))
      
    #for i in range(0, NofRANDOMS):
    #  mutant = createRandomBigCircuitMatrix(copy(createRandomValueVector()))
    #  hotGen.pool = np.append(hotGen.pool, deepcopy(mutant))

    #evaluate randoms
    results = cOS.dispatch(jobList=((PROBLEM, [hotGen.pool[i], generationNum, i, False]) for i in range(0, NofRANDOMS)), remote=True)
    #sort results
    results = np.array(results)
    hotGen.scores = np.append(hotGen.scores, copy(np.transpose(results[:,0])))
    #hotGen.matrixDensities = np.append(hotGen.matrixDensities, copy(np.transpose(results[:,1])))
   # hotGen.matrixQuaziIDs = np.append(hotGen.matrixQuaziIDs, copy(np.transpose(results[:,2])))

    #adding randoms and offspring together
    hotGen.pool = np.append(hotGen.pool, tempPool.pool[sortedTempPool_Indices[:(POP_SIZE-NofRANDOMS)]])
    hotGen.scores = np.append(hotGen.scores, tempPool.scores[sortedTempPool_Indices[:(POP_SIZE-NofRANDOMS)]])
    #hotGen.matrixDensities = np.append(hotGen.matrixDensities, tempPool.matrixDensities[sortedTempPool_Indices[:(POP_SIZE-NofRANDOMS)]])
    #hotGen.matrixQuaziIDs = np.append(hotGen.matrixQuaziIDs, tempPool.matrixQuaziIDs[sortedTempPool_Indices[:(POP_SIZE-NofRANDOMS)]])
    #sort all
    sortedPool_Indices = np.argsort(hotGen.scores, kind='mergesort')
    currentBestScore = float(hotGen.scores[sortedPool_Indices[0]])
    bestScoresList.append(currentBestScore)
    averageScoresList.append(np.average(np.sort(hotGen.scores)[:(POP_SIZE-NofRANDOMS)]))
    
    
    #if int(currentBestScore) > int(min(bestScoresList)): #Jao dzubre kaj si blesav!
    #  print "Something went wrong BADLY."
    #  print bestScoresList
    # input()
    
    #---Saving results, picle, ...  and so on----#

    #Evaluate best one together with results for plotting
    debug += 3
    WinnerResults = cOS.dispatch(jobList=((PROBLEM, [hotGen.pool[sortedPool_Indices[0]], generationNum, 0, False]) for i in range(1,2)), remote=True)
    WinnerResults = np.array(WinnerResults)[0]
    
    printer(WinnerResults, stw0, generationNum, problem = PROBLEMname) # prints the score and results summary
    debug -= 3


    #Write winner netlist to a directiory of current run for inspection and manual simulation 
    os.chdir("../_MAIN_data")
    os.chdir("data_" + startdate + "_" + starttime)
    fullRedMx = fullRedundancyBigCircuitMatrix(deepcopy(hotGen.pool[sortedPool_Indices[0]].BigCircuitMatrix))
    makeNetlist_netlister(deepcopy(hotGen.pool[sortedPool_Indices[0]]))
    os.chdir("../") 

    #DUMP results for plotting
    data = [hotGen, generationNum, bestScoresList, WinnerResults, sortedPool_Indices[0], datadirname, BigMatrixSize, POP_SIZE, averageScoresList]
    #with open("data.pkl","wb") as output:
    with open(datadirname + "/backdata.pkl","wb") as output:
      pickle.dump(data, output)
    output.close()
    #shutil.copy2('data.pkl', datadirname + '/backdata.pkl')	#copy current .pkl data to current run folder
    shutil.copy2(datadirname + "/backdata.pkl", working_directory_path + "data.pkl")
    
    os.chdir(working_directory_path)
    #End of evolution? :
    if (generationNum > (endingGenNum-1)) or (currentBestScore < minimalScore) or os.path.exists("./STOP"):
      DONE = 1
      makeNetlist_netlister(hotGen.pool[sortedPool_Indices[0]])
      #print hotGen.pool[sortedPool_Indices[0]].BigCircuitMatrix
      

  
  listener_thread.join(0)
  cOS.finalize()
  print("\n+++	 MATRIX EVOLUTIONS ENDED	+++")
  
