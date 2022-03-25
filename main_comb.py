#Main script
#Evolutionary algorythm for automatic topology synthesis
#2015-2016 Ziga Rojec - EDA, FE, UL 

import numpy as np
import random
from time import time, strftime
from copy import copy, deepcopy
import collections, shutil, pickle
import os, sys
from os import mkdir
from threading import Thread



#My modules
from globalVars import *
import AdamAndEve as AE #old adam style...
import adam_dict #The new Adam
from reproduction import *
from scoreFunctions import *
from paramOptimizer import *
from buildingBlocksBank import *

if globalVars.MOEA:
  from NSGAII import faster_nondominated_sort,crowding_distance_assignment,tournament_NSGA

from utils import dynamic_module_import

if __name__ == "__main__":
  module, PROBLEMCLASS = dynamic_module_import(PROBLEMpath + PROBLEMname, PROBLEMname)
  PROBLEM = getattr(PROBLEMCLASS, PROBLEMname)
  # module, PROBLEM = dynamic_module_import("scorefunctions.amplifiers.scoreCirc_commonEmitterAmp_resilenceMode", "scoreCirc_commonEmitterAmp_resilenceMode")



#Paralellism modules
from pyopus.parallel import base
from pyopus.parallel.cooperative import cOS
from pyopus.parallel.mpi import MPI

#Optimisation modules
from pyopus.optimizer.psade import ParallelSADE
from pyopus.optimizer.base import Reporter, CostCollector, RandomDelay


#Other settings
np.set_printoptions(threshold=sys.maxsize, linewidth=1000)	#print whole matrix when finished
#Settings
random.seed(seedN)	#Fixed seed
np.random.seed(seedN)

"""
problems = {'scoreCirc_CmosVoltageReference_2':scoreCirc_CmosVoltageReference_2,
	    'scoreCirc_ActiveFilter_2':scoreCirc_ActiveFilter_2,
	    'scoreCirc_PassiveBandPass':scoreCirc_PassiveBandPass,
	    'scoreCirc_HighPass':scoreCirc_HighPass,
	    'scoreCirc_commonEmitterAmp_resilenceMode':scoreCirc_commonEmitterAmp_resilenceMode,
	    } #set this also in paramOptimizer

PROBLEM = problems[PROBLEMname]
"""
MOEAMODE = globalVars.MOEA 



# Listener thread for intermediate STOPPING
def check_input():
    print("Starting listener thread. Type STOP if you want to end the algorithm gently.")
    while True:
        if os.path.exists("./STOP"):
          print("Exiting.")
          break
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
  starttime = strftime("%H-%M-%S")
  print("Starting date: ", startdate, ", starting time: ", starttime)
  
  # GET YOUR WORKING DIRECTORY!
  working_directory_path = os.getcwd()
  os.chdir("../_MAIN_data")

  datadirname = "data_" + startdate + "_" + starttime
  os.mkdir(datadirname)
  os.mkdir(datadirname + "/" + "diversityPlots")
  shutil.copy2(working_directory_path + '/globalVars.py', datadirname + '/globalVars.py')	#copy input scripts for logging reasons
  shutil.copy2(working_directory_path + '/buildingBlocksBank.py', datadirname + '/buildingBlocksBank.py')	
  #shutil.copy2(working_directory_path + '/scoreFunctions.py', datadirname + '/scoreFunctions.py')  
  shutil.copy2(working_directory_path + "/" + PROBLEMpath + PROBLEMname + ".py", datadirname + "/" + PROBLEMname + ".py")
  shutil.copy2(working_directory_path + "/" + PROBLEMpath + "/" + "runme.py", datadirname + "/" + "runme.py")


  output = open("data.pkl", "wb")
  os.chdir(working_directory_path)  


  # Set up MPI for parallel computing
  cOS.setVM(MPI(mirrorMap={
      #TODO set models in home folder for MPI. Fix that...
    'models_for_start.inc':'.',        # TODO TODO TODO remove this to a config file which is in gitignore!!
    'topdc_robust_commonemitter.cir':'.', 
    'topdc.cir':'.',
  }))
  
  generationNum = 0
  optimizedIndivids = [] 	#PSADE parameter-optimized individuals

  
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
  results=cOS.dispatch(jobList=((PROBLEM, [hotGen.pool[i], generationNum, i, globalVars.MOEA]) for i in range(0,len(hotGen.pool))), remote=True)
  results = np.array(results)
  
  if globalVars.MOEA:
    #Put together objects and objectivesScores
    for i in range(0,len(hotGen.pool)):
      hotGen.pool[i].objectivesScore = np.transpose(results[:,0])[i]


  stw1 = time()
  print("Evaluation of initial population lasted for %f s" %(stw1-stw0))

#Sort population
  if globalVars.MOEA:
    staSort = time()
    faster_nondominated_sort(hotGen)
    endSort = time()
    print("Fast_nondominated_sort finished in %.2f s." %(endSort-staSort))
    for i in range(0, len(hotGen.fronts)):
      crowding_distance_assignment(hotGen.fronts[i])
    print("Crowding_distance_assignment finished.")
  else:
    hotGen.scores = np.transpose(results[:,0])
    #hotGen.matrixDensities = np.transpose(results[:,1])
    #hotGen.matrixQuaziIDs = np.transpose(results[:,2])
    
    #Sort population
    sortedPool_Indices = np.argsort(hotGen.scores, kind='mergesort')
    currentBestScore = hotGen.scores[sortedPool_Indices[0]]
  
  #najboljsi v generaciji je...
  print(":::GENERATION %04.d - BEST ONE::: %f ::YEAH!::" %(generationNum,currentBestScore))
  
  if globalVars.MOEA:
    None
  else:
    printer(results[sortedPool_Indices[0]], stw0, generationNum, problem=PROBLEMname, resultspath = datadirname)
    bestScoresList.append(hotGen.scores[sortedPool_Indices[0]])
  averageScoresList.append(np.average(hotGen.scores))
  
  while globalVars.DONE == 0:
    generationNum += 1	#rise the generations counter

    oldGen = deepcopy(hotGen)
    hotGen = population() #hot generation (current) initialization
    
    if globalVars.MOEA:
      ##############
      #### MOEA ####
      aa = []
      bb = []
      tempPool = population()
      offspring = []
      for i in range(0, int(globalVars.POP_SIZE/2)):
        parent1 = tournament_NSGA(oldGen, globalVars.tournamentSize)
        #parent1 = deepcopy(oldGen.pool[np.random.randint(0,len(oldGen.pool))])
        parent2 = parent1
        while parent1.reduMtrxHash == parent2.reduMtrxHash:
          parent2 = tournament_NSGA(oldGen, globalVars.tournamentSize)
          #parent2 = deepcopy(oldGen.pool[np.random.randint(0,len(oldGen.pool))])
          #print "duplicant found 1x"
        bb.append(parent1)
        bb.append(parent2)
        twoChildren = geneticOperation3(parent1, parent2, generationNum)
        for j in range(0, len(twoChildren)):
          offspring.append(copy(twoChildren[j]))
    
      #if twoChildren[j].mtrxHash not in aa:
      #  aa.append(twoChildren[j].mtrxHash)
      #else:
      #  print("This topology is alredy in children.")

      #append all PSADE optimized to offspring
      if len(optimizedIndivids) > 0:
        for i in optimizedIndivids:
          offspring.append(i)


      #clean the offspring of duplicates immediately
      #	-cleaning by matrix - 		"mtrxHash"
      #	-cleaning by topology - 	"reduMtrxHash"
      #	-cleaning by circuit - 		"fullHash"
      print("Offspring len before cleaning:", len(offspring))
      offspring = removeDuplicatesFromArrayByAttribute(offspring, "fullHash")
      print("Offspring len after cleaning:", len(offspring))
      
      for i in range(0,len(offspring)):
        # WATCH Here add gennum ind num to individual object
        offspring[i].generationNum = generationNum
        offspring[i].individualNum = i


      #Evaluate the offspring - parents were already evaluated in ngen-1
      results=cOS.dispatch(jobList=((PROBLEM, [offspring[i], generationNum, i, True]) for i in range(0,len(offspring))), remote=True)
      results = np.array(results)
      
      #Put together objects and objectivesScores
      for i in range(0,len(offspring)):
        offspring[i].objectivesScore = np.transpose(results[:,0])[i]
          
      #Put together evaluated offspring and parents. After that, you perform non dominant sorting.
      for i in range(0, globalVars.POP_SIZE):
        tempPool.add_individual(copy(oldGen.pool[i]))
      for i in range(0, len(offspring)):
        tempPool.add_individual(copy(offspring[i]))
        #print(offspring[i].objectivesScore)
        #print(offspring[i].reduMtrxHash)
      if generationNum > 2:# and debug > 1:
        print("tempPool len before cleaning:", len(tempPool.pool))
        tempPool.pool = removeDuplicatesFromArrayByAttribute(tempPool.pool, "scoreHash")     
        print("tempPool len after cleaning:", len(tempPool.pool)) 
      
      
      #Sort population - fast-non-dominant sorting
      staSort = time()
      faster_nondominated_sort(tempPool)
      endSort = time()
      print("Fast_nondominated_sort finished in %.2fs." %(endSort-staSort))
      
      for i in range(0, len(tempPool.fronts)):
        #print "Front No:", i
        crowding_distance_assignment(tempPool.fronts[i])
      print("Crowding_distance_assignment finished.")
      #raw_input("Fast_nondominated_sort of gen"+str(generationNum)+" finished.")
      
      #Fill the new generation. First fronts first, better with greater crowding distance. 
      notFull = 1
      frontCount = 0
      while notFull:
        tempPool.fronts[frontCount].sort(key=lambda x: x.crow_dist, reverse=True)#sort the objects in fornt according to crowding distance
        for i in tempPool.fronts[frontCount]:
          hotGen.add_individual(copy(i))
          if len(hotGen.pool) == globalVars.POP_SIZE:
            notFull = 0
            break
        frontCount+=1
      #raw_input("Finished appending new population.")
      
      #for the case of diversity - checking hashes of all data
      a = []
      for i in tempPool.fronts[0]:
        a.append(np.array(i.objectivesScore).tostring())
      b = []
      for i in tempPool.fronts[0]:
        b.append(hash(i.fullRedundancyMatrix.tostring()))
      c = []
      for i in tempPool.fronts[0]:
        c.append(i.valuHash)
      d = []
      for i in tempPool.fronts[0]:
        d.append(i.fullHash)
      firstFlen = len(tempPool.fronts[0])
      #### END MOEA ####
      ##################

    else:
      ##################
      #### Basic EA ####

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
      #clean the offspring of duplicates immediately
      #	-cleaning by matrix - 		"mtrxHash"
      #	-cleaning by topology - 	"reduMtrxHash"
      #	-cleaning by circuit - 		"fullHash"
      print("TempPool len before cleaning:", len(tempPool))
      tempPool = removeDuplicatesFromArrayByAttribute(tempPool, "fullHash")
      print("TempPool len after cleaning:", len(tempPool))
      #----------------------------------------------------------------------
      
      #Evaluate the big temporary pool, sort results and append to current generation
      results=cOS.dispatch(jobList=((PROBLEM, [tempPool.pool[i], generationNum, i, False]) for i in range(0,len(tempPool.pool))), remote=True)
      results = np.array(results)
    
      tempPool.scores = np.append(tempPool.scores,np.transpose(results[:,0]))
      #tempPool.matrixDensities = np.append(tempPool.matrixDensities,np.transpose(results[:,1]))
      #tempPool.matrixQuaziIDs = np.append(tempPool.matrixQuaziIDs , np.transpose(results[:,2]))
      #Sort population
      sortedTempPool_Indices = np.argsort(tempPool.scores, kind='mergesort')

      #### END Basic EA ####
      ######################

    #Calculate time
    stw1 = time()		
    dtime = stw1-stw0
    m, s = divmod(dtime, 60)
    h, m = divmod(m, 60)

    if globalVars.MOEA:
      # Sort individulas according to distance from coordinate center
      objectivesList = []
      distxyz = []
      distxyzI = []
      for i in tempPool.pool:
        distxyz.append(np.sqrt(i.objectivesScore[0]**2 + i.objectivesScore[1]**2 + i.objectivesScore[2]**2))
      distxyzI = np.argsort(distxyz)
        
      currentBestScore = tempPool.pool[distxyzI[0]].objectivesScore
      bestScoresList.append(currentBestScore)

      #---OPTIMISE SOME OF BEST------------------------#
      if (globalVars.optimise == True) & (generationNum > 0) & ((not generationNum%10) or (generationNum == 1 & globalVars.insertAdam == 1)):
        bestToOptimise = [random.randint(0,globalVars.NofElite), random.randint(NofElite, 10*NofElite)]
    
        for i in bestToOptimise:
          topology = copy(tempPool.pool[distxyzI[i]].BigCircuitMatrix)
          values = copy(tempPool.pool[distxyzI[i]].ValueVector)
          print("PSADE on circuit ", i, "from gen", generationNum, "...")
          maxiter = 3000 if generationNum < 100 else 8000
          x, f = optimiseCircuit(topology, values, maxiter)
          optimizedIndivids.append(copy(circuit(topology, x)))
          
      #append those individuals to the NEXT generation!!
      else:
        optimizedIndivids = []
      #------------------------------------------------#  
    else: #not globalVars.MOEA    
      #--OPTIMISE SOME OF BEST AND TAKE THEM WITH YOU--#
      
      bSL_npA = np.array(bestScoresList)
      deltaScore = abs(np.average(bSL_npA[-10:])-tempPool.scores[sortedTempPool_Indices[0]])
      
      #if (optimise == True) & (currentBestScore < 300) & ((currentBestScore < 100) | (deltaScore < 1.0)) & (not generationNum%20):
      if (optimise == True) & (((generationNum > 30) & (deltaScore < 1e-1) & (not generationNum%20)) | (generationNum < 2)):
        bestToOptimise = [0, random.randint(1,NofElite), random.randint(NofElite+1, 10*NofElite)]
        if generationNum < 2:
          bestToOptimise = [0] #in first generation optimise just adam!
      
        for i in bestToOptimise:
          topology = copy(tempPool.pool[sortedTempPool_Indices[i]].BigCircuitMatrix)
          values = copy(tempPool.pool[sortedTempPool_Indices[i]].ValueVector)

          print("Circuit ", i, "from gen", generationNum, "...")
          maxiter = 300 if generationNum < 100 else 8000
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
        
      
    
    #---Saving results, picle, ...  and so on----#

    #Evaluate best one together with results for plotting
    if globalVars.MOEA:
      plotThisMany = 4
      WinnerResults = cOS.dispatch(jobList=((PROBLEM, [tempPool.pool[distxyzI[i]], generationNum, 0, True]) for i in range(0,plotThisMany)), remote=True)
      WinnerResults = np.array(WinnerResults)
    else: #not globalVars.MOEA
      WinnerResults = cOS.dispatch(jobList=((PROBLEM, [hotGen.pool[sortedPool_Indices[0]], generationNum, 0, False]) for i in range(1,2)), remote=True)
      WinnerResults = np.array(WinnerResults[0])
    

    printer(WinnerResults, stw0, generationNum, problem = PROBLEMname, resultspath = datadirname) # prints the score and results summary

    if globalVars.MOEA:
      print("\t - Unique of scores:",len(set(a)),"(/",firstFlen, "), of topologies:",len(set(b)),"(/",firstFlen, "), \n\t          of values:", len(set(c)),"(/",firstFlen, "), of individuals", len(set(d)),"(/",firstFlen, ") in 1st front.")
      print("\t- - - - - - - - - - - - - - - - - - - - - - - - - - ")

    #Write winner netlist to a directiory of current run for inspection and manual simulation 
    os.chdir("../_MAIN_data")
    os.chdir("data_" + startdate + "_" + starttime)
    #fullRedMx = fullRedundancyBigCircuitMatrix(deepcopy(hotGen.pool[sortedPool_Indices[0]].BigCircuitMatrix))
    if globalVars.MOEA:
      writeThisMany = 4
      for i in range(0, writeThisMany):
        makeNetlist_netlister(tempPool.pool[distxyzI[i]])

        allObjectiveValues = hotGen.pool[0].objectivesScore		#getting all objective values of hotGen for plotting
        for i in range(1, len(hotGen.pool)):
          allObjectiveValues = np.vstack((allObjectiveValues, hotGen.pool[i].objectivesScore))
        #DUMP results for plotting
        data = [tempPool, generationNum, bestScoresList, WinnerResults, None, datadirname, BigMatrixSize, globalVars.POP_SIZE, averageScoresList, allObjectiveValues]
    else:
      makeNetlist_netlister(deepcopy(hotGen.pool[sortedPool_Indices[0]]))

      #DUMP results for plotting
      data = [hotGen, generationNum, bestScoresList, WinnerResults, sortedPool_Indices[0], datadirname, BigMatrixSize, POP_SIZE, averageScoresList]

    os.chdir("../") 


    with open(datadirname + "/backdata.pkl","wb") as output:
      pickle.dump(data, output)
    output.close()
    #shutil.copy2('data.pkl', datadirname + '/backdata.pkl')	#copy current .pkl data to current run folder
    #shutil.copy2(datadirname + "/backdata.pkl", working_directory_path + "backdata.pkl")
    shutil.copy2(datadirname + "/backdata.pkl", "./backdata.pkl")
    
    os.chdir(working_directory_path)
    #End of evolution? :
    if globalVars.MOEA:
      stopCondition = (generationNum > (endingGenNum-1)) or (currentBestScore[0] < minimalScore) or os.path.exists("./STOP")
    else:
      stopCondition = (generationNum > (endingGenNum-1)) or (currentBestScore < minimalScore) or os.path.exists("./STOP")
    if stopCondition:
      globalVars.DONE = 1
      if globalVars.MOEA:
        makeNetlist_netlister(tempPool.fronts[0][0])
      else:
        makeNetlist_netlister(hotGen.pool[sortedPool_Indices[0]])
      #print hotGen.pool[sortedPool_Indices[0]].BigCircuitMatrix
      
  
  listener_thread.join(0)
  cOS.finalize()
  print("\n+++	 MATRIX EVOLUTIONS ENDED	+++")
  
