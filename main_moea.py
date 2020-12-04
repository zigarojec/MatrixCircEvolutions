#Main script for multi-objective search - using Non-Dominated Sorting GA (NSGA-II)
#Evolutionary algorythm for automatic topology synthesis
#2015-2016 Ziga Rojec - EDA, FE, UL 

import numpy as np
import random
from time import time, strftime
from copy import copy, deepcopy
import collections, pickle, shutil
import os, sys
from os import mkdir

#Paralellism modules
from pyopus.parallel import base
from pyopus.parallel.cooperative import cOS
from pyopus.parallel.mpi import MPI

#Optimisation modules
from pyopus.optimizer.psade import ParallelSADE
from pyopus.optimizer.base import Reporter, CostCollector, RandomDelay

#My modules
import AdamAndEve as AE
#from globalVars import *
import globalVars
from reproduction import *
from NSGAII import faster_nondominated_sort,crowding_distance_assignment,tournament_NSGA
from scoreFunctions import *
from paramOptimizer import *

#Other settings
np.set_printoptions(threshold='nan', linewidth=1000)	#print whole matrix

random.seed(globalVars.seedN)	#Fixed seed
np.random.seed(globalVars.seedN)
#-define the cost function--------------------------
#PROBLEM = scoreCirc_ActiveFilter_MOEA #IMPORTANT! Set this also in circuitUnderOptimiser function!

problems = {'scoreCirc_ActiveFilter_2':scoreCirc_ActiveFilter_2,
	    'scoreCirc_CmosVoltageReference_2':scoreCirc_CmosVoltageReference_2,
	    'scoreCirc_PassiveBandPass':scoreCirc_PassiveBandPass,
	    'scoreCirc_HighPass':scoreCirc_HighPass,
	    } #set this also in paramOptimizer

PROBLEM = problems[globalVars.PROBLEMname]
MOEAMODE = 1 # DO NOT CHANGE

sys.setrecursionlimit(5000)


if __name__=='__main__':
  t0 = time()

  print
  print "++++++++++++++++++++++++++++++++++++++++++++"
  print "+++	 MATRIX (R)EVOLUTIONS STARTED	+++\n"
  
  startdate = strftime("%Y_%m_%d")
  starttime = strftime("%H-%M")
  print "Starting date: ", startdate, ", starting time: ", starttime

  os.chdir("../_MAIN_data")
  
  datadirname = "data_" + startdate + "_" + starttime
  os.mkdir(datadirname)
  os.mkdir(datadirname + "/" + "diversityPlots")
    
  shutil.copy2('../_MAIN_work/globalVars.py', datadirname + '/globalVars.py')	#copy current globalVars script for logging reasons
  shutil.copy2('../_MAIN_work/scoreFunctions.py', datadirname + '/scoreFunctions.py')
  output = open("data.pkl", "wb")
  os.chdir("../_MAIN_work")
  
  # Set up MPI for parallel computing

  cOS.setVM(MPI(mirrorMap={
    'models.inc':'.', 
    'topdc.cir':'.',
    'topdc_psrr.cir':'.',
    'mosmm.inc':'.',
    'cmos180n.lib':'.',
    #'testTopCirc_hspice_AF.cir':'.', 
    #'testTopCirc_hspice_AF_outimped.cir':'.'
  }))
  
  generationNum = 0
  optimizedIndivids = [] 	#PSADE parameter-optimized individuals

  hotGen = population() #hot generation (current)
  oldGen = population() #previous generation (old)
  
  currentBestScore = np.Inf
  bestScoresList = []
  averageScoresList = []
  
  if globalVars.continuee:
    with open("backdata.pkl","r") as pkl_file:
      data = pickle.load(pkl_file)
    pkl_file.close()
    print "Ressurecting old population. Press any to proceed..."

    hotGen = data[0]
    generationNum = 0#data[1]
    bestScoresList = data[2]
    result = data[3]
    bestI = data[4]
    #datadirname = data[5]
    BigMatrixSize = data[6]
    POP_SIZE = data[7]
    averageScoresList = data[8]
    raw_input("...old population resurrected.")
  else:
    print "Creating initial population. Press any to proceed..."
    for i in range(0,globalVars.POP_SIZE):
      hotGen.add_individual(createRandomBigCircuitMatrix(copy(createRandomValueVector())))
    
    if globalVars.insertAdam:
      hotGen.pool[0] = AE.adam 	#Insert an already designed circuit to optimise.
      hotGen.pool[0].fullRedundancyMatrix = fullRedundancyBigCircuitMatrix(AE.adam.BigCircuitMatrix)
      print "Adam-circuit inserted into population."
    raw_input("...initial population created.")
  
  #---EVALUATE & SORT INITIAL POPULATION---# 
  print "EVALUATING GENERATION %d" %generationNum
  stw0 = time()
  results=cOS.dispatch(jobList=((PROBLEM, [hotGen.pool[i], generationNum, i, True]) for i in range(0,len(hotGen.pool))), remote=True)
  results = np.array(results)
  
  #Put together objects and objectivesScores
  for i in range(0,len(hotGen.pool)):
    hotGen.pool[i].objectivesScore = np.transpose(results[:,0])[i]
  
  stw1 = time()
  print "Evaluation of initial population lasted for %f s" %(stw1-stw0)
  
  #Sort population
  staSort = time()
  faster_nondominated_sort(hotGen)
  endSort = time()
  print "Fast_nondominated_sort finished in %.2f s." %(endSort-staSort)
  for i in range(0, len(hotGen.fronts)):
    crowding_distance_assignment(hotGen.fronts[i])
  print("Crowding_distance_assignment finished.")
  
  #najboljsi v generaciji je...
  print ":::GENERATION %04.d - BEST ONE::: %s ::YEAH!::" %(generationNum,hotGen.pool[0].objectivesScore)
  #bestScoresList.append(hotGen.scores[sortedPool_Indices[0]])
  averageScoresList.append(np.average(hotGen.scores))
  
  #DONE = 0	#replaced with DONE from globalVars.py
  while globalVars.DONE == 0:
    generationNum += 1	#rise the generations counter
    
    oldGen = copy(hotGen)
    hotGen = population() #hot (current) generation initialization
    
    aa = []
    bb = []
    tempPool = population()
    offspring = []
    for i in range(0, globalVars.POP_SIZE/2):
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
    print "Offspring len before cleaning:", len(offspring)
    offspring = removeDuplicatesFromArrayByAttribute(offspring, "fullHash")
    print "Offspring len after cleaning:", len(offspring)
    
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
      #print offspring[i].objectivesScore
      #print offspring[i].reduMtrxHash
    if generationNum > 2:# and debug > 1:
      print "tempPool len before cleaning:", len(tempPool.pool)    
      tempPool.pool = removeDuplicatesFromArrayByAttribute(tempPool.pool, "scoreHash")     
      print "tempPool len after cleaning:", len(tempPool.pool)     
    
    
    #Sort population - fast-non-dominant sorting
    staSort = time()
    faster_nondominated_sort(tempPool)
    endSort = time()
    print "Fast_nondominated_sort finished in %.2fs." %(endSort-staSort)    
    
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
    
    #---Saving results, picle, ...  and so on----#
    #Calculate time
    stw1 = time()		
    dtime = stw1-stw0
    m, s = divmod(dtime, 60)
    h, m = divmod(m, 60)
    
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
	  print "PSADE on circuit ", i, "from gen", generationNum, "..."
	  maxiter = 3000 if generationNum < 100 else 8000
	  x, f = optimiseCircuit(topology, values, maxiter)
	  optimizedIndivids.append(copy(circuit(topology, x)))
	  
	  #append those individuals to the NEXT generation!!
    else:
      optimizedIndivids = []
    #------------------------------------------------#          
    
    
    plotThisMany = 4
    #if len(tempPool.fronts[0]) < plotThisMany:
      #plotThisMany = len(tempPool.fronts[0])
      
    #Evaluate best one/few together with results for plotting
    WinnerResults = cOS.dispatch(jobList=((PROBLEM, [tempPool.pool[distxyzI[i]], generationNum, 0, True]) for i in range(0,plotThisMany)), remote=True)
    WinnerResults = np.array(WinnerResults)
    
    printer(WinnerResults[0], stw0, generationNum, problem = globalVars.PROBLEMname) #print results to screen
    
    print "\t - Unique of scores:",len(set(a)),"(/",firstFlen, "), of topologies:",len(set(b)),"(/",firstFlen, "), \n\t          of values:", len(set(c)),"(/",firstFlen, "), of individuals", len(set(d)),"(/",firstFlen, ") in 1st front."
    print "\t- - - - - - - - - - - - - - - - - - - - - - - - - - "
    
    #Write winner netlist to a directiory of current run for inspection and manual simulation 
    os.chdir("../_MAIN_data")
    os.chdir("data_" + startdate + "_" + starttime)

    writeThisMany = 4
    for i in range(0, writeThisMany):
      makeNetlist(tempPool.pool[distxyzI[i]], 0, i, tempPool.pool[distxyzI[i]].fullRedundancyMatrix)
    os.chdir("../") 

    #DUMP results for plotting
    allObjectiveValues = hotGen.pool[0].objectivesScore		#getting all objective values of hotGen for plotting
    for i in range(1, len(hotGen.pool)):
      allObjectiveValues = np.vstack((allObjectiveValues, hotGen.pool[i].objectivesScore))
    
    data = [tempPool, generationNum, bestScoresList, WinnerResults, None, datadirname, BigMatrixSize, globalVars.POP_SIZE, averageScoresList, allObjectiveValues]
    with open(datadirname + "/backdata.pkl","wb") as output:
      pickle.dump(data, output)
    output.close()
    #shutil.copy2("../_MAIN_work/" + 'data.pkl', datadirname + '/backdata.pkl')	#copy current .pkl data to current run folder WARNING: not healthy when running multithread #2018_04_03
    shutil.copy2(datadirname + "/backdata.pkl", "../_MAIN_work/" + "data.pkl")
    
    os.chdir("../_MAIN_work")
    #End of evolution? :
    if generationNum > (globalVars.endingGenNum-1):
      globalVars.DONE = 1
      #makeNetlist(tempPool.fronts[0][0], generationNum, 0, fullRedundancyBigCircuitMatrix(tempPool.fronts[0][0].BigCircuitMatrix))
      
  
  cOS.finalize()
  print "\n+++	 MATRIX EVOLUTIONS ENDED	 +++"
  print "+++++++++++++++++++++++++++++++++++++++++++"
  