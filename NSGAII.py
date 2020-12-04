"""Another try. MOEA non-dominating sorting."""
import numpy as np
import random
from globalVars import *
from time import time

def fast_nondominated_sort(popu):
  """
  Takes a population object and assigns non-dominant fronts to the population. 
  It works on arbitrarily long objectivesScore vector. It works with numpy arrays only!
  """
  
  stw0 = time()
  popu.fronts = []
  popu.fronts.append([])
  
  c1 = 0
  c2 = 0
  c3 = 0
  for individual in popu:
    individual.isDominatedByThatMany = 0
    #individual.dominatesThose = np.array([])#set()
    individual.dominatesThose = [] #TESTING
	
    for other_individual in popu:   
      if (min(individual.objectivesScore <= other_individual.objectivesScore)==1 and
	  sum(individual.objectivesScore < other_individual.objectivesScore)>=1):
	#individual.dominatesThose = np.append(individual.dominatesThose, other_individual) 
	individual.dominatesThose.append(other_individual) #TESTING
	#print "individual dominates"
	c1 += 1
      elif (min(other_individual.objectivesScore <= individual.objectivesScore)==1 and
	sum(other_individual.objectivesScore < individual.objectivesScore)>=1):
	individual.isDominatedByThatMany += 1
	c2 += 1
	#print "other dominates"
      c3 += 1
    if individual.isDominatedByThatMany == 0:
      #print popu.fronts
      popu.fronts[0].append(individual)
      individual.rank = 0
  stw1 = time() 
  print "Internal time calculation fast_nondominated_sort: %f s" %(stw1-stw0),
  print c1, c2, c3
  i = 0
  while len(popu.fronts[i]) > 0:
    temp = []
    for individual in popu.fronts[i]:
      for other_individual in individual.dominatesThose:
	other_individual.isDominatedByThatMany -= 1
	if other_individual.isDominatedByThatMany == 0:
	  other_individual.rank = i+1
	  temp.append(other_individual)
    i = i+1
    popu.fronts.append(temp)
   

def faster_nondominated_sort(popu):
  """
  Takes a population object and assigns non-dominant fronts to the population. 
  This version is derived from fast_nondominated_sort and is suited for 3 objectives and optimized for speed.
  Native python arrays are 3 times faster than numpy arrays. All comparisons are 3 times faster. 
  This whole procedure is then up to 20 times faster than fast_nondominated_sort.
  """
  
  stw0 = time()
  popu.fronts = []
  popu.fronts.append([])
  
  c1 = 0
  c2 = 0
  c3 = 0
  for individual in popu:
    individual.isDominatedByThatMany = 0
    individual.dominatesThose = [] #TESTING
    #TODO replace hardcoded comparison with implicit loop like:
    #any([i<=j for i,j in zip(a,b)])
    
    for other_individual in popu:
      
      if ((individual.objectivesScore[0] <= other_individual.objectivesScore[0] and 	#---
	   individual.objectivesScore[1] <= other_individual.objectivesScore[1] and 	#  |
	   individual.objectivesScore[2] <= other_individual.objectivesScore[2]) and	#---
	  (individual.objectivesScore[0] < other_individual.objectivesScore[0] or 	#---
	   individual.objectivesScore[1] < other_individual.objectivesScore[1] or 	#  |
	   individual.objectivesScore[2] < other_individual.objectivesScore[2])): 	#---
	individual.dominatesThose.append(other_individual) #TESTING
	#print "individual dominates"
	c1 += 1
      elif ((other_individual.objectivesScore[0] <= individual.objectivesScore[0] and   #---
	   other_individual.objectivesScore[1] <= individual.objectivesScore[1] and     #  |
	   other_individual.objectivesScore[2] <= individual.objectivesScore[2]) and    #---
	  (other_individual.objectivesScore[0] < individual.objectivesScore[0] or       #---
	   other_individual.objectivesScore[1] < individual.objectivesScore[1] or       #  |
	   other_individual.objectivesScore[2] < individual.objectivesScore[2])):       #---
	individual.isDominatedByThatMany += 1
	c2 += 1
	#print "other dominates"
      c3 += 1
    if individual.isDominatedByThatMany == 0:
      #print popu.fronts
      popu.fronts[0].append(individual)
      individual.rank = 0
  stw1 = time() 
  print "Internal time calculation faster_nondominated_sort: %f s" %(stw1-stw0),
  print c1, c2, c3
  i = 0
  while len(popu.fronts[i]) > 0:
    temp = []
    for individual in popu.fronts[i]:
      for other_individual in individual.dominatesThose:
	other_individual.isDominatedByThatMany -= 1
	if other_individual.isDominatedByThatMany == 0:
	  other_individual.rank = i+1
	  temp.append(other_individual)
    i = i+1
    popu.fronts.append(temp)
   



#-------

def faster_nondominated_sort2(popu):
  """
  Takes a population object and assigns non-dominant fronts to the population.
  This version is derived from fast_nondominated_sort and is suited for 3 objectives and optimized for $
  UPDATE: Suited for arbitrarily long abjectives list.
  Native python arrays are 3 times faster than numpy arrays. All comparisons are 3 times faster.
  This whole procedure is then up to 20 times faster than fast_nondominated_sort.

  UNDER CONSTRUCTION
  """

  stw0 = time()
  popu.fronts = []
  popu.fronts.append([])

  c1 = 0
  c2 = 0
  c3 = 0
  for individual in popu:
    individual.isDominatedByThatMany = 0
    individual.dominatesThose = [] #TESTING
    #TODO replace hardcoded comparison with implicit loop like:
    #any([i<=j for i,j in zip(a,b)])
    #it is 2 times slower than hardcoded one

    for other_individual in popu:
      if (all([i<=oi for i, oi in zip(individual.objectivesScore, other_individual.objectivesScore)])
	  and
	  any([i<oi for i, oi in zip(individual.objectivesScore, other_individual.objectivesScore)])):
        individual.dominatesThose.append(other_individual)
      elif ((other_individual.objectivesScore[0] <= individual.objectivesScore[0] and   #---
           other_individual.objectivesScore[1] <= individual.objectivesScore[1] and     #  |
           other_individual.objectivesScore[2] <= individual.objectivesScore[2]) and    #---
          (other_individual.objectivesScore[0] < individual.objectivesScore[0] or       #---
           other_individual.objectivesScore[1] < individual.objectivesScore[1] or       #  |
           other_individual.objectivesScore[2] < individual.objectivesScore[2])):       #---
        individual.isDominatedByThatMany += 1
        c2 += 1
        #print "other dominates"
      c3 += 1
    if individual.isDominatedByThatMany == 0:
      #print popu.fronts
      popu.fronts[0].append(individual)
      individual.rank = 0
  stw1 = time()
  print "Internal time calculation faster_nondominated_sort: %f s" %(stw1-stw0),
  print c1, c2, c3
  i = 0
  while len(popu.fronts[i]) > 0:
    temp = []
    for individual in popu.fronts[i]:
      for other_individual in individual.dominatesThose:
        other_individual.isDominatedByThatMany -= 1
        if other_individual.isDominatedByThatMany == 0:
          other_individual.rank = i+1
          temp.append(other_individual)
    i = i+1
    popu.fronts.append(temp)

    
def crowding_distance_assignment(front):
  """"""
  l = len(front)
  if l>0: #if any members in front
    for i in front:
      i.crow_dist = 0
    for m in range(0,len(front[0].objectivesScore)):
      #sort the objects in a front according to one of objectives (m)
      front.sort(key=lambda x: x.objectivesScore[m], reverse=False)
      front[0].crow_dist = np.Inf	#side objects of the front are assigned crowding distance being Inf
      front[l-1].crow_dist = np.Inf   
      for i in range(2, l-1):
	#print front[i+1].objectivesScore[m], front[i-1].objectivesScore[m], front[l-1].objectivesScore[m], front[0].objectivesScore[m]
	if front[l-1].objectivesScore[m]-front[0].objectivesScore[m] != 0:
	  front[i].crow_dist += ((front[i+1].objectivesScore[m]-front[i-1].objectivesScore[m])/
			      (front[l-1].objectivesScore[m]-front[0].objectivesScore[m]))
	  #Calculate crowding distance by summing crow_dist from each of objectives
	else:
	  if(debug>5):
	    print "\tcrowding_distance_assignment --> Division by zero avoided!"
	  front[i].crow_dist += np.Inf
	
    
def crowding_comparator(individual, other_individual):
  """"""
  if (individual.rank < other_individual.rank) or \
      ((individual.rank == other_individual.rank) and (individual.crow_dist > other_individual.crow_dist)):
    return 1
  else:
    return 0


def tournament_NSGA(generation, tournamentSize):
  """NonDominant CrowdingDistance awared tournament."""
  
  #print "\t-------Tournament"
  tour_Pool = []
  best = None
  
  #print len(generation.pool)
  #select several
  for i in range(0,tournamentSize):
    tour_Pool.append(generation.pool[random.randint(0,len(generation.pool)-1)])  
  #print "\t Individuals:", tour_Pool
  for individual in tour_Pool:
    if best is None or crowding_comparator(individual, best) == 1:
      best = individual
  #print "\t Tournament best:", best
  
  #print "\t-------Tournament end"
  return best

