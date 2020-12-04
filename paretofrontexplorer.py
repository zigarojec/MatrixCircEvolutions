#Pareto front explorer
"""
Use this script only for crawling ower obtained results. 
A .pkl file with all the data should exist in ../_MAIN_data/backdata.pkl. 


"""


import pickle
import numpy as np
from scoreFunctions import scoreCirc_CmosVoltageReference_2
from copy import copy
from utils import printer
import matplotlib.pylab as plt

#open backdata.pkl
with open("../_MAIN_data/backdata.pkl","r") as pkl_file:
  data = pickle.load(pkl_file)
pkl_file.close()

generation = data[0]
generationNum = data[1]
bestScoresList = data[2]
result = data[3]
bestI = data[4]
datadirname = data[5]
BigMatrixSize = data[6]
POP_SIZE = data[7]
averageScore = data[8]

datadirname = '../_MAIN_data/' + datadirname

gMOEA = 1


def plotTemp(results):
  x = results[1]['vout_vdd_scale']['nominal']
  y = results[1]['vout_vdd_temp1']['nominal']
  plt.hold(True)
  plt.plot(x, y,linewidth=2)
  y = results[1]['vout_vdd_temp2']['nominal']
  plt.plot(x, y,linewidth=2)
  y = results[1]['vout_vdd_temp3']['nominal']
  plt.plot(x, y,linewidth=2)
  plt.grid(True)
  plt.hold(False)
  plt.show()
  
def plotTran(results):
  x = results[1]['t']['nominal']
  yin = results[1]['vin_psrr']['nominal']
  yout = results[1]['vout_psrr']['nominal']
  plt.hold(True)
  plt.plot(x, yin)
  plt.plot(x, yout)
  plt.grid(True)
  plt.hold(False)
  plt.show()

if __name__=='__main__':
  setOfChosen = []
  print datadirname
  if gMOEA: 
    #plt.hold(True)
    for individual in generation.pool:
      if (individual.objectivesScore[2] < 0.0019 #oP (power)
	  and
	  individual.objectivesScore[1] < 0.015 #oPsrr (P-P value)
	  and
	  individual.objectivesScore[0] < 0.19 #oMediana (abs(x - VREF) at 3 temps @11V)
	  ):
	setOfChosen.append(individual)
	print individual
	print individual.objectivesScore
	lastBestIndividual = copy(individual)

	if len(setOfChosen):
	  results = scoreCirc_CmosVoltageReference_2(lastBestIndividual, 1234, 1234, 1)
	
	  if 1:
	    printer(results, 0, generationNum, problem = 'scoreCirc_CmosVoltageReference_2')
	    plotTemp(results)
	    #plotTran(results)
    #plt.hold(False)
    #plt.show() 

  else:
    bestIndeces = np.argsort(generation.scores)
    
    #plt.hold(True)
    
    for i in bestIndeces:
      results = scoreCirc_CmosVoltageReference_2(generation.pool[i], 1234, 1234, 0)
      printer(results, 0, generationNum, problem = 'scoreCirc_CmosVoltageReference_2')
      plotTemp(results)
      #plotTran(results)

    #plt.hold(False)
    #plt.show()  
    
  print "End."
  print "This is ", len(setOfChosen), " individuals."
  
  
