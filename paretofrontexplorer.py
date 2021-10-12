#Pareto front explorer
"""
Use this script only for crawling ower obtained results. 
A .pkl file with all the data should exist in ../_MAIN_data/backdata.pkl. 


"""


import pickle
import numpy as np

from copy import copy
from utils import printer
import matplotlib.pylab as plt
import globalVars as GLOBAL 

from scorefunctions.arithmetic.scoreCirc_squareroot_resilenceMode import scoreCirc_squareroot_resilenceMode as PROBLEM

#open backdata.pkl
with open("../_MAIN_data/backdata.pkl","rb") as pkl_file:
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


def plot_sqrt(results):
  fig = plt.figure(1, figsize=(14, 11), dpi=100, facecolor='w', edgecolor='k')
  x = results[1][1]['scale']['nominal']
  for r in results[1]:
    x = r['scale']['nominal']
    y = r['vout']['nominal']
    plt.plot(x, y)


  plt.plot(x, np.sqrt(x), "*")
  plt.grid(True)

  #plt.show()
  name = datadirname + "/" + "diversityPlots/" + "test_pareto_search" + ".png" #ploting every generation in separate file
  plt.savefig(name)  

if __name__=='__main__':
  setOfChosen = []
  print(datadirname)

  eps1 = 2
  if gMOEA: 
    #plt.hold(True)
    for individual in generation.pool:
      if (
        #individual.objectivesScore[2] < 10 # score_array.std()
        # and
        #(individual.objectivesScore[1] < 250+eps1 and individual.objectivesScore[1] > 250-eps1) # score_array.sum()
        #and
        (individual.objectivesScore[0]) < 11 #+eps1 and individual.objectivesScore[0] > 10-eps1) # score_array[0]
        ):
          setOfChosen.append(individual)
          print(individual)
          print(individual.objectivesScore)
    
    print("This is ", len(setOfChosen), " individuals. Evaluating now.")
    for chosen in setOfChosen:
      print(chosen)
      print(chosen.objectivesScore)
      results = PROBLEM(chosen, 1234, 1234, 1)  
      print(results)
      input("Tisni!")

      #plot_sqrt(results)

      fig = plt.figure(1, figsize=(14, 11), dpi=100, facecolor='w', edgecolor='k')
      
      for r in results[1]:
        x = r['scale']['nominal']
        y = r['vout']['nominal']
        if x is None or y is None:
          print(x)
          print(y)
        plt.plot(x, y)



    
    x = results[1][1]['scale']['nominal']
    plt.plot(x, np.sqrt(x), "*")
    plt.grid(True)
    
    #plt.show()
    name = datadirname + "/" + "diversityPlots/" + "test_pareto_search" + ".png" #ploting every generation in separate file
    plt.savefig(name)  
    
  print("End.")
  
  
