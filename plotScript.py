#Analyse - plot - MatrixRevolutions
"""
Plotting. Pretty dirty. 
"""

import matplotlib.pylab as plt
import matplotlib.gridspec as gridspec
import pickle
import numpy as np
from reproduction import fullRedundancyBigCircuitMatrix

with open("../_MAIN_data/data.pkl","r") as pkl_file:
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

def diversityPlot_pLP(generation, generationNum, bestScoresList, result, bestI):
  """Plots a diversity in generation in a shape of a matrix."""
  
  progress = bestScoresList
  
  #sum a generation into single matrix
  Msum = np.diag(np.zeros(BigMatrixSize,dtype=int),0)	
  for i in generation.pool:
	  Msum = Msum + i.BigCircuitMatrix
  Msum = Msum - np.diag(POP_SIZE*np.ones(BigMatrixSize,dtype=int),0)
  # Plot window
  fig = plt.figure(1, figsize=(14, 11), dpi=100, facecolor='w', edgecolor='k')
  
  # Create 4 subplots, 2x2
  IDs=			plt.subplot2grid((3,4), (0,1))
  mtrxDiver = 		plt.subplot2grid((3,4), (0,0))
  scoresPlt=		plt.subplot2grid((3,4), (0,2))
  mtrxDensPlt=		plt.subplot2grid((3,4), (0,3))
  
  progressPlt=		plt.subplot2grid((3,4), (1,0), colspan=4)
  
  filterPlt=		plt.subplot2grid((3,4), (2,0), colspan=3)
  filterMtrxPlot=	plt.subplot2grid((3,4), (2,3))
  
  plt.subplots_adjust(hspace=0.3)

  vectorValuesHASHES = []
  BigCircuitMatrixHASHES = []
  for i in range(0, POP_SIZE):
    vectorValuesHASHES.append(hash(generation.pool[i].ValueVector.tostring()))
    BigCircuitMatrixHASHES.append(hash(generation.pool[i].BigCircuitMatrix.tostring()))

  IDs.hold(True)
  IDs.set_title('Topology/values\n uniquiness')
  #IDs.plot(range(0,POP_SIZE), sorted(generation.matrixQuaziIDs), '-', label='mtrx qID')
  IDs.plot(range(0,POP_SIZE), sorted(BigCircuitMatrixHASHES), 'k-', label='mtrx qID')
  IDs.plot(range(0,POP_SIZE), sorted(vectorValuesHASHES), 'g-', label='values qID')
  IDs.grid(True)
  IDs.hold(False)

  mtrxDiver.hold(True)
  mtrxDiver.set_title('Last generation\n diversity')
  mtrxDiver.imshow(Msum, interpolation='nearest', cmap=plt.cm.OrRd)#Blues)
  mtrxDiver.hold(False)
  
  scoresPlt.hold(True)
  scoresPlt.set_title('Last generation\ncost values')
  scoresPlt.semilogy(range(0,POP_SIZE), sorted(generation.scores), 'x', label='GenScores')
  scoresPlt.grid(True)
  scoresPlt.hold(False)

  mtrxDensPlt.hold(True)
  mtrxDensPlt.set_title('Last generation\n matrix fillfactors')
  #mtrxDensPlt.plot(range(0,POP_SIZE), sorted(generation.matrixDensities), '-', label='mtrx density')
  mtrxDensPlt.grid(True)
  mtrxDensPlt.hold(False)

  #plot results (FILTER)
  filterPlt.hold(True)
  try:
    #print result['freq_scale']['default']
    #print result['tf']['default']
    #filterPlt.semilogx(kozaLPF()[1], kozaLPF()[0], '-', label='Koza response')
    filterPlt.set_title('Evolving filter bode plot')
    filterPlt.set_ylabel('gain [dB]')
    filterPlt.set_xlabel('frequency [Hz]')
    filterPlt.semilogx(result[1]['x']['nominal'], result[1]['y']['nominal'], '-', label='Frequency response')
    #filterPlt.semilogx(result['freq_scale']['default'], result['tf']['default'], '-', label='Frequency response')
    #filterPlt.semilogx(result['freq_scale']['default'], result['idealFilter'], '-', label='Ideal')
  except:
    print "No results to plot."
  filterPlt.grid(True)
  filterPlt.hold(False)

  progressPlt.hold(True)  
  progressPlt.set_title('Evolution progress - cost function over generations')
  progressPlt.set_ylabel('Cost value')  
  
  a, = progressPlt.semilogy(range(0,len(progress)), progress, '-', label='best')
  b, = progressPlt.semilogy(range(0,len(progress)), np.array(averageScore)/1, '-', label='average\n(w/o randoms)')
  #first_legend = progressPlt.legend(handles=[a,b], loc=1) 
  progressPlt.grid(True)
  progressPlt.hold(False)	
  
  filterMtrxPlot.hold(True)
  filterMtrxPlot.set_title("Evolving filter \nconnection matrix")
  filterMtrxPlot.imshow(generation.pool[bestI].BigCircuitMatrix, interpolation='nearest', cmap=plt.cm.Blues)
  filterMtrxPlot.hold(False)

  
  name = datadirname + "/" + "diversityPlots/gen_" + str(generationNum) + ".png" #ploting every generation in separate file
  plt.savefig(name)
  name = datadirname + "/" + "generationPlot" + ".png"
  plt.savefig(name)

def diversityPlot(generation, generationNum, bestScoresList, result, bestI):
  """Plots a diversity in generation in a shape of a matrix."""
  
  progress = bestScoresList
  
  #sum a generation into single matrix
  Msum = np.diag(np.zeros(BigMatrixSize,dtype=int),0)	
  for i in generation.pool:
	  Msum = Msum + i.BigCircuitMatrix
  Msum = Msum - np.diag(POP_SIZE*np.ones(BigMatrixSize,dtype=int),0)
  # Plot window
  fig = plt.figure(1, figsize=(14, 11), dpi=100, facecolor='w', edgecolor='k')
  
  # Create 4 subplots, 2x2
  IDs=			plt.subplot2grid((3,4), (0,1))
  mtrxDiver = 		plt.subplot2grid((3,4), (0,0))
  scoresPlt=		plt.subplot2grid((3,4), (0,2))
  mtrxDensPlt=		plt.subplot2grid((3,4), (0,3))
  
  progressPlt=		plt.subplot2grid((3,4), (1,0), colspan=4)
  
  vdd_sweep =		plt.subplot2grid((3,4), (2,0))
  rload_sweep =		plt.subplot2grid((3,4), (2,1))
  temp_sweep =		plt.subplot2grid((3,4), (2,2))
  
  filterMtrxPlot=	plt.subplot2grid((3,4), (2,3))
  
  plt.subplots_adjust(hspace=0.3)

  vectorValuesHASHES = []
  BigCircuitMatrixHASHES = []
  for i in range(0, POP_SIZE):
    vectorValuesHASHES.append(hash(generation.pool[i].ValueVector.tostring()))
    BigCircuitMatrixHASHES.append(hash(generation.pool[i].BigCircuitMatrix.tostring()))

  IDs.hold(True)
  IDs.set_title('Topology/values\n uniquiness')
  #IDs.plot(range(0,POP_SIZE), sorted(generation.matrixQuaziIDs), '-', label='mtrx qID')
  IDs.plot(range(0,POP_SIZE), sorted(BigCircuitMatrixHASHES), 'k-', label='mtrx qID')
  IDs.plot(range(0,POP_SIZE), sorted(vectorValuesHASHES), 'g-', label='values qID')
  IDs.grid(True)
  IDs.hold(False)

  mtrxDiver.hold(True)
  mtrxDiver.set_title('Last generation\n diversity')
  mtrxDiver.imshow(Msum, interpolation='nearest', cmap=plt.cm.OrRd)#Blues)
  mtrxDiver.hold(False)
  
  scoresPlt.hold(True)
  scoresPlt.set_title('Last generation\ncost values')
  scoresPlt.semilogy(range(0,POP_SIZE), sorted(generation.scores), 'x', label='GenScores')
  scoresPlt.grid(True)
  scoresPlt.hold(False)

  mtrxDensPlt.hold(True)
  mtrxDensPlt.set_title('Last generation\n matrix fillfactors')
  #mtrxDensPlt.plot(range(0,POP_SIZE), sorted(generation.matrixDensities), '-', label='mtrx density')
  mtrxDensPlt.grid(True)
  mtrxDensPlt.hold(False)

  #plot results (VOLTAGE REFERENCE)
  try:
    vdd_sweep.set_title('Vout(Vdd)')
    vdd_sweep.set_ylabel('[V]')
    vdd_sweep.set_xlabel('[V]')
    vdd_sweep.plot(result[1]['vout_vdd_scale']['nominal'], result[1]['vout_vdd']['nominal'], '-')
    
    rload_sweep.set_title('Vout(Rload)')
    rload_sweep.set_ylabel('[V]')
    rload_sweep.set_xlabel('[Ohm]')
    rload_sweep.semilogx(result[1]['vout_rload_scale']['nominal'], result[1]['vout_rload']['nominal'], '-')   
    
    temp_sweep.set_title('Vout(temp)')
    temp_sweep.set_ylabel('[V]')
    temp_sweep.set_xlabel('[deg C]')
    temp_sweep.plot(result[1]['vout_temp_scale']['nominal'], result[1]['vout_temp']['nominal'], '-')    
    
  except:
    print "No results to plot."
  vdd_sweep.grid(True)
  vdd_sweep.hold(False)
  
  rload_sweep.grid(True)
  rload_sweep.hold(False)
  
  temp_sweep.grid(True)
  temp_sweep.hold(False)

  progressPlt.hold(True)  
  progressPlt.set_title('Evolution progress - cost function over generations')
  progressPlt.set_ylabel('Cost value')  
  
  a, = progressPlt.semilogy(range(0,len(progress)), progress, '-', label='best')
  #b, = progressPlt.semilogy(range(0,len(progress)), np.array(averageScore)/1, '-', label='average\n(w/o randoms)')
  #first_legend = progressPlt.legend(handles=[a,b], loc=1) 
  progressPlt.grid(True)
  progressPlt.hold(False)	
  
  filterMtrxPlot.hold(True)
  filterMtrxPlot.set_title("Evolving circuit \nconnection matrix")
  filterMtrxPlot.imshow(generation.pool[bestI].BigCircuitMatrix, interpolation='nearest', cmap=plt.cm.Blues)
  filterMtrxPlot.hold(False)

  
  name = datadirname + "/" + "diversityPlots/gen_" + str(generationNum) + ".png" #ploting every generation in separate file
  plt.savefig(name)
  name = datadirname + "/" + "generationPlot" + ".png"
  plt.savefig(name)
  
def cmosVoltRef_evolutionPlot(generation, generationNum, bestScoresList, result, bestI):
  """Draw the plots of cmosVoltage reference evolution."""
  
  progress = bestScoresList
  
  #sum a generation into single matrix
  Msum = np.diag(np.zeros(BigMatrixSize,dtype=int),0)	
  for i in generation.pool:
	  Msum = Msum + i.BigCircuitMatrix
  Msum = Msum - np.diag(POP_SIZE*np.ones(BigMatrixSize,dtype=int),0)
  # Plot window
  fig = plt.figure(1, figsize=(14, 11), dpi=100, facecolor='w', edgecolor='k')
  
  # Create 4 subplots, 2x2
  IDs=			plt.subplot2grid((3,4), (0,1))
  mtrxDiver = 		plt.subplot2grid((3,4), (0,0))
  scoresPlt=		plt.subplot2grid((3,4), (0,2))
  mtrxDensPlt=		plt.subplot2grid((3,4), (0,3))
  
  progressPlt=		plt.subplot2grid((3,4), (1,0), colspan=3)
  filterMtrxPlot=	plt.subplot2grid((3,4), (1,3))
  
  
  vdd_sweep1 =		plt.subplot2grid((3,4), (2,0), colspan=2)
  vdd_sweep2 =		plt.subplot2grid((3,4), (2,2), colspan=2)  
  
  
  plt.subplots_adjust(hspace=0.3)

  vectorValuesHASHES = []
  BigCircuitMatrixHASHES = []
  for i in range(0, POP_SIZE):
    vectorValuesHASHES.append(hash(generation.pool[i].ValueVector.tostring()))
    BigCircuitMatrixHASHES.append(hash(generation.pool[i].BigCircuitMatrix.tostring()))

  IDs.hold(True)
  IDs.set_title('Topology/values\n uniquiness')
  #IDs.plot(range(0,POP_SIZE), sorted(generation.matrixQuaziIDs), '-', label='mtrx qID')
  IDs.plot(range(0,POP_SIZE), sorted(BigCircuitMatrixHASHES), 'k-', label='mtrx qID')
  IDs.plot(range(0,POP_SIZE), sorted(vectorValuesHASHES), 'g-', label='values qID')
  IDs.grid(True)
  IDs.hold(False)

  mtrxDiver.hold(True)
  mtrxDiver.set_title('Last generation\n diversity')
  mtrxDiver.imshow(Msum, interpolation='nearest', cmap=plt.cm.OrRd)#Blues)
  mtrxDiver.hold(False)
  
  scoresPlt.hold(True)
  scoresPlt.set_title('Last generation\ncost values')
  scoresPlt.semilogy(range(0,POP_SIZE), sorted(generation.scores), 'x', label='GenScores')
  scoresPlt.grid(True)
  scoresPlt.hold(False)

  mtrxDensPlt.hold(True)
  mtrxDensPlt.set_title('Last generation\n matrix fillfactors')
  #mtrxDensPlt.plot(range(0,POP_SIZE), sorted(generation.matrixDensities), '-', label='mtrx density')
  mtrxDensPlt.grid(True)
  mtrxDensPlt.hold(False)

  #plot results (VOLTAGE REFERENCE)
  try:
    vdd_sweep1.set_title('Vout(Vdd)(@Rloads)')
    vdd_sweep1.set_ylabel('[V]')
    vdd_sweep1.set_xlabel('[V]')
    vdd_sweep1.plot(result[1]['vout_vdd_res_scale']['nominal'], result[1]['vout_vdd_res1']['nominal'], '-')
    vdd_sweep1.plot(result[1]['vout_vdd_res_scale']['nominal'], result[1]['vout_vdd_res2']['nominal'], '-')
    vdd_sweep1.plot(result[1]['vout_vdd_res_scale']['nominal'], result[1]['vout_vdd_res3']['nominal'], '-')
                          
    vdd_sweep2.set_title('Vout(Vdd)(@temps)')
    vdd_sweep2.set_ylabel('[V]')
    vdd_sweep2.set_xlabel('[V]')
    vdd_sweep2.plot(result[1]['vout_vdd_scale']['nominal'], result[1]['vout_vdd_temp1']['nominal'], '-')
    vdd_sweep2.plot(result[1]['vout_vdd_scale']['nominal'], result[1]['vout_vdd_temp2']['nominal'], '-')
    vdd_sweep2.plot(result[1]['vout_vdd_scale']['nominal'], result[1]['vout_vdd_temp3']['nominal'], '-')
    
  except:
    print "No results to plot."
  vdd_sweep1.grid(True)
  vdd_sweep1.hold(False)
  
  vdd_sweep2.grid(True)
  vdd_sweep2.hold(False)
  
  progressPlt.hold(True)  
  progressPlt.set_title('Evolution progress - cost function over generations')
  progressPlt.set_ylabel('Cost value')  
  
  a, = progressPlt.semilogy(range(0,len(progress)), progress, '-', label='best')
  #b, = progressPlt.semilogy(range(0,len(progress)), np.array(averageScore)/1, '-', label='average\n(w/o randoms)')
  #first_legend = progressPlt.legend(handles=[a,b], loc=1) 
  progressPlt.grid(True)
  progressPlt.hold(False)	
  
  filterMtrxPlot.hold(True)
  filterMtrxPlot.set_title("Evolving circuit \nconnection matrix")
  filterMtrxPlot.imshow(generation.pool[bestI].BigCircuitMatrix, interpolation='nearest', cmap=plt.cm.Blues)
  filterMtrxPlot.hold(False)

  
  name = datadirname + "/" + "diversityPlots/gen_" + str(generationNum) + ".png" #ploting every generation in separate file
  plt.savefig(name)
  name = datadirname + "/" + "generationPlot" + ".png"
  plt.savefig(name)

def cmosVoltRef_evolutionPlot_MOEA(generation, generationNum, bestScoresList, result, bestI):
  """Draw the plots of cmosVoltage reference evolution."""
  
  progress = np.array(bestScoresList)
  allObjectiveScores = data[9]
  
  ##sum a generation into single matrix
  #Msum = np.diag(np.zeros(BigMatrixSize,dtype=int),0)	
  #for i in generation.pool:
#	  Msum = Msum + i.BigCircuitMatrix
  #Msum = Msum - np.diag(POP_SIZE*np.ones(BigMatrixSize,dtype=int),0)
  
  #sum a generation into single matrix
  Msum = np.diag(np.zeros(BigMatrixSize,dtype=int),0)	
  for i in generation.pool:
	  Msum = Msum + i.BigCircuitMatrix
  Msum = Msum - np.diag(len(generation.pool)*np.ones(BigMatrixSize,dtype=int),0)  
  
  # Plot window
  fig = plt.figure(1, figsize=(14, 11), dpi=100, facecolor='w', edgecolor='k')
  
  # Create 4 subplots, 2x2
  IDs=			plt.subplot2grid((3,4), (0,1))
  mtrxDiver = 		plt.subplot2grid((3,4), (0,0))
  scoresPlt=		plt.subplot2grid((3,4), (0,2))
  mtrxDensPlt=		plt.subplot2grid((3,4), (0,3))
  
  progressPlt=		plt.subplot2grid((3,4), (1,0), colspan=3)
  filterMtrxPlot=	plt.subplot2grid((3,4), (1,3))
  
  
  vdd_sweep1 =		plt.subplot2grid((3,4), (2,0), colspan=1)
  vdd_sweep2 =		plt.subplot2grid((3,4), (2,1), colspan=1)
  vdd_sweep3 =		plt.subplot2grid((3,4), (2,2), colspan=1)
  vdd_sweep4 =		plt.subplot2grid((3,4), (2,3), colspan=1)  
  
  plt.subplots_adjust(hspace=0.3)

  vectorValuesHASHES = []
  BigCircuitMatrixHASHES = []
  for i in range(0, POP_SIZE):
    vectorValuesHASHES.append(hash(generation.pool[i].ValueVector.tostring()))
    BigCircuitMatrixHASHES.append(hash(generation.pool[i].BigCircuitMatrix.tostring()))

  IDs.hold(True)
  IDs.set_title('Topology/values\n uniquiness')
  #IDs.plot(range(0,POP_SIZE), sorted(generation.matrixQuaziIDs), '-', label='mtrx qID')
  IDs.plot(range(0,POP_SIZE), sorted(BigCircuitMatrixHASHES), 'k-', label='mtrx qID')
  IDs.plot(range(0,POP_SIZE), sorted(vectorValuesHASHES), 'g-', label='values qID')
  IDs.grid(True)
  IDs.hold(False)

  mtrxDiver.hold(True)
  mtrxDiver.set_title('Last generation\n diversity')
  mtrxDiver.imshow(Msum, interpolation='nearest', cmap=plt.cm.OrRd)#Blues)
  mtrxDiver.hold(False)
  
  scoresPlt.hold(True)
  #scoresPlt.set_title('Last generation\ncost values')
  scoresPlt.set_title('Last generation\nVref/PSRR objectives')
  scoresPlt.plot(allObjectiveScores[:,0], allObjectiveScores[:,1], '*', label='Vref/PSRR')
  scoresPlt.grid(True)
  scoresPlt.hold(False)

  mtrxDensPlt.hold(True)
  #mtrxDensPlt.set_title('Last generation\n matrix fillfactors')
  mtrxDensPlt.set_title('Last generation\npower/PSRR objective')
  mtrxDensPlt.semilogy(allObjectiveScores[:,0], allObjectiveScores[:,2], '*', label='power/PSRR')
  mtrxDensPlt.grid(True)
  mtrxDensPlt.hold(False)

  #plot results (VOLTAGE REFERENCE)
  try:
    vdd_sweep1.set_title('Vout(Vdd)(@Temps)')
    vdd_sweep1.set_ylabel('[V]')
    vdd_sweep1.set_xlabel('[V]')
    vdd_sweep1.plot(result[0][1]['vout_vdd_scale']['nominal'], result[0][1]['vout_vdd_temp1']['nominal'], '-', linewidth=2.0)
    vdd_sweep1.plot(result[0][1]['vout_vdd_scale']['nominal'], result[0][1]['vout_vdd_temp2']['nominal'], '-', linewidth=2.0)
    vdd_sweep1.plot(result[0][1]['vout_vdd_scale']['nominal'], result[0][1]['vout_vdd_temp3']['nominal'], '-', linewidth=2.0)
                 
    vdd_sweep2.set_title('Vout(Vdd)(@Temps)')
    #vdd_sweep2.set_ylabel('[V]')
    vdd_sweep2.set_xlabel('[V]')
    vdd_sweep2.plot(result[1][1]['vout_vdd_scale']['nominal'], result[1][1]['vout_vdd_temp1']['nominal'], '-', linewidth=2.0)
    vdd_sweep2.plot(result[1][1]['vout_vdd_scale']['nominal'], result[1][1]['vout_vdd_temp2']['nominal'], '-', linewidth=2.0)
    vdd_sweep2.plot(result[1][1]['vout_vdd_scale']['nominal'], result[1][1]['vout_vdd_temp3']['nominal'], '-', linewidth=2.0)
    
    vdd_sweep3.set_title('Vout(Vdd)(@Temps)')
    #vdd_sweep3.set_ylabel('[V]')
    vdd_sweep3.set_xlabel('[V]')
    vdd_sweep3.plot(result[2][1]['vout_vdd_scale']['nominal'], result[2][1]['vout_vdd_temp1']['nominal'], '-', linewidth=2.0)
    vdd_sweep3.plot(result[2][1]['vout_vdd_scale']['nominal'], result[2][1]['vout_vdd_temp2']['nominal'], '-', linewidth=2.0)
    vdd_sweep3.plot(result[2][1]['vout_vdd_scale']['nominal'], result[2][1]['vout_vdd_temp3']['nominal'], '-', linewidth=2.0)
    
    vdd_sweep4.set_title('Vout(Vdd)(@Temps)')
    #vdd_sweep4.set_ylabel('[V]')
    vdd_sweep4.set_xlabel('[V]')
    vdd_sweep4.plot(result[3][1]['vout_vdd_scale']['nominal'], result[3][1]['vout_vdd_temp1']['nominal'], '-', linewidth=2.0)
    vdd_sweep4.plot(result[3][1]['vout_vdd_scale']['nominal'], result[3][1]['vout_vdd_temp2']['nominal'], '-', linewidth=2.0)
    vdd_sweep4.plot(result[3][1]['vout_vdd_scale']['nominal'], result[3][1]['vout_vdd_temp3']['nominal'], '-', linewidth=2.0)    
    
  except:
    print "No results to plot."
  vdd_sweep1.grid(True)
  vdd_sweep1.hold(False)
  
  vdd_sweep2.grid(True)
  vdd_sweep2.hold(False)
  
  vdd_sweep3.grid(True)
  vdd_sweep3.hold(False)
  
  vdd_sweep4.grid(True)
  vdd_sweep4.hold(False)  
  
  
  progressPlt.hold(True)  
  progressPlt.set_title('Evolution progress - objectives over generations')
  progressPlt.set_ylabel('Cost value')  

  progressPlt.semilogy(range(0,len(progress[2:,0])), progress[2:,0]/max(progress[2:,0]), '-', label='bestO1')
  progressPlt.semilogy(range(0,len(progress[2:,1])), progress[2:,1]/max(progress[2:,1]), '-', label='bestO2')
  progressPlt.semilogy(range(0,len(progress[2:,2])), progress[2:,2]/max(progress[2:,2]), '-', label='bestO3')
  #progressPlt.set_yscale('symlog')  
  
  #a, = progressPlt.semilogy(range(0,len(progress)), progress, '-', label='best')
  #b, = progressPlt.semilogy(range(0,len(progress)), np.array(averageScore)/1, '-', label='average\n(w/o randoms)')
  #first_legend = progressPlt.legend(handles=[a,b], loc=1) 
  progressPlt.grid(True)
  progressPlt.hold(False)	
  
  filterMtrxPlot.hold(True)
  filterMtrxPlot.set_title("Evolving circuit \nconnection matrix")
  filterMtrxPlot.imshow(generation.fronts[0][0].BigCircuitMatrix, interpolation='nearest', cmap=plt.cm.Blues)
  filterMtrxPlot.hold(False)

  
  name = datadirname + "/" + "diversityPlots/gen_" + str(generationNum) + ".png" #ploting every generation in separate file
  plt.savefig(name)
  name = datadirname + "/" + "generationPlot" + ".eps"
  plt.savefig(name)

def filterActiveLP_MOEA(generation, generationNum, bestScoresList, result, bestI):
  """Plots a diversity in generation and evolution progress of Multi-Objective Evolutionary Process."""
  
  progress = np.array(bestScoresList)
  allObjectiveScores = data[9]
  
  #sum a generation into single matrix
  Msum = np.diag(np.zeros(BigMatrixSize,dtype=int),0)	
  for i in generation.pool:
	  Msum = Msum + i.BigCircuitMatrix
  Msum = Msum - np.diag(len(generation.pool)*np.ones(BigMatrixSize,dtype=int),0)
  # Plot window
  fig = plt.figure(1, figsize=(14, 11), dpi=60, facecolor='w', edgecolor='k')
  
  # Create 4 subplots, 2x2
  IDs=			plt.subplot2grid((3,4), (0,1))
  mtrxDiver = 		plt.subplot2grid((3,4), (0,0))
  scoresPlt=		plt.subplot2grid((3,4), (0,2))
  mtrxDensPlt=		plt.subplot2grid((3,4), (0,3))
  
  progressPlt=		plt.subplot2grid((3,4), (1,0), colspan=4)
  
  filterPlt=		plt.subplot2grid((3,4), (2,0), colspan=3)
  filterMtrxPlot=	plt.subplot2grid((3,4), (2,3))
  
  plt.subplots_adjust(hspace=0.3)

  vectorValuesHASHES = []
  BigCircuitMatrixHASHES = []
  for i in range(0, POP_SIZE):
    vectorValuesHASHES.append(generation.pool[i].valuHash)
    BigCircuitMatrixHASHES.append(generation.pool[i].reduMtrxHash)

  IDs.hold(True)
  IDs.set_title('Topology/values\n uniquiness')
  #IDs.plot(range(0,POP_SIZE), sorted(generation.matrixQuaziIDs), '-', label='mtrx qID')
  IDs.plot(range(0,POP_SIZE), sorted(BigCircuitMatrixHASHES), 'k-', label='mtrx qID')
  IDs.plot(range(0,POP_SIZE), sorted(vectorValuesHASHES), 'g-', label='values qID')
  IDs.grid(True)
  IDs.hold(False)

  mtrxDiver.hold(True)
  mtrxDiver.set_title('Last generation\n diversity')
  mtrxDiver.imshow(Msum, interpolation='nearest', cmap=plt.cm.OrRd)#Blues)
  mtrxDiver.hold(False)
  
  scoresPlt.hold(True)
  #scoresPlt.set_title('Last generation\ncost values')
  scoresPlt.set_title('Last generation\nRipple/Slope objectives')
  scoresPlt.plot(allObjectiveScores[:,0], allObjectiveScores[:,1], '*', label='Ripple/Slope')
  scoresPlt.grid(True)
  scoresPlt.hold(False)

  mtrxDensPlt.hold(True)
  #mtrxDensPlt.set_title('Last generation\n matrix fillfactors')
  mtrxDensPlt.set_title('Last generation\nBw/Slope objective')
  mtrxDensPlt.plot(allObjectiveScores[:,0], allObjectiveScores[:,2], '*', label='Bw/Slope')
  mtrxDensPlt.grid(True)
  mtrxDensPlt.hold(False)

  #plot results (FILTER)
  filterPlt.hold(True)
  try:
    #print result['freq_scale']['default']
    #print result['tf']['default']
    #filterPlt.semilogx(kozaLPF()[1], kozaLPF()[0], '-', label='Koza response')
    filterPlt.set_title('Evolving filter Bode plot')
    filterPlt.set_ylabel('gain [dB]')
    filterPlt.set_xlabel('frequency [Hz]')
    
    for i in range(0, len(result)):
      filterPlt.semilogx(result[i][1]['x']['nominal'], result[i][1]['y']['nominal'], '-', label='Frequency response')
    #filterPlt.semilogx(result['freq_scale']['default'], result['tf']['default'], '-', label='Frequency response')
    #filterPlt.semilogx(result['freq_scale']['default'], result['idealFilter'], '-', label='Ideal')
  except:
    print "No results to plot."
  filterPlt.grid(True)
  filterPlt.hold(False)

  progressPlt.hold(True)  
  progressPlt.set_title('Evolution progress - cost function over generations')
  progressPlt.set_ylabel('Cost value')  
  
  progressPlt.plot(range(0,len(progress[:,0])), progress, '-', label='bestO1')
  progressPlt.plot(range(0,len(progress[:,1])), progress, '-', label='bestO2')
  progressPlt.plot(range(0,len(progress[:,2])), progress, '-', label='bestO3')
  #b, = progressPlt.semilogy(range(0,len(progress)), np.array(averageScore)/1, '-', label='average\n(w/o randoms)')
  #first_legend = progressPlt.legend(handles=[a,b, c], loc=1) 
  progressPlt.grid(True)
  progressPlt.hold(False)	
  
  filterMtrxPlot.hold(True)
  filterMtrxPlot.set_title("Evolving filter \nconnection matrix")
  filterMtrxPlot.imshow(generation.fronts[0][0].BigCircuitMatrix, interpolation='nearest', cmap=plt.cm.Blues)
  filterMtrxPlot.hold(False)

  
  
  name = datadirname + "/" + "diversityPlots/gen_" + str(generationNum) + ".png" #ploting every generation in separate file
  plt.savefig(name)
  name = datadirname + "/" + "generationPlot" + ".png"
  plt.savefig(name)

#diversityPlot(generation, generationNum, bestScoresList, result, bestI)
#diversityPlot_pLP(generation, generationNum, bestScoresList, result, bestI)
cmosVoltRef_evolutionPlot(generation, generationNum, bestScoresList, result, bestI)
#cmosVoltRef_evolutionPlot_MOEA(generation, generationNum, bestScoresList, result, bestI)

#filterActiveLP_MOEA(generation, generationNum, bestScoresList, result, bestI)
