import matplotlib.pylab as plt
import matplotlib.gridspec as gridspec
import matplotlib.ticker 
import pickle
import numpy as np
from reproduction import fullRedundancyBigCircuitMatrix
import globalVars as GLOBAL

#with open("../_MAIN_data/backdata.pkl","rb") as pkl_file:
with open("./backdata.pkl","rb") as pkl_file:
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

if GLOBAL.MOEA:
  allObjectiveScores = data[9]

datadirname = '../_MAIN_data/' + datadirname
#datadirname = "."


def logcirc_MOEA(generation, generationNum, bestScoresList, result, bestI):
  """Plots a diversity in generation and evolution progress of Multi-Objective Evolutionary Process."""


  progress = np.array(bestScoresList)
  allObjectiveScores = data[9]
  
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
  
  progressPlt=		plt.subplot2grid((3,4), (1,0), colspan=4)
  
  errPlt=		plt.subplot2grid((3,4), (2,0), colspan=2)
  voutPlt=	plt.subplot2grid((3,4), (2,2), colspan=2)
  
  plt.subplots_adjust(hspace=0.3)

  vectorValuesHASHES = []
  BigCircuitMatrixHASHES = []
  for i in range(0, POP_SIZE):
    vectorValuesHASHES.append(generation.pool[i].valuHash)
    BigCircuitMatrixHASHES.append(generation.pool[i].reduMtrxHash)


  IDs.set_title('Topology/values\n uniquiness')
  #IDs.plot(range(0,POP_SIZE), sorted(generation.matrixQuaziIDs), '-', label='mtrx qID')
  IDs.plot(range(0,POP_SIZE), sorted(BigCircuitMatrixHASHES), 'k-', label='mtrx qID')
  IDs.plot(range(0,POP_SIZE), sorted(vectorValuesHASHES), 'g-', label='values qID')
  IDs.grid(True)


  mtrxDiver.set_title('Last generation\n diversity')
  mtrxDiver.imshow(Msum, interpolation='nearest', cmap=plt.cm.OrRd)#Blues)
  
  #scoresPlt.set_title('Last generation\ncost values')
  scoresPlt.set_title('Last generation\nFullDesign/Nominal objectives')
  scoresPlt.semilogy(allObjectiveScores[:,0], allObjectiveScores[:,1], '*', label='FullDesign/Nominal')
  scoresPlt.grid(True)

  #mtrxDensPlt.set_title('Last generation\n matrix fillfactors')
  mtrxDensPlt.set_title('Last generation\nFailuresSTD/Nominal objective')
  mtrxDensPlt.plot(allObjectiveScores[:,0], allObjectiveScores[:,2], '*', label='FailuresSTD/Nominal')
  mtrxDensPlt.grid(True)

  
  errPlt.set_title('Err (Vin)')
  errPlt.set_ylabel('[V]')
  errPlt.set_xlabel('[V]')

  
  if not GLOBAL.robustMode:
    target = 2*np.log(result[1]['scale']['nominal'] + 1)
    errPlt.plot(result[1]['scale']['nominal'], target-result[1]['vout']['nominal'], '-')
  else:
    target = 2*np.log(result[0][1][0]['scale']['nominal'] + 1)
    first = True
    for r in result[0][1]:
      if first:
        errPlt.plot(r['scale']['nominal'], target-r['vout']['nominal'], '*')
        first = False
      else:
        errPlt.plot(r['scale']['nominal'], target-r['vout']['nominal'], '-')

  errPlt.grid(True)

  voutPlt.set_title('Vout(Vin)')
  voutPlt.set_ylabel('[V]')
  voutPlt.set_xlabel('[V]')
  
  if not GLOBAL.robustMode:
    voutPlt.plot(result[1]['scale']['nominal'], target, '.')
    voutPlt.plot(result[1]['scale']['nominal'], result[1]['vout']['nominal'], '-')
  else:
    voutPlt.plot(result[0][1][0]['scale']['nominal'], 2*np.log(result[0][1][0]['scale']['nominal'] + 1), '.')
    for r in result[0][1]:
      voutPlt.plot(r['scale']['nominal'], r['vout']['nominal'], '-')

  voutPlt.grid(True)

  progressPlt.set_title('Evolution progress - cost function over generations')
  progressPlt.set_ylabel('Cost value')  
  
  progressPlt.semilogy(range(0,len(progress[:,0])), progress, '-', label='bestO1')
  progressPlt.semilogy(range(0,len(progress[:,1])), progress, '-', label='bestO2')
  progressPlt.semilogy(range(0,len(progress[:,2])), progress, '-', label='bestO3')
  #b, = progressPlt.semilogy(range(0,len(progress)), np.array(averageScore)/1, '-', label='average\n(w/o randoms)')
  #first_legend = progressPlt.legend(handles=[a,b, c], loc=1) 
  progressPlt.grid(True)

  

  
  
  name = datadirname + "/" + "diversityPlots/gen_" + str(generationNum) + ".png" #ploting every generation in separate file
  plt.savefig(name)
  #name = datadirname + "/" + "generationPlot" + ".png"
  #plt.savefig(name)


def abscirc_MOEA(generation, generationNum, bestScoresList, result, bestI):
  """Plots a diversity in generation and evolution progress of Multi-Objective Evolutionary Process."""


  progress = np.array(bestScoresList)
  allObjectiveScores = data[9]
  
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
  
  progressPlt=		plt.subplot2grid((3,4), (1,0), colspan=4)
  
  errPlt=		plt.subplot2grid((3,4), (2,0), colspan=2)
  voutPlt=	plt.subplot2grid((3,4), (2,2), colspan=2)
  
  plt.subplots_adjust(hspace=0.3)

  vectorValuesHASHES = []
  BigCircuitMatrixHASHES = []
  for i in range(0, POP_SIZE):
    vectorValuesHASHES.append(generation.pool[i].valuHash)
    BigCircuitMatrixHASHES.append(generation.pool[i].reduMtrxHash)


  IDs.set_title('Topology/values\n uniquiness')
  #IDs.plot(range(0,POP_SIZE), sorted(generation.matrixQuaziIDs), '-', label='mtrx qID')
  IDs.plot(range(0,POP_SIZE), sorted(BigCircuitMatrixHASHES), 'k-', label='mtrx qID')
  IDs.plot(range(0,POP_SIZE), sorted(vectorValuesHASHES), 'g-', label='values qID')
  IDs.grid(True)


  mtrxDiver.set_title('Last generation\n diversity')
  mtrxDiver.imshow(Msum, interpolation='nearest', cmap=plt.cm.OrRd)#Blues)
  
  #scoresPlt.set_title('Last generation\ncost values')
  scoresPlt.set_title('Last generation\nFullDesign/Nominal objectives')
  scoresPlt.semilogy(allObjectiveScores[:,0], allObjectiveScores[:,1], '*', label='FullDesign/Nominal')
  scoresPlt.grid(True)

  #mtrxDensPlt.set_title('Last generation\n matrix fillfactors')
  mtrxDensPlt.set_title('Last generation\nFailuresSTD/Nominal objective')
  mtrxDensPlt.plot(allObjectiveScores[:,0], allObjectiveScores[:,2], '*', label='FailuresSTD/Nominal')
  mtrxDensPlt.grid(True)

  
  errPlt.set_title('Err (Vin)')
  errPlt.set_ylabel('[V]')
  errPlt.set_xlabel('[V]')

  
  if not GLOBAL.robustMode:
    target = np.abs(result[1]['scale']['nominal'])
    errPlt.plot(result[1]['scale']['nominal'], target-result[1]['vout']['nominal'], '-')
  else:
    target = np.abs(result[0][1][0]['scale']['nominal'])
    first = True
    for r in result[0][1]:
      if first:
        errPlt.plot(r['scale']['nominal'], target-r['vout']['nominal'], '*')
        first = False
      else:
        errPlt.plot(r['scale']['nominal'], target-r['vout']['nominal'], '-')

  errPlt.grid(True)

  voutPlt.set_title('Vout(Vin)')
  voutPlt.set_ylabel('[V]')
  voutPlt.set_xlabel('[V]')
  
  if not GLOBAL.robustMode:
    voutPlt.plot(result[1]['scale']['nominal'], target, '.')
    voutPlt.plot(result[1]['scale']['nominal'], result[1]['vout']['nominal'], '-')
  else:
    voutPlt.plot(result[0][1][0]['scale']['nominal'], np.abs(result[0][1][0]['scale']['nominal']), '.')
    for r in result[0][1]:
      voutPlt.plot(r['scale']['nominal'], r['vout']['nominal'], '-')

  voutPlt.grid(True)

  progressPlt.set_title('Evolution progress - cost function over generations')
  progressPlt.set_ylabel('Cost value')  
  
  progressPlt.semilogy(range(0,len(progress[:,0])), progress, '-', label='bestO1')
  progressPlt.semilogy(range(0,len(progress[:,1])), progress, '-', label='bestO2')
  progressPlt.semilogy(range(0,len(progress[:,2])), progress, '-', label='bestO3')
  #b, = progressPlt.semilogy(range(0,len(progress)), np.array(averageScore)/1, '-', label='average\n(w/o randoms)')
  #first_legend = progressPlt.legend(handles=[a,b, c], loc=1) 
  progressPlt.grid(True)

  

  
  
  name = datadirname + "/" + "diversityPlots/gen_" + str(generationNum) + ".png" #ploting every generation in separate file
  plt.savefig(name)
  #name = datadirname + "/" + "generationPlot" + ".png"
  #plt.savefig(name)





def duoplot_log(generation, generationNum, bestScoresList, result, bestI):
    """
    """

    fig, (voutPlt, errPlt) = plt.subplots(nrows=1, ncols=2, sharex=True,
                                        figsize=(12, 6))

    errPlt.set_title('Err (Vin)')
    errPlt.set_ylabel('[V]')
    errPlt.set_xlabel('[V]')

    
    if not GLOBAL.robustMode:
        target = 2*np.log(result[1]['scale']['nominal'] + 1)
        errPlt.plot(result[1]['scale']['nominal'], target-result[1]['vout']['nominal'], '-')
    else:
        target = 2*np.log(result[0][1][0]['scale']['nominal'] + 1)
        first = True
        for r in result[0][1]:
          if first:
              errPlt.plot(r['scale']['nominal'], target-r['vout']['nominal'], '*')
              first = False
          else:
              errPlt.plot(r['scale']['nominal'], target-r['vout']['nominal'], '-')

    errPlt.grid(True)

    voutPlt.set_title('Vout(Vin)')
    voutPlt.set_ylabel('[V]')
    voutPlt.set_xlabel('[V]')
    
    if not GLOBAL.robustMode:
        voutPlt.plot(result[1]['scale']['nominal'], target, '.')
        voutPlt.plot(result[1]['scale']['nominal'], result[1]['vout']['nominal'], '-')
    else:
        voutPlt.plot(result[0][1][0]['scale']['nominal'], 2*np.log(result[0][1][0]['scale']['nominal'] + 1), '.')
        for r in result[0][1]:
          voutPlt.plot(r['scale']['nominal'], r['vout']['nominal'], '-')

    voutPlt.grid(True)

    name = datadirname + "/" + "diversityPlots/gen_" + str(generationNum) + "_duoplot_vout.png" #ploting every generation in separate file
    plt.savefig(name)
    name = datadirname + "/" + "generationPlot" + ".png"
    plt.savefig(name)    

def duoplot_sqrt(generation, generationNum, bestScoresList, result, bestI):
    """
    """
    def errorfun(reality, target):
      #reality = reality+(target[0]-reality[0])
      reality = reality + OFFSET
      #target = target-target[0]
      return ((reality-target)/np.maximum.reduce([abs(reality), abs(target)]))*100

    fig, (voutPlt, errPlt) = plt.subplots(nrows=1, ncols=2, sharex=True,
                                        figsize=(12, 4))

    #errPlt.set_title('Err(Vin)=sqrt(Vin)-Vout(Vin)')
    errPlt.set_title('RelErr($V_{in}$)')
    #errPlt.set_ylabel('Err [V]')
    errPlt.set_ylabel('RelErr [%]')
    errPlt.set_xlabel('$V_{in} [V]$')
    errPlt.set_ylim([-30, 80])

    
    if not GLOBAL.robustMode:
        target = np.sqrt(result[1]['scale']['nominal'])
        errPlt.plot(result[1]['scale']['nominal'], target-result[1]['vout']['nominal'], '-')
    else:
        target = np.sqrt(result[0][1][0]['scale']['nominal'])
        OFFSET = target[0] - result[0][1][0]['vout']['nominal'][0]
        print(OFFSET)
        first = True
        nexxt = None # 2 or 3 , 2st, 3nd
        count = 0
        for r in result[0][1]:
          if first:
              errPlt.plot(r['scale']['nominal'], errorfun(r['vout']['nominal'], target), 
              color = 'black', linestyle = '-', marker='', label='nominal', linewidth=3)
              first = False
              nexxt = 2
          elif nexxt == 2:
            count+=1
            if count <= 2:
              errPlt.plot(r['scale']['nominal'], errorfun(r['vout']['nominal'], target), 
              color = 'gray', linestyle = '--', marker='', label="failure_himp")
            else:
              errPlt.plot(r['scale']['nominal'], errorfun(r['vout']['nominal'], target), 
              color = 'gray', linestyle = '--', marker='')
            nexxt = 3
          elif nexxt == 3:
            errPlt.plot(r['scale']['nominal'], errorfun(r['vout']['nominal'], target), 
            color = 'gray', linestyle = ':', marker='')
            count+=1
            if count <= 2:
              errPlt.plot(r['scale']['nominal'], errorfun(r['vout']['nominal'], target), 
              color = 'gray', linestyle = ':', marker='', label= 'failure_sck')
            else:
              errPlt.plot(r['scale']['nominal'], errorfun(r['vout']['nominal'], target), 
              color = 'gray', linestyle = ':', marker='')
            nexxt = 2
          else:
              errPlt.plot(r['scale']['nominal'], errorfun(r['vout']['nominal'], target), 
              color = 'gray', linestyle = '--', marker='', label='failure')

    errPlt.grid(True)
    errPlt.legend()

    voutPlt.set_title('$V_{out}(V_{in})$')
    voutPlt.set_ylabel('$V_{out} [V]$')
    voutPlt.set_xlabel('$V_{in} [V]$')
    
    if not GLOBAL.robustMode:
        voutPlt.plot(result[1]['scale']['nominal'], target, '.')
        voutPlt.plot(result[1]['scale']['nominal'], result[1]['vout']['nominal'], '-')
    else:
        voutPlt.plot(result[0][1][0]['scale']['nominal'], 
                    np.sqrt(result[0][1][0]['scale']['nominal'] ), 
                    color = 'red', linestyle = '-.', marker='', label='sqrt($V_{in}$)')
        first = True
        nexxt = None # 2 or 3 , 2st, 3nd
        count = 0
        multidim = np.empty((0,len(r['scale']['nominal'])), float)
        for r in result[0][1]:
          #print(r['vout']['nominal'], len(r['vout']['nominal']))
          multidim = np.row_stack((multidim, r['vout']['nominal']))

          if first:
            voutPlt.plot(r['scale']['nominal'], r['vout']['nominal'], 
            color = 'black', linestyle = '-', marker='', label='nominal', linewidth=3)
            first = False
            nexxt = 2
          elif nexxt == 2:
            count+=1
            if count <= 2:
              None #voutPlt.plot(r['scale']['nominal'], r['vout']['nominal'], 
              #color = 'gray', linestyle = '--', marker='', label="failure_himp")
            else:
              None #voutPlt.plot(r['scale']['nominal'], r['vout']['nominal'], 
              #color = 'gray', linestyle = '--', marker='')
            nexxt = 3
          elif nexxt == 3:
            None #voutPlt.plot(r['scale']['nominal'], r['vout']['nominal'], 
            #color = 'gray', linestyle = ':', marker='')
            count+=1
            if count <= 2:
              None #voutPlt.plot(r['scale']['nominal'], r['vout']['nominal'], 
              #color = 'gray', linestyle = ':', marker='', label= 'failure_sck')
            else:
              None #voutPlt.plot(r['scale']['nominal'], r['vout']['nominal'], 
              #color = 'gray', linestyle = ':', marker='')
            nexxt = 2
          else:
            None #voutPlt.plot(r['scale']['nominal'], r['vout']['nominal'], 
            #color = 'gray', linestyle = '', marker='o', label='failure')

        errmax = np.amax(multidim, axis=0)
        errmin = np.amin(multidim, axis=0)
        voutPlt.fill_between(r['scale']['nominal'], errmin, errmax, hatch="||||", alpha=0.4, label="Failure range")

    voutPlt.grid(True)

    voutPlt.legend()

    name = datadirname + "/" + "diversityPlots/gen_" + str(generationNum) + "_duoplot_vout.eps" #ploting every generation in separate file
    plt.savefig(name)
    name = datadirname + "/" + "diversityPlots/gen_" + str(generationNum) + "_duoplot_vout.pdf" #ploting every generation in separate file
    plt.savefig(name)    
    #name = datadirname + "/" + "generationPlot" + ".png"
    #plt.savefig(name) 
    
    plt.close()


def duoplot_abs(generation, generationNum, bestScoresList, result, bestI):
    """
    """
    def errorfun(reality, target):
      #reality = reality+(target[0]-reality[0])
      reality = reality + OFFSET
      #target = target-target[0]
      return ((reality-target)/np.maximum.reduce([abs(reality), abs(target)]))*100

    fig, (voutPlt, errPlt) = plt.subplots(nrows=1, ncols=2, sharex=True,
                                        figsize=(12, 4))

    #errPlt.set_title('Err(Vin)=sqrt(Vin)-Vout(Vin)')
    errPlt.set_title('RelErr($V_{in}$)')
    #errPlt.set_ylabel('Err [V]')
    errPlt.set_ylabel('RelErr [%]')
    errPlt.set_xlabel('$V_{in} [V]$')
    errPlt.set_ylim([-30, 80])

    
    if not GLOBAL.robustMode:
        target = np.abs(result[1]['scale']['nominal'])
        errPlt.plot(result[1]['scale']['nominal'], target-result[1]['vout']['nominal'], '-')
    else:
        target = np.abs(result[0][1][0]['scale']['nominal'])
        OFFSET = target[0] - result[0][1][0]['vout']['nominal'][0]
        print(OFFSET)
        first = True
        nexxt = None # 2 or 3 , 2st, 3nd
        count = 0
        for r in result[0][1]:
          if first:
              errPlt.plot(r['scale']['nominal'], errorfun(r['vout']['nominal'], target), 
              color = 'black', linestyle = '-', marker='', label='nominal', linewidth=3)
              first = False
              nexxt = 2
          elif nexxt == 2:
            count+=1
            if count <= 2:
              errPlt.plot(r['scale']['nominal'], errorfun(r['vout']['nominal'], target), 
              color = 'gray', linestyle = '--', marker='', label="failure_himp")
            else:
              errPlt.plot(r['scale']['nominal'], errorfun(r['vout']['nominal'], target), 
              color = 'gray', linestyle = '--', marker='')
            nexxt = 3
          elif nexxt == 3:
            errPlt.plot(r['scale']['nominal'], errorfun(r['vout']['nominal'], target), 
            color = 'gray', linestyle = ':', marker='')
            count+=1
            if count <= 2:
              errPlt.plot(r['scale']['nominal'], errorfun(r['vout']['nominal'], target), 
              color = 'gray', linestyle = ':', marker='', label= 'failure_sck')
            else:
              errPlt.plot(r['scale']['nominal'], errorfun(r['vout']['nominal'], target), 
              color = 'gray', linestyle = ':', marker='')
            nexxt = 2
          else:
              errPlt.plot(r['scale']['nominal'], errorfun(r['vout']['nominal'], target), 
              color = 'gray', linestyle = '--', marker='', label='failure')

    errPlt.grid(True)
    errPlt.legend()

    voutPlt.set_title('$V_{out}(V_{in})$')
    voutPlt.set_ylabel('$V_{out} [V]$')
    voutPlt.set_xlabel('$V_{in} [V]$')
    
    if not GLOBAL.robustMode:
        voutPlt.plot(result[1]['scale']['nominal'], target, '.')
        voutPlt.plot(result[1]['scale']['nominal'], result[1]['vout']['nominal'], '-')
    else:
        voutPlt.plot(result[0][1][0]['scale']['nominal'], 
                    np.abs(result[0][1][0]['scale']['nominal'] ), 
                    color = 'red', linestyle = '-.', marker='', label='sqrt($V_{in}$)')
        first = True
        nexxt = None # 2 or 3 , 2st, 3nd
        count = 0
        multidim = np.empty((0,len(r['scale']['nominal'])), float)
        for r in result[0][1]:
          #print(r['vout']['nominal'], len(r['vout']['nominal']))
          multidim = np.row_stack((multidim, r['vout']['nominal']))

          if first:
            voutPlt.plot(r['scale']['nominal'], r['vout']['nominal'], 
            color = 'black', linestyle = '-', marker='', label='nominal', linewidth=3)
            first = False
            nexxt = 2
          elif nexxt == 2:
            count+=1
            if count <= 2:
              None #voutPlt.plot(r['scale']['nominal'], r['vout']['nominal'], 
              #color = 'gray', linestyle = '--', marker='', label="failure_himp")
            else:
              None #voutPlt.plot(r['scale']['nominal'], r['vout']['nominal'], 
              #color = 'gray', linestyle = '--', marker='')
            nexxt = 3
          elif nexxt == 3:
            None #voutPlt.plot(r['scale']['nominal'], r['vout']['nominal'], 
            #color = 'gray', linestyle = ':', marker='')
            count+=1
            if count <= 2:
              None #voutPlt.plot(r['scale']['nominal'], r['vout']['nominal'], 
              #color = 'gray', linestyle = ':', marker='', label= 'failure_sck')
            else:
              None #voutPlt.plot(r['scale']['nominal'], r['vout']['nominal'], 
              #color = 'gray', linestyle = ':', marker='')
            nexxt = 2
          else:
            None #voutPlt.plot(r['scale']['nominal'], r['vout']['nominal'], 
            #color = 'gray', linestyle = '', marker='o', label='failure')

        errmax = np.amax(multidim, axis=0)
        errmin = np.amin(multidim, axis=0)
        voutPlt.fill_between(r['scale']['nominal'], errmin, errmax, hatch="||||", alpha=0.4, label="Failure range")

    voutPlt.grid(True)

    voutPlt.legend()

    name = datadirname + "/" + "diversityPlots/gen_" + str(generationNum) + "_duoplot_vout.eps" #ploting every generation in separate file
    plt.savefig(name)
    name = datadirname + "/" + "diversityPlots/gen_" + str(generationNum) + "_duoplot_vout.pdf" #ploting every generation in separate file
    plt.savefig(name)    
    #name = datadirname + "/" + "generationPlot" + ".png"
    #plt.savefig(name) 
    
    plt.close()



def duoplot_pareto(generation, generationNum, bestScoresList, result, bestI):
    """
    """

    fig, (pareto1, pareto2) = plt.subplots(nrows=1, ncols=2, sharex=True,
                                        figsize=(10, 4))

    #scoresPlt.set_title('Last generation\ncost values')
    pareto1.set_title('Last generation\nsum(F) vs. fnom objective')
    pareto1.plot(allObjectiveScores[:,0], allObjectiveScores[:,1], 'o', label='sigma(F)/f_{nom}', color="gray")
    pareto1.grid(True)

    pareto1.set_ylabel('sum(F)')
    pareto1.set_xlabel('fnom')

    pareto1.get_yaxis().set_major_formatter(matplotlib.ticker.FuncFormatter(lambda x, p: format(int(x), ',')))



    #mtrxDensPlt.set_title('Last generation\n matrix fillfactors')
    pareto2.set_title('Last generation\nstd(F) vs. fnom objective')
    pareto2.plot(allObjectiveScores[:,0], allObjectiveScores[:,2], 'o', label='FailuresSTD/Nominal', color="gray")
    pareto2.grid(True)

    pareto2.set_ylabel('std(F)')
    pareto2.set_xlabel('fnom')

    name = datadirname + "/" + "diversityPlots/gen_" + str(generationNum) + "_pareto.eps" 
    plt.savefig(name)
    name = datadirname + "/" + "diversityPlots/gen_" + str(generationNum) + "_pareto.pdf" 
    plt.savefig(name)

    plt.close()

#duoplot_sqrt(generation, generationNum, bestScoresList, result, bestI)
#duoplot_pareto(generation, generationNum, bestScoresList, result, bestI)
duoplot_abs(generation, generationNum, bestScoresList, result, bestI)
abscirc_MOEA(generation, generationNum, bestScoresList, result, bestI)