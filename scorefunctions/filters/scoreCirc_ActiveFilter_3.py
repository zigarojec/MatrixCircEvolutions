def scoreCirc_ActiveFilter_3(circuit, gen, indi, makeRedundancyInMatrix):
  """
  makeRedundancyInMatrix == False - if redundancy matrix already exists in circuit.BigCircuitMatrix
  makeRedundancyInMatrix == True  - if redundancy matrix has to be calculated.
  """

  FullBigCircuitMatrix = deepcopy(circuit.fullRedundancyMatrix)

  rowsR,columnsR,columnsC,rowsC = sortedNonZeroIndices(FullBigCircuitMatrix)

  matrixDensity = float(len(rowsR))/float((BigMatrixSize*BigMatrixSize/2))	#(ones/(all/2))
  matrixQuaziID = sum(rowsR)+sum(columnsR)-BigMatrixSize*(BigMatrixSize-1)
  OcSc, IcNc, SelfConnElm = checkConnsConnected(FullBigCircuitMatrix) #Outer connections Short cut, Inner connections Not connected
  #print "Kratkih stikov zunanjih povezav:", OcSc
  
  score = 0
  results = None
  if OcSc > 1:
    score = 1e4*np.exp(OcSc)
  else:
    makeNetlist(circuit, gen, indi, FullBigCircuitMatrix)
    results = runme2.evaluateActiveFilter_2(gen, indi)

    disfCount = 0
    
    ripple = np.array(results['ripple']['nominal'], dtype=float)
    if np.isnan(ripple):
      disfCount = disfCount + 1
      r = 0      
    else:
      r = abs(ripple - 0.5) if ripple > 0.5 else 0
      
    damping = np.array(results['damping']['nominal'], dtype=float)
    if np.isnan(damping):
      disfCount = disfCount + 1
      d = 0
    else:
      d = abs(20 - damping) if damping < 20 else 0
      
    gain = np.array(results['gain']['nominal'], dtype=float)
    if np.isnan(gain):
      disfCount = disfCount + 1
      g = 0
    else:
      g = abs(gain - 0)# if gain < 10 else 0.01
      
    THD_Lf = np.array(results['THD_Lf']['nominal'], dtype=float)
    if np.isnan(THD_Lf):
      disfCount = disfCount + 1
      thd_lf = 0
    else:
      thd_lf = THD_Lf-1 if THD_Lf > 1 else 0
      
    THD_Hf = np.array(results['THD_Hf']['nominal'], dtype=float)
    if np.isnan(THD_Hf):
      disfCount = disfCount + 1
      thd_hf = 0
    else:
      thd_hf = THD_Hf-1 if THD_Hf > 1 else 0
      
    #RIN = np.array(results['rin_meas']['nominal'], dtype=float) #--------not in use
    #if np.isnan(RIN):
    #  disfCount = disfCount + 1
    #  rin = 0
    #else:
    #  rin = 1/RIN*1e6 if RIN < 1e7 else 0

    isLP = np.array(results['is_LP']['nominal'], dtype=float)
    if np.isnan(isLP):
      disfCount = disfCount + 1
      islp = 0
    else:
      islp = 0 if isLP>0 else 100# np.abs(isLP)
      
    #slope = np.array(results['maxDampingSlope']['nominal'], dtype=float)
    #print slope
    #if np.isnan(slope):
    #  disfCount = disfCount + 1
    #  slo = 0
    #else:
    #  slo = 0 if slope>60 else 60-slope
    
    maxSlope = results['maxDampingSlope']['nominal']
    if type(np.nan) == type(maxSlope) or type(None) == type(maxSlope):
      disfCount = disfCount + 2
      slo = 0
      slof = 0 
    else:
      if len(maxSlope)==2:
        slo = 0 if maxSlope[0]>60 else 60-maxSlope[0]
        slof = np.log10(abs(maxSlope[1]-1000))
      else:
        slo = 0
        slof = 0
        disfCount = disfCount + 1    
    
    
    bandwidth = np.array(results['bw']['nominal'], dtype=float)
    if np.isnan(bandwidth):
      #disfCount = disfCount + 1
      bandwidth = 0
    bw = abs(bandwidth-1000)
    
    StaticOut = not results['isOutVNonStationary']['nominal']
    score = 10*slo + 10*r + (100*StaticOut + 10*(thd_lf + thd_hf) + 1*islp + g)#rin!

    #print disfCount
    if disfCount > 0:
      score = 0 + np.exp(disfCount) * 1e3
      #print "disfCount was there"

    #score = score + (IcNc+1)# + abs(BW-bw)*1e2 + abs(CUTOFF-cutoff)*1e2 #add small punishment if not all nodes connected and bw and cutoff are off

      
    print("\t\t\t\tG_" + str(gen) + "_I_" + str(indi) + " SCORE:", score)
    #cleanup current subcircuit
    filename = "g_" + str(gen) + "_i_" + str(indi) + "_subckt.cir"
    os.remove(filename)
    #print ".",
    #circuit.objectivesScore = copy(score)	#id does not work with mpirun since mpirun works with copies
    #circuit.matrixDensity = matrixDensity
  return score, matrixDensity, matrixQuaziID, results

