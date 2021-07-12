def scoreCirc_PassiveFilter_2(circuit, gen, indi, makeRedundancyInMatrix):
  """
  makeRedundancyInMatrix == False - if redundancy matrix already exists in circuit.BigCircuitMatrix
  makeRedundancyInMatrix == True  - if redundancy matrix has to be calculated.
  
  #TODO:
  Commentary: The score function ineed produces high slope, but rather a notch filter at 1kHz than LP filter!! FIX THAT!
  
  
  """
  #Calculate density and uniquiness (as in makeNetlist)
  if makeRedundancyInMatrix == True:
    FullBigCircuitMatrix = deepcopy(fullRedundancyBigCircuitMatrix(circuit.BigCircuitMatrix))
  else:
    FullBigCircuitMatrix = deepcopy(circuit.BigCircuitMatrix)

  rowsR,columnsR,columnsC,rowsC = sortedNonZeroIndices(FullBigCircuitMatrix)

  matrixDensity = float(len(rowsR))/float((BigMatrixSize*BigMatrixSize/2))	#(ones/(all/2))
  matrixQuaziID = sum(rowsR)+sum(columnsR)-BigMatrixSize*(BigMatrixSize-1)
  OcSc, IcNc, SelfConnElm = checkConnsConnected(FullBigCircuitMatrix) #Outer connections Short cut, Inner connections Not connected
  #print "Kratkih stikov zunanjih povezav:", OcSc
  
  results = None
  if OcSc > 1:
    score = 1e4*np.exp(OcSc)
  else:
    makeNetlist(circuit, gen, indi, FullBigCircuitMatrix)
    results = runme2.evaluatePassiveFilter_2(gen, indi)
    
    disfCount = 0

    gain = np.array(results['gain']['nominal'], dtype=float)
    if np.isnan(gain):
      disfCount = disfCount + 1
      g = 0
    else:
      g = abs(gain - 0) if gain < 0 else 0

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
      d = abs(40 - damping) if damping < 40 else 0

    slope = np.array(results['dumpingSlope']['nominal'], dtype=float)
    if np.isnan(slope):
      disfCount = disfCount + 1
      slo = 0
    else:
      slo = 0 if slope>60 else 60-slope
    
    bandwidth = np.array(results['bw']['nominal'], dtype=float)
    if np.isnan(bandwidth):
      disfCount = disfCount + 1
      bw = 0
    else:
      bw = abs(bandwidth-1000)/100
      
    #THD = np.array(results['THD']['nominal'], dtype=float)
    #if np.isnan(THD):
    #  disfCount = disfCount + 1
    #  thd = 0
    #else:
    #  thd = THD-1 if THD > 1 else 0
    #print 10*r, g, d, slo, bw
    score = 10*r + g + d + slo + bw

    if disfCount > 0:
      score += np.exp(disfCount) * 1e3

    #score = score + (IcNc*IcNc+1)# + abs(BW-bw)*1e2 + abs(CUTOFF-cutoff)*1e2 #add small punishment if not all nodes connected and bw and cutoff are off

    #print "\t\t\t\t\tG_" + str(gen) + "_I_" + str(indi) + " SCORE:", score
    #cleanup current subcircuit
    filename = "g_" + str(gen) + "_i_" + str(indi) + "_subckt.cir"
    os.remove(filename)
  return score, matrixDensity, matrixQuaziID, results


