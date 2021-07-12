
#------------------------------------------------------------
#Almost obsolete stuff:--------------------------------------
#Functions below are left there for some history learning.---
#------------------------------------------------------------
def scoreCirc_LP_PassiveFilter(circuit, gen, indi):
  """To write."""
  #Calculate density and uniquiness (as in makeNetlist)
  FullBigCircuitMatrix = deepcopy(fullRedundancyBigCircuitMatrix(circuit.BigCircuitMatrix))

  rowsR,columnsR,columnsC,rowsC = sortedNonZeroIndices(FullBigCircuitMatrix)

  matrixDensity = float(len(rowsR))/float((BigMatrixSize*BigMatrixSize/2))	#(ones/(all/2))
  matrixQuaziID = sum(rowsR)+sum(columnsR)-BigMatrixSize*(BigMatrixSize-1)

  OcSc, IcNc, SelfConnElm = checkConnsConnected(FullBigCircuitMatrix) #Outer connections Short cut, Inner connections Not connected
  #print "Kratkih stikov zunanjih povezav:", OcSc

  results = None
  if OcSc > 1:
    score = 1e5*np.exp(OcSc)
  else:
    makeNetlist(circuit, gen, indi, FullBigCircuitMatrix)
    results = runme2.evaluatePassiveFilter(gen, indi)		#Evaluate a circuit (diff stage) and gather results
    
    
    disfCount = 0
    tf = results['tf']['default']
    if tf != None:
      freq_scale = results['freq_scale']['default']
      idealFilter = results['idealFilter']
      maxOfPartPass, minOfPartPass, deltaPass, deltaMeanRefPass, gradientPass = minmaxMeas(tf, freq_scale, 1, BW) #passband
      maxOfPartStop, minOfPartStop, deltaStop, deltaMeanRefStop, gradientStop = minmaxMeas(tf, freq_scale, CUTOFF, 1e4 - 1)#stopband
      maxOfPartTrans, x, x, x, gradientTrans = minmaxMeas(tf, freq_scale, BW, CUTOFF)#transition
      """  
      #evklidova razdalja med tf in ideal_tf
      #score = np.linalg.norm(idealFilter-tf)
      order = 2
      score = sum(np.abs((idealFilter-tf)*10)**order)**(1./order)	#razdalja 3. reda bo mogoce zgladila spice, ki se pojavijo
      """
      score =  (-1)*np.divide(np.arctan(SelfConnElm),10.0) + deltaMeanRefPass*1e3 + abs(PARTPASS-deltaMeanRefStop)*2e2 + abs(maxOfPartTrans-1)*1e2 + deltaPass*1e1 + deltaStop*1e2 + gradientPass*1e1 + gradientStop*1e1
      
      
      bw = results['bw']['default']
      if bw == None:
        disfCount = disfCount + 1
      cutoff = results['cutoff']['default']
      if cutoff == None:
        disfCount = disfCount + 1  
    else:
      disfCount = disfCount + 1
      score = 1e4
    
    
    
    #score =	abs(BW - float(bw or 0))/1000 + 	\	
    #	abs(CUTOFF - float(cutoff or 0))/1000
    
    if disfCount > 0:
      score = (np.exp(disfCount)) * 1e3
    
    try:
      score = score *(IcNc*IcNc+1) #add small punishment if not all nodes connected
    except:
      score = score *(IcNc*IcNc+1)
    
    #print "\t\t\t\t\tG_" + str(gen) + "_I_" + str(indi) + " SCORE:", score
    #cleanup current subcircuit
    filename = "G_" + str(gen) + "_I_" + str(indi) + "_subckt.cir"
    os.remove(filename)
  """
  try:
    print "Is results."
    return score, matrixDensity, matrixQuaziID, results
    print "Is results drugic."
  except:
    print "Is not results"
    return score, matrixDensity, matrixQuaziID, None
  """
  return score, matrixDensity, matrixQuaziID, results

def scoreCirc_PassiveFilterOLDER(circuit, gen, indi):
  """To write."""
  #Calculate density and uniquiness (as in makeNetlist)
  FullBigCircuitMatrix = deepcopy(fullRedundancyBigCircuitMatrix(circuit.BigCircuitMatrix))

  rowsR,columnsR,columnsC,rowsC = sortedNonZeroIndices(FullBigCircuitMatrix)

  matrixDensity = float(len(rowsR))/float((BigMatrixSize*BigMatrixSize/2))	#(ones/(all/2))
  matrixQuaziID = sum(rowsR)+sum(columnsR)-BigMatrixSize*(BigMatrixSize-1)

  OcSc, IcNc = checkConnsConnected(FullBigCircuitMatrix) #Outer connections Short cut, Inner connections Not connected
  #print "Kratkih stikov zunanjih povezav:", OcSc

  results = None
  if OcSc > 1:
    score = 1e5*np.exp(OcSc)
  else:
    makeNetlist(circuit, gen, indi, FullBigCircuitMatrix)
    results = runme2.evaluatePassiveFilter(gen, indi)		#Evaluate a circuit (diff stage) and gather results
    
    
    disfCount = 0
    tf = results['tf']['default']
    if tf != None:
      freq_scale = results['freq_scale']['default']
      idealFilter = results['idealFilter']
      
      #evklidova razdalja med tf in ideal_tf
      #score = np.linalg.norm(idealFilter-tf)
      order = 2
      score = sum(np.abs((idealFilter-tf)*10)**order)**(1./order)	#razdalja 3. reda bo mogoce zgladila spice, ki se pojavijo

      
      bw = results['bw']['default']
      if bw == None:
        disfCount = disfCount + 1
      cutoff = results['cutoff']['default']
      if cutoff == None:
        disfCount = disfCount + 1  
    else:
      disfCount = disfCount + 1
      score = 1e4
    
    
    
    #score =	abs(BW - float(bw or 0))/1000 + 	\	
    #	abs(CUTOFF - float(cutoff or 0))/1000
    
    if disfCount > 0:
      score = (np.exp(disfCount)) * 1e3
    
    try:
      score = score *(IcNc*IcNc+1) + abs(BW-bw)/float(1000) + abs(CUTOFF-cutoff)/float(1000) #add small punishment if not all nodes are connected and bw and cutoff are off
    except:
      score = score *(IcNc*IcNc+1)
    
    #print "\t\t\t\t\tG_" + str(gen) + "_I_" + str(indi) + " SCORE:", score
    #cleanup current subcircuit
    filename = "G_" + str(gen) + "_I_" + str(indi) + "_subckt.cir"
    os.remove(filename)
  """
  try:
    print "Is results."
    return score, matrixDensity, matrixQuaziID, results
    print "Is results drugic."
  except:
    print "Is not results"
    return score, matrixDensity, matrixQuaziID, None
  """
  return score, matrixDensity, matrixQuaziID, results

def scoreCirc_PassiveFilterOLD(circuit, gen, indi, makeRedundancyInMatrix):
  """To write."""
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
    score = 1e5*np.exp(OcSc)
  else:
    makeNetlist(circuit, gen, indi, FullBigCircuitMatrix)
    results = runme2.evaluatePassiveFilter(gen, indi)		#Evaluate a circuit (diff stage) and gather results
    
    
    disfCount = 0
    tf = results['tf']['default']
    if tf != None:
      freq_scale = results['freq_scale']['default']
      idealFilter = results['idealFilter']
      #Passive low pass filter:
      maxOfPartPass, minOfPartPass, deltaPass, deltaMeanRefPass, gradientPass = minmaxMeas(tf, freq_scale, 1, BW)#1000, 1e4 - 1) #passband
      maxOfPartStop, minOfPartStop, deltaStop, deltaMeanRefStop, gradientStop = minmaxMeas(tf, freq_scale, CUTOFF, 1e4-1)#1, 800)#stopband
      maxOfPartTrans, x, x, x, gradientTrans = minmaxMeas(tf, freq_scale, BW, CUTOFF)#transition
      #Cost function:
      score =  (-1)*np.divide(np.arctan(SelfConnElm),10.0) + deltaMeanRefPass*1e3 + abs(PARTPASS-deltaMeanRefStop)*2e2 + abs(maxOfPartTrans-1)*1e2 + deltaPass*1e1 + deltaStop*1e1 + gradientPass*1e1 + gradientStop*1e1
      
      ##Passive notch filter:
      #maxOfPartPassL, minOfPartPassL, deltaPassL, deltaMeanRefPassL, gradientPassL = minmaxMeas(tf, freq_scale, 1, 45)#1000, 1e4 - 1) #passbandL
      #maxOfPartPassH, minOfPartPassH, deltaPassH, deltaMeanRefPassH, gradientPassH = minmaxMeas(tf, freq_scale, 55, 1e5-1)#1000, 1e4 - 1) #passbandH
      #maxOfPartStop, minOfPartStop, deltaStop, deltaMeanRefStop, gradientStop = minmaxMeas(tf, freq_scale, 45, 55)#1, 800)#stopband
      ##Cost function:
      #score =  (-1)*np.divide(np.arctan(SelfConnElm),10.0) + deltaMeanRefPassL*1e3 + deltaMeanRefPassH*1e3 + deltaPassL*1e2 + deltaPassH*1e2 + abs(maxOfPartStop)*1e3 + abs(1-deltaStop)*1e3
      
      
      """  
      #evklidova razdalja med tf in ideal_tf
      #score = np.linalg.norm(idealFilter-tf)
      order = 2
      score = sum(np.abs((idealFilter-tf)*10)**order)**(1./order)	#razdalja 3. reda bo mogoce zgladila spice, ki se pojavijo
      """

      
      #bw = results['bw']['default']
      #if bw == None:
        #disfCount = disfCount + 1
      #cutoff = results['cutoff']['default']
      #if cutoff == None:
        #disfCount = disfCount + 1  
    else:
      disfCount = disfCount + 1
      score = 1e4
    
    
    
    #score =	abs(BW - float(bw or 0))/1000 + 	\	
    #	abs(CUTOFF - float(cutoff or 0))/1000
    #print disfCount
    #if disfCount > 0:
    #  score = (np.exp(disfCount)) * 1e3
    
    try:
      score = score *(IcNc*IcNc+1)# + abs(BW-bw)*1e2 + abs(CUTOFF-cutoff)*1e2 #add small punishment if not all nodes connected and bw and cutoff are off
    except:
      score = score *(IcNc*IcNc+1)
      
    #print "\t\t\t\t\tG_" + str(gen) + "_I_" + str(indi) + " SCORE:", score
    #cleanup current subcircuit
    filename = "G_" + str(gen) + "_I_" + str(indi) + "_subckt.cir"
    os.remove(filename)
  """
  try:
    print "Is results."
    return score, matrixDensity, matrixQuaziID, results
    print "Is results drugic."
  except:
    print "Is not results"
    return score, matrixDensity, matrixQuaziID, None
  """
  return score, matrixDensity, matrixQuaziID, results

def scoreCirc_ActiveFilter(circuit, gen, indi, makeRedundancyInMatrix):#TODO
  """
  makeRedundancyInMatrix == False - if redundancy matrix already exists in circuit.BigCircuitMatrix
  makeRedundancyInMatrix == True  - if redundancy matrix has to be calculated.
  """
  #Calculate density and uniquiness (as in makeNetlist)
  if makeRedundancyInMatrix == True:
    #FullBigCircuitMatrix = deepcopy(fullRedundancyBigCircuitMatrix(circuit.BigCircuitMatrix))
    FullBigCircuitMatrix = deepcopy(circuit.fullRedundancyMatrix)
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
    results = runme2.evaluateActiveFilter_SUHAD(gen, indi)#TODO
    
    
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
      d = abs(40 - damping) if damping < 40 else 0
      
    gain = np.array(results['gain']['nominal'], dtype=float)
    if np.isnan(gain):
      disfCount = disfCount + 1
      g = 0
    else:
      g = abs(gain - 10) if gain < 10 else 0
      
    THD = np.array(results['THD']['nominal'], dtype=float)
    if np.isnan(THD):
      disfCount = disfCount + 1
      thd = 0
    else:
      thd = THD-1 if THD > 1 else 0
	  
    StaticOut = not results['isOutVNonStationary']['nominal']
      
    score = 5*r + 4*d + 2*g + (100*StaticOut + 10*thd)

    #print disfCount
    if disfCount > 0:
      score = np.exp(disfCount) * 1e3
    
    ##add a little salt!
    #score = score + random.uniform(0.0, 1)

    score = score + (IcNc*IcNc+1)# + abs(BW-bw)*1e2 + abs(CUTOFF-cutoff)*1e2 #add small punishment if not all nodes connected and bw and cutoff are off

      
    #print "\t\t\t\t\tG_" + str(gen) + "_I_" + str(indi) + " SCORE:", score
    #cleanup current subcircuit
    filename = "g_" + str(gen) + "_i_" + str(indi) + "_subckt.cir"
    os.remove(filename)
    #print ".",
  return score, matrixDensity, matrixQuaziID, results

def scoreCirc_PassiveFilter(circuit, gen, indi, makeRedundancyInMatrix):#TODO
  """
  makeRedundancyInMatrix == False - if redundancy matrix already exists in circuit.BigCircuitMatrix
  makeRedundancyInMatrix == True  - if redundancy matrix has to be calculated.
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
    results = runme2.evaluatePassiveFilter_SUHAD(gen, indi)#TODO
    
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
      r = abs(ripple - 0.5)# if ripple > 0.5 else 0
      
    damping = np.array(results['damping']['nominal'], dtype=float)
    if np.isnan(damping):
      disfCount = disfCount + 1
      d = 0
    else:
      d = abs(40 - damping)# if damping < 60 else 0
      
    #THD = np.array(results['THD']['nominal'], dtype=float)
    #if np.isnan(THD):
    #  disfCount = disfCount + 1
    #  thd = 0
    #else:
    #  thd = THD-1 if THD > 1 else 0
    
    score = 10*r + g + 10*d

    if disfCount > 0:
      score = np.exp(disfCount) * 1e3

    score = score + (IcNc*IcNc+1)# + abs(BW-bw)*1e2 + abs(CUTOFF-cutoff)*1e2 #add small punishment if not all nodes connected and bw and cutoff are off

    #print "\t\t\t\t\tG_" + str(gen) + "_I_" + str(indi) + " SCORE:", score
    #cleanup current subcircuit
    filename = "g_" + str(gen) + "_i_" + str(indi) + "_subckt.cir"
    os.remove(filename)
  return score, matrixDensity, matrixQuaziID, results

def scoreCirc_VoltageReference(circuit, gen, indi, makeRedundancyInMatrix):
  """Calculates the score of the voltage reference circuit given in argument. """
  #----------#
  VREF = 1.5
  #----------#
  
  FullBigCircuitMatrix = deepcopy(circuit.fullRedundancyMatrix)
  rowsR,columnsR,columnsC,rowsC = sortedNonZeroIndices(FullBigCircuitMatrix)

  matrixDensity = float(len(rowsR))/float((BigMatrixSize*BigMatrixSize/2))	#(ones/(all/2))
  matrixQuaziID = sum(rowsR)+sum(columnsR)-BigMatrixSize*(BigMatrixSize-1)
  OcSc, IcNc, SelfConnElm = checkConnsConnected(FullBigCircuitMatrix) #Outer connections Short cut, Inner connections Not connected
  
  results = None
  badSweep = 0
  if OcSc > 1:
    score = 1e4*np.exp(OcSc)
  else:
    makeNetlist(circuit, gen, indi, FullBigCircuitMatrix)
    results = runme2.evaluateVoltageRef(gen, indi)
    disfCount = 0
    
    vdd_sweep = np.array(results['vout_vdd']['nominal'], dtype=float) #This line changes Nones to np.nans
    vdd_sweep_scale = np.array(results['vout_vdd_scale']['nominal'], dtype=float)
    # if measurement is empty 
    if np.any(np.isnan(vdd_sweep)):
      disfCount = disfCount + 1
      vdd_s = 0
      vdd_s_d = 0
      #print "tukej!", vdd_sweep_scale
    else:
      x = np.median(vdd_sweep)
      vdd_s =  abs(x - VREF) #if x > VREF else 0
      vdd_s_d = np.max(vdd_sweep) - np.min(vdd_sweep)
      #if sweep did not finish completely - add to score
      #check last scale value in runme2!!
      #print "tukiii", vdd_sweep_scale
      if (vdd_sweep_scale[-1]<20): #20V
        badSweep = badSweep + 1
     
    rload_sweep = np.array(results['vout_rload']['nominal'], dtype=float)
    rload_sweep_scale = np.array(results['vout_rload_scale']['nominal'], dtype=float)
    # if measurement is empty
    if np.any(np.isnan(rload_sweep)):
      disfCount = disfCount + 1
      rload_s = 0
      rload_s_d = 0
    else:
      x = np.median(rload_sweep)
      rload_s =   abs(x - VREF) #if x > VREF else 0
      rload_s_d = np.max(rload_sweep) - np.min(rload_sweep)
      #if sweep did not finish completely - add to score
      #check last scale value in runme2!!
      if (rload_sweep_scale[-1]<100e3): #100kOhm
        badSweep = badSweep + 1
      
    temp_sweep = np.array(results['vout_temp']['nominal'], dtype=float)
    temp_sweep_scale = np.array(results['vout_temp_scale']['nominal'], dtype=float)
    # if measurement is empty OR sweep did not finish completely - check last scale value in runme2!!
    if np.any(np.isnan(temp_sweep)):
      disfCount = disfCount + 1
      temp_s = 0
      temp_s_d = 0
    else:
      x = np.median(temp_sweep)
      temp_s =  abs(x - VREF) #if x > VREF else 0
      temp_s_d = np.max(temp_sweep) - np.min(temp_sweep)
      if (temp_sweep_scale[-1]<120): #120 deg celsius
        badSweep = badSweep + 1
      
    power = results['power']['nominal']
    if np.isnan(np.array(power, dtype=float)):
      disfCount = disfCount + 1
      powe = 0
    else:
      powe = power
    
    #---COST FUNCTION DEFINITION---#
    score = (vdd_s) + (vdd_s_d) + 5*(rload_s) + 5*(rload_s_d) + (100*temp_s) + (100*temp_s_d) + (100*powe) + badSweep*100

    #print disfCount
    if disfCount > 0:
      score = np.exp(disfCount) * 1e3
    if np.isnan(score):
      score = 2e4
    score = score + (IcNc+1) #add small punishment if not all nodes connected

    #print "\t\t\t\t\tG_" + str(gen) + "_I_" + str(indi) + " SCORE:", score
    #print vdd_s, vdd_s_d, rload_s, rload_s_d, temp_s, temp_s_d, powe
    #print vdd_s, vdd_s_d, rload_s, rload_s_d, 100*temp_s, 100*temp_s_d, 100*powe
    
    filename = "g_" + str(gen) + "_i_" + str(indi) + "_subckt.cir"
    os.remove(filename) #cleanup current subcircuit

  return score, matrixDensity, matrixQuaziID, results
  
def scoreCirc_CmosVoltageReference(circuit, gen, indi, makeRedundancyInMatrix): #TODO 6.9.2016 napisi cost function ki se sklada z evaluateCmosVoltageRef
  """Calculates the score of the voltage reference circuit given in argument. """
  #----------#
  VREF = 1.5
  #----------#
  
  FullBigCircuitMatrix = deepcopy(circuit.fullRedundancyMatrix)
  rowsR,columnsR,columnsC,rowsC = sortedNonZeroIndices(FullBigCircuitMatrix)

  matrixDensity = float(len(rowsR))/float((BigMatrixSize*BigMatrixSize/2))	#(ones/(all/2))
  matrixQuaziID = sum(rowsR)+sum(columnsR)-BigMatrixSize*(BigMatrixSize-1)
  OcSc, IcNc, SelfConnElm = checkConnsConnected(FullBigCircuitMatrix) #Outer connections Short cut, Inner connections Not connected
  
  results = None
  if OcSc > 1:
    score = 1e4*np.exp(OcSc)
  else:
    makeNetlist(circuit, gen, indi, FullBigCircuitMatrix)
    results = runme2.evaluateCmosVoltageRef(gen, indi)
    disfCount = 0
    
    
    #Vdd sweeps on 3 temperatures - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    # -20 deg
    vdd_sweep_scale = np.array(results['vout_vdd_scale']['nominal'], dtype=float)
    vdd_sweep_t1 = np.array(results['vout_vdd_temp1']['nominal'], dtype=float) #This line changes Nones to np.nans
    # if measurement is empty 
    if np.any(np.isnan(vdd_sweep_t1)):
      disfCount = disfCount + 1
      vdd_s_t1 = 0
      vdd_s_t1_d = 0
    else:
      x = np.median(vdd_sweep_t1)
      vdd_s_t1 =  abs(x - VREF) #if x > VREF else 0
      vdd_s_t1_d = np.max(vdd_sweep_t1) - np.min(vdd_sweep_t1)
    
    
    # 25 deg
    vdd_sweep_scale = np.array(results['vout_vdd_scale']['nominal'], dtype=float)
    vdd_sweep_t2 = np.array(results['vout_vdd_temp2']['nominal'], dtype=float) #This line changes Nones to np.nans
    # if measurement is empty 
    if np.any(np.isnan(vdd_sweep_t2)):
      disfCount = disfCount + 1
      vdd_s_t2 = 0
      vdd_s_t2_d = 0
    else:
      x = np.median(vdd_sweep_t2)
      vdd_s_t2 =  abs(x - VREF) #if x > VREF else 0
      vdd_s_t2_d = np.max(vdd_sweep_t2) - np.min(vdd_sweep_t2)    
    
    # 120 deg
    vdd_sweep_scale = np.array(results['vout_vdd_scale']['nominal'], dtype=float)
    vdd_sweep_t3 = np.array(results['vout_vdd_temp3']['nominal'], dtype=float) #This line changes Nones to np.nans
    # if measurement is empty 
    if np.any(np.isnan(vdd_sweep_t3)):
      disfCount = disfCount + 1
      vdd_s_t3 = 0
      vdd_s_t3_d = 0
    else:
      x = np.median(vdd_sweep_t3)
      vdd_s_t3 =  abs(x - VREF) #if x > VREF else 0
      vdd_s_t3_d = np.max(vdd_sweep_t3) - np.min(vdd_sweep_t3)     
    
    #Vdd sweeps on 3 loads - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    # 10e6 Ohm
    vdd_sweep_scale = np.array(results['vout_vdd_res_scale']['nominal'], dtype=float)
    vdd_sweep_r1 = np.array(results['vout_vdd_res1']['nominal'], dtype=float) #This line changes Nones to np.nans
    # if measurement is empty 
    if np.any(np.isnan(vdd_sweep_r1)):
      disfCount = disfCount + 1
      vdd_s_r1 = 0
      vdd_s_r1_d = 0
    else:
      x = np.median(vdd_sweep_r1)
      vdd_s_r1 =  abs(x - VREF) #if x > VREF else 0
      vdd_s_r1_d = np.max(vdd_sweep_r1) - np.min(vdd_sweep_r1)
     
    # 10e4 Ohm
    vdd_sweep_scale = np.array(results['vout_vdd_res_scale']['nominal'], dtype=float)
    vdd_sweep_r2 = np.array(results['vout_vdd_res2']['nominal'], dtype=float) #This line changes Nones to np.nans
    # if measurement is empty 
    if np.any(np.isnan(vdd_sweep_r2)):
      disfCount = disfCount + 1
      vdd_s_r2 = 0
      vdd_s_r2_d = 0
    else:
      x = np.median(vdd_sweep_r2)
      vdd_s_r2 =  abs(x - VREF) #if x > VREF else 0
      vdd_s_r2_d = np.max(vdd_sweep_r2) - np.min(vdd_sweep_r2)    
    
    # 10e2 Ohm
    vdd_sweep_scale = np.array(results['vout_vdd_res_scale']['nominal'], dtype=float)
    vdd_sweep_r3 = np.array(results['vout_vdd_res3']['nominal'], dtype=float) #This line changes Nones to np.nans
    # if measurement is empty 
    if np.any(np.isnan(vdd_sweep_r3)):
      disfCount = disfCount + 1
      vdd_s_r3 = 0
      vdd_s_r3_d = 0
    else:
      x = np.median(vdd_sweep_r3)
      vdd_s_r3 =  abs(x - VREF) #if x > VREF else 0
      vdd_s_r3_d = np.max(vdd_sweep_r3) - np.min(vdd_sweep_r3)     
      
    power = results['power']['nominal']
    if np.isnan(np.array(power, dtype=float)):
      disfCount = disfCount + 1
      powe = 0
    else:
      powe = power
    
    #---COST FUNCTION DEFINITION---#
    score =    vdd_s_t1 + vdd_s_t1_d + \
	       vdd_s_t2 + vdd_s_t2_d + \
	       vdd_s_t3 + vdd_s_t3_d + \
	       vdd_s_r1 + vdd_s_r1_d + \
	    vdd_s_r2 + vdd_s_r2_d + \
	    vdd_s_r3 + vdd_s_r3_d + \
	    (100*powe)

    #print disfCount
    if disfCount > 0:
      score = np.exp(disfCount) * 1e3
    if np.isnan(score):
      score = 2e4
    score = score + (IcNc+1) #add small punishment if not all nodes connected

    #print "\t\t\t\t\tG_" + str(gen) + "_I_" + str(indi) + " SCORE:", score
    
    filename = "g_" + str(gen) + "_i_" + str(indi) + "_subckt.cir"
    os.remove(filename) #cleanup current subcircuit

  return score, matrixDensity, matrixQuaziID, results

