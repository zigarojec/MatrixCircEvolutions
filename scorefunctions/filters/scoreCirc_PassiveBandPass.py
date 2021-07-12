def scoreCirc_PassiveBandPass(circuit, gen, indi, MOEAMODE):#LGA+MOEA
  """
  Anchestor: scoreCirc_ActiveFilter_2
  This function takes:
    - 	a circuit object, 
    - 	generetion number,
    - 	individual number,
    -	0 - if score function is to be used in single-objective mode (score is a single float)
	1 - if score function is to be used in multi-objective mode (score is a vector of three objectives)
	
  This function returns:
    - 	score (float value (or array of floats) of the objective(s))
    -	results (dictionary from runme2 script)
  
  This function measures no THD beacuse of being passive circuit.  
  """  
  if debug > 2:
    print("\t\tG_" + str(gen) + "_I_" + str(indi))

  #---------------------------------------------------------BigMatrix stuff, check short-circuits, matrix density, matrix identifier (obsolete)
  FullBigCircuitMatrix = circuit.fullRedundancyMatrix
  OcSc, IcNc, SelfConnElm = checkConnsConnected(FullBigCircuitMatrix) #Outer connections Short cut, Inner connections Not connected
  #---------------------------------------------------------
  
  score = np.array([0.0,0.0,0.0]) if MOEAMODE == 1 else 0.0
  
  score += 2e4*np.exp(OcSc)
  results = None  
  if OcSc > 1:
    score += 1e4*np.exp(OcSc)
    # for the sake of problems in sorting when scores are equal, add a small random difference
    score += np.array([random.random()*10,random.random()*10,random.random()*10]) if MOEAMODE == 1 else random.random()*10
    
  else:
    #----------------------------------------------------------Try to make netlist and evaluate the individual
    makeNetlist(circuit, gen, indi, FullBigCircuitMatrix)
    results = runme2.evaluateActiveF_2_PB(gen, indi)   
    #----------------------------------------------------------Start of results analysis and objectives creation
    disfCount = 0
    
    ripple = np.array(results['ripple']['nominal'], dtype=float)#---------------ripple
    if np.isnan(ripple):
      disfCount = disfCount + 1
      r = 0      
    else:
      r = abs(ripple - 0.5) if ripple > 0.5 else 0      

    dampingL = np.array(results['damping_L']['nominal'], dtype=float)#----------damping
    if np.isnan(dampingL):
      disfCount = disfCount + 1
      dL = 0
    else:
      dL = abs(60 + dampingL) if dampingL > -60 else 0
      
    dampingH = np.array(results['damping_H']['nominal'], dtype=float)
    if np.isnan(dampingH):
      disfCount = disfCount + 1
      dH = 0
    else:
      dH = abs(60 + dampingH) if dampingH > -60 else 0      
       
    gain = np.array(results['gain']['nominal'], dtype=float)#-------------------gain
    if np.isnan(gain):
      disfCount = disfCount + 1
      g = 0
    else:
      g = abs(gain - 0) if gain < 0 else 0
	
    bandwidthL = np.array(results['bw_start']['nominal'], dtype=float)#---------bandwidth
    if np.isnan(bandwidthL):
      disfCount = disfCount + 0.2
      bwL = 0
    else:
      bwL = abs(np.log10(bandwidthL) - np.log10(abs(1000)))

    bandwidthH = np.array(results['bw_stop']['nominal'], dtype=float)
    if np.isnan(bandwidthH):
      disfCount = disfCount + 0.2
      bwH = 0
    else:
      bwH = abs(np.log10(bandwidthH) - np.log10(abs(30000)))      
      

    if debug > 2:      
      print(r, dL, dH, g, bwL, bwH)
    
    #----------------------------------------------------------Score function SINGLE-OBJECTIVE
    if MOEAMODE == 0:
      score =  1*r + (dH + dL) + 100*(bwL + bwH) + g*10

      if disfCount > 0:
        score += np.exp(disfCount) * 1e3 + random.random()*10
      
    #----------------------------------------------------------Score function MULTI-OBJECTIVE
    else: #MOEAMODE == 1
      score = np.array([(bwL + bwH), r, (dH + dL)]) + 10*g

      if disfCount > 0:
        score += (np.array([0,0,0])+ disfCount * 1e3) + np.array([random.random()*10,random.random()*10,random.random()*10])
      
      #print "score64", score
      if gen > 30:
        for i in range(0, len(score)): #round to 5 decimal points
            score[i] = np.float16(score[i]) 
        #print "score16", score
      score = np.ndarray.tolist(score) #for faster non-dominant sorting let the score be python array instead of np.array 5.7.2017
    #-------------------------------------------------------------------
    if debug > 2:    
      print("\t\tG_" + str(gen) + "_I_" + str(indi) + " SCORE:", score)
    
    filename = "g_" + str(gen) + "_i_" + str(indi) + "_subckt.cir"
    os.remove(filename) #cleanup current subcircuit
    
  return score, results  #removed matrixQuaziID and matrixDensity!!
