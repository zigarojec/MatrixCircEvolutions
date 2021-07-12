def scoreCirc_HighPass(circuit, gen, indi, MOEAMODE):#LGA+MOEA
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
    results = runme2.evaluateActiveF_2_HP(gen, indi)   
    #----------------------------------------------------------Start of results analysis and objectives creation
    disfCount = 0
    
    ripple = np.array(results['ripple']['nominal'], dtype=float)#---------------ripple
    if np.isnan(ripple):
      disfCount = disfCount + 1
      r = 0      
    else:
      r = abs(ripple - 1) if ripple > 1 else 0        
       
    gain = np.array(results['gain']['nominal'], dtype=float)#-------------------gain
    if np.isnan(gain):
      disfCount = disfCount + 1
      g = 0
      gain = -99
    else:
      g = abs(gain - 0) if gain < 0 else 0
      
    dampingL = np.array(results['damping_L']['nominal'], dtype=float)#----------damping
    if np.isnan(dampingL):
      disfCount = disfCount + 1
      dL = 0
    else:
      dL = abs(60 + dampingL - gain) if (dampingL-gain) > -60 else 0  
      
    bandwidthL = np.array(results['bw_start']['nominal'], dtype=float)#---------bandwidth
    if np.isnan(bandwidthL):
      disfCount = disfCount + 0.2
      bwL = 0
    else:
      bwL = abs(np.log10(bandwidthL) - np.log10(abs(8000)))
      if bandwidthL < 200:
        bwL = 1e6*(1/np.exp(bandwidthL/20))
      

    if debug > 2:      
      print(r, dL, g, bwL)
    
    #----------------------------------------------------------Score function SINGLE-OBJECTIVE
    if MOEAMODE == 0:
      #score =  7*r + 20*(dL) + 50*(bwL) + 20*g
      score =  15*r + 10*(dL) + 5*(bwL) + 4*g

      if disfCount > 0:
        score += np.exp(disfCount) * 1e3 + random.random()*10
	
    #-EXPERIMENTAL---------------------------------------------Score function for PSADE
    # NO !!! NO GOOD!
    elif MOEAMODE == 2:
      score =  5*r + 3*(dL) + 200*(bwL) + 3*g

      if disfCount > 0:
        score += np.exp(disfCount) * 1e3 + random.random()*10      
      
    #----------------------------------------------------------Score function MULTI-OBJECTIVE
    elif MOEAMODE == 1:
      score = np.array([(bwL), r, (dL)]) + 10*g

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
