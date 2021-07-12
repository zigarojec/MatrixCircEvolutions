

#Hybrid LGA+MOEA score functions

def scoreCirc_ActiveFilter_2(circuit, gen, indi, MOEAMODE):#LGA+MOEA
  """
  Anchestor: scoreCirc_ActiveFilter_MOEA
  This function takes:
    -   a circuit object, 
    -   generetion number,
    -   individual number,
    -   0 - if score function is to be used in single-objective mode (score is a single float)
        1 - if score function is to be used in multi-objective mode (score is a vector of three objectives)
    
  This function returns:
    - 	score (float value (or array of floats) of the objective(s))
    -	results (dictionary from runme2 script)
  
  This function measures THD at low freq default, and at bw that is measured via evaluation function.  
  """  
  if debug > 2:
    print("\t\tG_" + str(gen) + "_I_" + str(indi))

  #---------------------------------------------------------BigMatrix stuff, check short-circuits, matrix density, matrix identifier (obsolete)
  FullBigCircuitMatrix = circuit.fullRedundancyMatrix
  OcSc, IcNc, SelfConnElm = checkConnsConnected(FullBigCircuitMatrix) #Outer connections Short cut, Inner connections Not connected
  #---------------------------------------------------------
  
  score = np.array([0,0,0], dtype="float64") if MOEAMODE == 1 else 0
  
  score += 2e4*np.exp(OcSc)
  results = None  
  if OcSc > 1:
    score += 1e4*np.exp(OcSc)
    # for the sake of problems in sorting when scores are equal, add a small random difference
    score += np.array([random.random()*10,random.random()*10,random.random()*10]) if MOEAMODE == 1 else random.random()*10
    
  else:
    #----------------------------------------------------------Try to make netlist and evaluate the individual
    makeNetlist(circuit, gen, indi, FullBigCircuitMatrix)
    results = runme2.evaluateActiveFilter_2(gen, indi)   
    #----------------------------------------------------------Start of results analysis and objectives creation
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
      d = abs(60 - damping) if damping < 60 else 0
      
    gain = np.array(results['gain']['nominal'], dtype=float)
    if np.isnan(gain):
      disfCount = disfCount + 1
      g = 0
    else:
      #g = abs(gain - 0) if gain < 0 else 0.01
      g = abs(gain - 10)

    THD_Lf = np.array(results['THD_Lf']['nominal'], dtype=float)
    if np.isnan(THD_Lf):
      disfCount = disfCount + 1
      thd_lf = 0
    else:
      thd_lf = THD_Lf-1 if THD_Lf > 1 else 0	

    #THD_Hf = np.array(results['THD_Hf']['nominal'], dtype=float)
    #if np.isnan(THD_Hf):
      #disfCount = disfCount + 1
      #thd_hf = 0
    #else:
      #thd_hf = THD_Hf-1 if THD_Hf > 1 else 0
      
    isLP = np.array(results['is_LP']['nominal'], dtype=float)
    if np.isnan(isLP):
      disfCount = disfCount + 1
      islp = 0
    else:
      islp = 0 if isLP>0 else 100# np.abs(isLP)
      
    maxSlope = results['maxDampingSlope']['nominal']
    if type(np.nan) == type(maxSlope) or type(None) == type(maxSlope):
      disfCount = disfCount + 2
      slo = 0
      slof = 0 
    else:
      if len(maxSlope)==2:
        slo = 0 if maxSlope[0]>80 else 80-maxSlope[0]
        slof = abs(np.log10(maxSlope[1])-np.log10(1000))
      else:
        slo = 0
        slof = 0
        disfCount = disfCount + 1
	
    bandwidth = np.array(results['bw']['nominal'], dtype=float)
    results['THD_Hf'] = {}
    results['THD_Hf']['nominal'] = None
    if np.isnan(bandwidth):
      disfCount = disfCount + 1
      bw = 0
      thd_hf = 0
    else:
      bw = abs(np.log10(bandwidth) - np.log10(abs(1000)))
      resultsHf = runme2.evaluateActiveFilter_2_Hf(gen, indi, bandwidth)#-------------check THD at calculated (not fixed!) high freq 
      THD_Hf = np.array(resultsHf['THD_Hf']['nominal'], dtype=float)
      results['THD_Hf'] = {}
      results['THD_Hf']['nominal'] = resultsHf['THD_Hf']['nominal']#------------------put the Hf result in results
      if np.isnan(THD_Hf):
        disfCount = disfCount + 1
        thd_hf = 0
      else:
        thd_hf = THD_Hf-1 if THD_Hf > 1 else 0
            
    StaticOut = not results['isOutVNonStationary']['nominal']
    
    
    inImped = np.array(results['inimped']['nominal'], dtype=float)
    if np.isnan(inImped):
      disfCount = disfCount + 1
      inimped = 0
      outimped = 0
    else:
      outImped = np.array(results['outimped']['nominal'], dtype=float) #if inImped is measurable, then also outImped should be
      inimped = 0 if inImped > 1e5 else 1e5-inImped				# more is better (for voltage signals)
      outimped = outImped/1e2							# less is better (for voltage signals)
      
    
    if debug > 2:      
      print(slo, r, bw, d, slof,  100*StaticOut, thd_lf, thd_hf, g)
    
    #----------------------------------------------------------Score function SINGLE-OBJECTIVE
    if MOEAMODE == 0:
      score =  5*slo + 10*r + d + g*10 + 20*(2*bw + slof) + 100*StaticOut + 10*(thd_lf + thd_hf) + 1e-3*(inimped + outimped)  #+abs(slof-bw)/100)#rin! #TODO: slof ima tezave!!! poglej zgoraj, zakometirano 

      if disfCount > 0:
        score = 0 + np.exp(disfCount) * 1e3 + random.random()*10
      
    #----------------------------------------------------------Score function MULTI-OBJECTIVE
    else: #MOEAMODE == 1
      score = np.array([(d), r, (bw+slof)]) + (100*StaticOut + 10*(thd_lf + thd_hf) + 1*islp + g*10)  + 1e-3*(inimped + outimped) #+ slo#+abs(slof-bw)/100)#rin!

      if disfCount > 0:
        score = (np.array([0,0,0])+np.exp(disfCount) * 1e3) + np.array([random.random()*10,random.random()*10,random.random()*10])
      
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
