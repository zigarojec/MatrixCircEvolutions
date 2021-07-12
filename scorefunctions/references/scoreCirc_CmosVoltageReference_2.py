def scoreCirc_CmosVoltageReference_2(circuit, gen, indi, MOEAMODE):
  """Calculates the score of the voltage reference circuit given in argument. """
  
  if debug > 2:
    print("\t\tG_" + str(gen) + "_I_" + str(indi))
  #----------#
  VREF = 1.5
  #----------#

  #---------------------------------------------------------BigMatrix stuff, check short-circuits, matrix density, matrix identifier (obsolete)  
  FullBigCircuitMatrix = copy(circuit.fullRedundancyMatrix)
  OcSc, IcNc, SelfConnElm = checkConnsConnected(FullBigCircuitMatrix) #Outer connections Short cut, Inner connections Not connected
  #---------------------------------------------------------  
  
  score = np.array([0,0,0], dtype="float64") if MOEAMODE == 1 else 0
  
  score += 2e4*np.exp(OcSc)
  results = None
  if OcSc > 1:
    score += 1e4*np.exp(OcSc)
  else:
    #----------------------------------------------------------Try to make netlist and evaluate the individual
    makeNetlist(circuit, gen, indi, FullBigCircuitMatrix)
    results = runme2.evaluateCmosVoltageRef(gen, indi)
    #----------------------------------------------------------Start of results analysis and objectives creation
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
      
    psrr = results['psrr']['nominal']
#    if np.isnan(np.array(psrr, dtype=float)):
#      disfCount = disfCount + 1
#      psr = 0
#    else:
#      psr = 1.0/psrr #abs(90 - psrr) if psrr < 90 else 0 #tole kot objective ni ok. ker je opravljena meritev samo pri vdd=15 je to precej stala.


    #----------------------------------------------------------Score function SINGLE-OBJECTIVE
    if MOEAMODE == 0:
      score =(vdd_s_t1 + 5*vdd_s_t1_d +
	      2*vdd_s_t2 + 2*vdd_s_t2_d +
	      vdd_s_t3 + 5*vdd_s_t3_d +
	      #vdd_s_r1 + 2*vdd_s_r1_d +
	      #vdd_s_r2 + 2*vdd_s_r2_d + 
	      #vdd_s_r3 + 2*vdd_s_r3_d + 
	      (100*powe)
      )
      if disfCount > 0:
        score = 0 + np.exp(disfCount) * 1e3
	
    #----------------------------------------------------------Score function MULTI-OBJECTIVE	
    else: #MOEAMODE == 1:
      oMediana = vdd_s_t1 + vdd_s_t2 + vdd_s_t3
      oPsrr = vdd_s_t1_d + vdd_s_t2_d + vdd_s_t3_d	#DC rejection
      #oPsrr = psr
      oP = powe
					      #add constraints
      score = (np.array([oMediana, oPsrr, oP]) 	+ (oMediana if oMediana > 4 else 0) + 
						#+ (oPsrr*1000 if oPsrr > 1.0/40 else 0) +
						+ (oPsrr if oPsrr > 3 else 0) +
						+ (oP if oP > 1e-1 else 0)
      )
      if disfCount > 0:
        score = (np.array([0,0,0])+np.exp(disfCount) * 1e3) + random.randint(0, 200)

    #-------------------------------------------------------------------
    if debug > 2:    
      print("\t\tG_" + str(gen) + "_I_" + str(indi) + " SCORE:", score)

    filename = "g_" + str(gen) + "_i_" + str(indi) + "_subckt.cir"
    #os.remove(filename) #cleanup current subcircuit
    
    
    # TRIGGER STOP SIGNAL if:
    if (vdd_s_t2 <= 0.001 and 
	psrr >= 80 and 
	powe <= 1e-5):
      globalVars.DONE = 1 # End evolution, feasible solution evolved.
    

  return score, results  
