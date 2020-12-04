"""
scoreFunctions.py
Ziga Rojec, EDA group, Faculty of Electrical Engineering, University of Ljubljana, Slovenia. Jan 2018.

In scoreFunctions.py you will find COST functions that distinguish good circuits from bad circuits. 

In May 2018, four cist functions are provided:
scoreCirc_ActiveFilter_2
scoreCirc_PassiveBandPass
scoreCirc_HighPass
scoreCirc_CmosVoltageReference_2

Altering those you can define, what kind of a circuit do you want at the end. Cost functions are in strong relation to script runme2.py. 
"""


from utils import *
from reproduction import makeNetlist
import runme2
import os, sys
import globalVars
#import signal

debug = 0




#Hybrid LGA+MOEA score functions

def scoreCirc_ActiveFilter_2(circuit, gen, indi, MOEAMODE):#LGA+MOEA
  """
  Anchestor: scoreCirc_ActiveFilter_MOEA
  This function takes:
    - 	a circuit object, 
    - 	generetion number,
    - 	individual number,
    -	0 - if score function is to be used in single-objective mode (score is a single float)
	1 - if score function is to be used in multi-objective mode (score is a vector of three objectives)
	
  This function returns:
    - 	score (float value (or array of floats) of the objective(s))
    -	results (dictionary from runme2 script)
  
  This function measures THD at low freq default, and at bw that is measured via evaluation function.  
  """  
  if debug > 2:
    print "\t\tG_" + str(gen) + "_I_" + str(indi)

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
      print slo, r, bw, d, slof,  100*StaticOut, thd_lf, thd_hf, g
    
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
      print "\t\tG_" + str(gen) + "_I_" + str(indi) + " SCORE:", score
    
    filename = "g_" + str(gen) + "_i_" + str(indi) + "_subckt.cir"
    os.remove(filename) #cleanup current subcircuit
    
  return score, results  #removed matrixQuaziID and matrixDensity!!


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
    print "\t\tG_" + str(gen) + "_I_" + str(indi)

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
      print r, dL, dH, g, bwL, bwH
    
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
      print "\t\tG_" + str(gen) + "_I_" + str(indi) + " SCORE:", score
    
    filename = "g_" + str(gen) + "_i_" + str(indi) + "_subckt.cir"
    os.remove(filename) #cleanup current subcircuit
    
  return score, results  #removed matrixQuaziID and matrixDensity!!

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
    print "\t\tG_" + str(gen) + "_I_" + str(indi)

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
      print r, dL, g, bwL
    
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
      print "\t\tG_" + str(gen) + "_I_" + str(indi) + " SCORE:", score
    
    filename = "g_" + str(gen) + "_i_" + str(indi) + "_subckt.cir"
    os.remove(filename) #cleanup current subcircuit
    
  return score, results  #removed matrixQuaziID and matrixDensity!!


def scoreCirc_CmosVoltageReference_2(circuit, gen, indi, MOEAMODE):
  """Calculates the score of the voltage reference circuit given in argument. """
  
  if debug > 2:
    print "\t\tG_" + str(gen) + "_I_" + str(indi)
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
      print "\t\tG_" + str(gen) + "_I_" + str(indi) + " SCORE:", score

    filename = "g_" + str(gen) + "_i_" + str(indi) + "_subckt.cir"
    #os.remove(filename) #cleanup current subcircuit
    
    
    # TRIGGER STOP SIGNAL if:
    if (vdd_s_t2 <= 0.001 and 
	psrr >= 80 and 
	powe <= 1e-5):
      globalVars.DONE = 1 # End evolution, feasible solution evolved.
    

  return score, results  






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

      
    print "\t\t\t\tG_" + str(gen) + "_I_" + str(indi) + " SCORE:", score
    #cleanup current subcircuit
    filename = "g_" + str(gen) + "_i_" + str(indi) + "_subckt.cir"
    os.remove(filename)
    #print ".",
    #circuit.objectivesScore = copy(score)	#id does not work with mpirun since mpirun works with copies
    #circuit.matrixDensity = matrixDensity
  return score, matrixDensity, matrixQuaziID, results

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


