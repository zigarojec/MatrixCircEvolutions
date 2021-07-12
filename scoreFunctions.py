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
from buildingBlocksBank import *

from reproduction import makeNetlist, makeNetlist_netlister
import runme2
import os, sys
import globalVars
#import signal


debug = 0

# TODO Make every scoreFunction a single file. This is getting to damn big!!

 
def scoreCirc_commonEmitterAmp_resilenceMode(circuit, gen, indi, MOEAMODE):
    """
    # Slice the current scoreFunction into checking the BigCircuitMatrix here
    # When creating the netlist, go to scoreCirc
    
    # There is no need to complicate actually.
    # Just create an interator in the else statement in scorecirc that would work with or without existing modelScheme...
    """
    
    def evaluateNetlistAppendResults(circuitFilename, results_list):
        """
        Internal function for netlist evaluation call. 
        """
        #makeNetlist_netlister(circuit)
        results = runme2.evaluate_CommonEmitterAmp(circuit.filename)
        results_list.append(results)
        #print(results)
        disfCount = 0
        dcgain = np.array(results['DCgain']['nominal'], dtype=float)
        if np.isnan(dcgain):
            disfCount += 1
            _dcgain = 0      
        else:
            _dcgain = abs(dcgain - 73.6) if dcgain < 73.6 else 0          
        
        dcvout_rmse = np.array(results['dcvout_rmse']['nominal'], dtype=float)
        if np.isnan(dcvout_rmse):
            disfCount += 1
            _dcvout_rmse = 0      
        else:
            _dcvout_rmse = abs(dcvout_rmse-7.6e-3)        # Nominal value deducted. 
        
        maxpower = np.array(results['maxpower']['nominal'], dtype=float)
        if np.isnan(maxpower):
            disfCount += 1
            _maxpower = 0      
        else:
            _maxpower = abs(maxpower-(7.1e-3))          # Nominal value deducted. 
    
        gain_stddev_norm = np.array(results['gain_stddev_norm']['nominal'], dtype=float)
        if np.isnan(gain_stddev_norm):
            disfCount += 1
            _gain_stddev_norm = 0      
        else:   # If stddev lower than 1% it does not count. 
            _gain_stddev_norm = abs(gain_stddev_norm) if gain_stddev_norm > 0.01 else 0   
                    
        
        return disfCount, _dcgain, _dcvout_rmse, _maxpower, _gain_stddev_norm
    
    
    if debug > 2:
        print("\t\tG_" + str(gen) + "_I_" + str(indi))
    
    # Checking the connection matrix health
    OcSc, IcNc, SelfConnElm = checkConnsConnected(circuit.BigCircuitMatrix) #Outer connections Short cut, Inner connections Not connected   

    score = np.array([0,0,0], dtype="float64") if MOEAMODE == 1 else 0
    score_array = np.array([0,0,0,0,0], dtype="float64") # Watch now: this is the real numbers of all objectives. Optimize that code.
    # [disfCount, _dcgain, _dcvout_rmse, _maxpower, _gain_stddev_norm]
    
    results = None  
    if OcSc > 1:    # if outer connections are short-circuited, punish the individual based on the matrix health and exit.
        score += 1e4*np.exp(OcSc)
        # for the sake of problems in sorting when scores are equal, add a small random difference
        score += np.array([random.random()*10,random.random()*10,random.random()*10]) if MOEAMODE == 1 else random.random()*10    
    else:
        # Create netlist and calculate the cost function. 
        # Based on all models that are being used, create a series of evaluations, using the models listed in MODELSCHEME!
        if robustMode:
            results_list = []
            
            makeNetlist_netlister(circuit)  # Nominal models evaluation.
            disfCount, _dcgain, _dcvout_rmse, _maxpower, _gain_stddev_norm = evaluateNetlistAppendResults(circuit.filename, results_list)
            os.remove(circuit.filename) #cleanup current subcircuit
            score_array += np.array([disfCount, _dcgain, _dcvout_rmse, _maxpower, _gain_stddev_norm])
            
            for elm in MODELSCHEME:
                elm_num = elm.split("_")
                for mod in MODELSCHEME[elm][1:]:  # Consider only anti-models. The all nominal is already evaluated. 
                    #print("\t", mod, elm_num[0], elm_num[1])
                    # TODO Solving the repeated nominal model analysis!
                    # Duh... A little complicated now. 
                    makeNetlist_netlister(circuit, 
                                          Element = elm_num[0], 
                                          ElementNo =  int(elm_num[1]),
                                          ModelName = mod)
                    # ----------------------------------- #
                    disfCount, _dcgain, _dcvout_rmse, _maxpower, _gain_stddev_norm = evaluateNetlistAppendResults(circuit.filename, results_list)
                    score_array += np.array([disfCount, _dcgain, _dcvout_rmse, _maxpower, _gain_stddev_norm])
                    # ----------------------------------- # 

                    if debug > 2:    
                        print("\t\tG_" + str(gen) + "_I_" + str(indi) + " SCORE:", score)
                        
                    #input() # Wait and inspect every netlist please PROCEED AND TEST HERE
                    os.remove(circuit.filename) #cleanup current subcircuit

            disfCount = score_array[0] # TOTAL disfCount
            #----------------------------------------------------------Score function SINGLE-OBJECTIVE
            if MOEAMODE == 0:
                score +=  sum(score_array[1:]) # disfCount are left out

                if disfCount > 0:
                    score = 0 + np.exp(disfCount) * 1e3 + random.random()*10
                
                if debug > 2:
                    print(_dcgain, _dcgain, _dcvout_rmse, _maxpower, _gain_stddev_norm)
            
            #----------------------------------------------------------Score function MULTI-OBJECTIVE
            else: #MOEAMODE == 1
                score = np.array([_dcgain, _maxpower, _gain_stddev_norm]) + 0 # + _dcvout_rmse
                if disfCount > 0:
                    score = (np.array([0,0,0])+np.exp(disfCount) * 1e3) + np.array([random.random()*10,random.random()*10,random.random()*10])
                
                #print "score64", score
                if gen > 30:
                    for i in range(0, len(score)): #round to 5 decimal points
                        score[i] = np.float16(score[i]) 
                    #print "score16", score
                score = np.ndarray.tolist(score) #for faster non-dominant sorting let the score be python array instead of np.array 5.7.2017
            #-------------------------------------------------------------------

            
            return score, results_list
        
        else:
            print("TODO!")
            # What TODO when no robustMode ... 

    return score, results















