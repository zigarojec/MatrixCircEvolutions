
from utils import *
from buildingBlocksBank import *
from reproduction import makeNetlist, makeNetlist_netlister
#import ... runme2
from  globalVars import *

import scorefunctions.arithmetic.runme as runme

import os, sys

def scoreCirc_squareroot_resilenceMode(circuit, gen, indi, MOEAMODE):
    """
    scoreCirc_squarerootcirc_resilenceMode is a wrapper function which calls a runme file. 
    In resilenceMode it adds multiple scores (costs) together. 
    ResilenceMode means running the analysis over failure-corners (special failure modes of each individual component).
    """

    debug = 0

    def evaluateNetlistAppendResults(circuitFilename, results_list):
        """
        Internal function for netlist evaluation call. 
        """
        #makeNetlist_netlister(circuit)
        cost, results = runme.evaluate_squarerootcirc(circuit.filename)
        results_list.append(results)
        return cost

    if debug > 2:
        print("\t\tG_" + str(gen) + "_I_" + str(indi))
    
    # Checking the connection matrix health
    OcSc, IcNc, SelfConnElm = checkConnsConnected(circuit.BigCircuitMatrix) #Outer connections Short cut, Inner connections Not connected   

    score = np.array([0,0,0], dtype="float64") if MOEAMODE == 1 else 0.0    
    score_array = np.array([], dtype="float64")
    
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
            
            # Nominal models evaluation.
            makeNetlist_netlister(circuit)  
            score, results = runme.evaluate_squarerootcirc(circuit.filename, nominalTopology=True)
            results_list.append(results)
            #input() # Wait and inspect every netlist please PROCEED AND TEST HERE
            os.remove(circuit.filename) #cleanup current subcircuit
            score_array = np.append(score_array, score)
            
            # Anti-models evaluation.
            for elm in MODELSCHEME:
                elm_num = elm.split("_")
                for mod in MODELSCHEME[elm][1:]:  # Consider only anti-models. The all nominal is already evaluated. 
                    #print("\t", mod, elm_num[0], elm_num[1])
                    # Duh... A little complicated now. 
                    makeNetlist_netlister(circuit, 
                                          Element = elm_num[0], 
                                          ElementNo =  int(elm_num[1]),
                                          ModelName = mod)
                    # ----------------------------------- #
                    score, results = runme.evaluate_squarerootcirc(circuit.filename)
                    results_list.append(results)
                    score_array = np.append(score_array, score)
                    # ----------------------------------- # 

                    if debug > 2:    
                        print("\t\tG_" + str(gen) + "_I_" + str(indi) + " SCORE:", score)
                        
                    #input("Here!") # Wait and inspect every netlist please PROCEED AND TEST HERE
                    os.remove(circuit.filename) #cleanup current subcircuit

            #----------------------------------------------------------Score function SINGLE-OBJECTIVE
            if MOEAMODE == 0:
                #score =  (1 + score_array.std())*score_array.sum()  # Mind the gap between device evaluation scores
                score = score_array.sum()

            #----------------------------------------------------------Score function MULTI-OBJECTIVE
            else: #MOEAMODE == 1
                score = np.array([score_array[0], score_array.max(), score_array.std()])

                #print "score64", score
                if gen > 30:
                    for i in range(0, len(score)): #round to 5 decimal points
                        score[i] = np.float16(score[i]) 
                    #print "score16", score
                score = np.ndarray.tolist(score) #for faster non-dominant sorting let the score be python array instead of np.array 5.7.2017
            #-------------------------------------------------------------------

            
            return score, results_list
        
        else:
            # Nominal models evaluation.
            makeNetlist_netlister(circuit)  
            score, results = runme.evaluate_squarerootcirc(circuit.filename, nominalTopology=True)
            os.remove(circuit.filename) #cleanup current subcircuit
    return score, results
