#!/usr/bin/env python

from pyopus.evaluator.performance import PerformanceEvaluator
from pyopus.evaluator.aggregate import formatParameters, Nbelow, Nabove, Rworst, Rexcluded, Aggregator
from pyopus.evaluator.auxfunc import paramList
import numpy as np

from pyopus.simulator.hspice import ipath

# ROBUST CIRCUIT EVOLUTION

def evaluate_CommonEmitterAmp(filename, **kwargs):
    """Evaluation script for pyopus. It evaluates performance of simple amplifier. 
    
    """    
    heads={
        'opus': {
            'simulator' : 'SpiceOpus',
            'moddefs': {
                #'def':     { 'file': '/spice/commonemitter/g_' + str(generation) + '_i_' + str(individual) + '_subckt.cir' },
                'def':     { 'file': str(filename) },
                'tb':      { 'file': 'topdc_robust_commonemitter.cir'},#'testTopCircuit_hspice.cir' }, 
                'models':   {'file': 'models_for_start.inc'} # those files are !for now! to be accessible in main.py homefolder
                }, 
            'settings':{
                'debug': 0,
                'timeout' : 5, # Preizkusimo timeout 
                    },
            'params':{
                'temperature': 25,
                'iref': 50e-6
                    },
            'options': {
                'method': 'trap', 
                'noautoconv': True,		#prevent Spice to try to converge for every price
                    }   
            }
        }
    analyses={
        'op': {
            'head': 'opus', 
            'modules': [ 'def', 'tb', 'models' ], 
            'params': {'rl':1000e3}, 
            'command': "op()",
        }, 
        'dc_sweep_tf': {
            'head': 'opus',
            'modules': [ 'def', 'tb', 'models' ], 
            'params': {'rl':1000e3}, 
            'saves': [ ],  
            'command': "dc(-100e-6, 100e-6, 'lin', 100, 'iin', 'dc')" #if stop point changed, change also in reproduction cost function. 
            },            
        }

    # This is a dictionary of Python variables that can be accessed in the 
    # expressions as var['name']
    variables={
        #'saveInst': [ 'xmn2', 'xmn3', 'xmn1' ], 
        #'saveProp': [ 'vgs', 'vth', 'vds', 'vdsat' ]
    }

            
    measures={
        'DCgain' : {
            'analysis' : 'dc_sweep_tf',
            'corners' : [ 'nominal' ],
            'expression': """
gain = abs(m.DCgain(v("vout"), scale()))
__result = m.gain2dB(gain, unit="db20")
                """,
            'vector' : False,
        },         
        'dcvout_rmse' : {
            'analysis' : 'dc_sweep_tf',
            'corners' : [ 'nominal' ],
            'script': """
targets = 20e3*scale() + 0.66 #linear gain of 4.8kV/A with offset 0.66
outputs = v("vout")
__result = np.sqrt((((outputs-outputs[0]) - (targets-targets[0])) ** 2).mean()) # Offset excluded!
""",
            'vector' : False,
        },           
        'maxpower':{
            'analysis' : 'dc_sweep_tf',
            'corners' : [ 'nominal' ],
            'expression': '(-v("vdd")*i("vdd")).max()',
            'vector' : False
        },
        # Normalised gain standard deviation (checking if gain is a constant...!)
        'gain_linearity':{    
            'analysis' : 'dc_sweep_tf',
            'corners' : [ 'nominal' ],
            'expression': """
deriv = m.gain2dB(abs(m.dYdX(v("vout"), scale())))
derivdiff = deriv.max()-deriv.min()
#__result =  abs((np.std(deriv, dtype="float64")/maxderiv)*100) # norm the std deviation with maximum gain
__result = derivdiff
""",
            'vector' : True 
            },      
        
        }
    corners = {
        'nominal': {
            'heads': [ 'opus' ], 
            'modules': [ ],            # modules must be added in corners now!
            'params': {
                'temperature': 25, 
                'vcc': 15, 
                }
            }
        }		

    # Definition
    definition = [
        {
            'measure': 'DCgain',
            'norm': Nabove(70, failure=10000.0), 	
            'reduce': Rworst()
        },
        {
            'measure': 'dcvout_rmse',
            'norm': Nbelow(1e-3, failure=10000.0), 	
            'reduce': Rworst()
        },        
        {
            'measure': 'maxpower',
            'norm': Nbelow(10e-3, failure=10000.0), 	
            'reduce': Rexcluded()
        },     
        {
            'measure': 'gain_linearity',
            'norm': Nbelow(0.1, failure=10000.0), 	
            'reduce': Rworst()
        },   
    ]
    # End definition


    # Add definitions an analyses for nominal topology!
    if(len(kwargs)):
        if(kwargs["nominalTopology"] == True):
            # Append definition 
            definition.append(
                {
                    'measure': 'transistorActive',
                    'norm': Nabove(0.8, 1e-3, failure=10000.0), 	
                    'reduce': Rworst()
                }
            )
            # Append measure on nominal topology
            measures['transistorActive'] = {    # Count transisor that are "alive"  THIS SHOULD BE RUN ONLY ON NOMINAL TOPOLOGY NETLIST 
            'analysis' : 'dc_sweep_tf',
            'corners' : [ 'nominal' ],
            #'expression': "v('q1:xnpn_1:xcirc#collector')", # WHAT TODO TO to access subnodes?? 
            # Ja. Ne dela zato, ker se zunanji nodi podvezja preimenujejo! Notranji nodi pa se ne.  
            
            # In this expression the followng is defined. 
            # We check Vce voltage across a transistor. Also Ice current.
            # Then we calculate a derivative over the input signal. 
            # If both dice/dscale and dvce/dscale produce at least 1/100 of maximal derivative, 
            # then a transistor is somehow connected to the input signal, or "is alive". 
            # TODO 1: Return a vector of active transistors. Calculate the procent in another analysis. 
            # TODO 2: transistorList has to be imported from buildingBlocksBank somehow...


            'script': """   
transistorList = ['xnpns_1', 'xnpns_2', 'xnpns_3', 'xnpns_4', 'xpnps_1', 'xpnps_2', 'xpnps_3', 'xpnps_4']  #
transistorActive = []
for t in transistorList:
    #vce = v('cint:' + t + ':xcirc', 'eint:' + t + ':xcirc') #Vce 
    vce = v(ipath('cint', [t, 'xcirc']), ipath('eint', [t, 'xcirc'])) #Vce 
    #ipath() is universal pyopus function to access instance (node) e.g. cint in subcircuit (device) t that is part of another subcircuit xcirc. 
    #ice = i('vc:' + t + ':xcirc')    #Ice
    ice = i(ipath('vc', [t, 'xcirc']))    #Ice
    x = scale()

    if(vce.mean() < 1e-12):     # Added against numerical noise of unconnected elements
        transistorActive.append(0.0)
    else:
        derV = m.dYdX(x, vce)   # This works better than m.dYdX(vce, x) !
        derI = m.dYdX(x, ice)
        transistorActive.append(int(abs(min(derV))>abs(1/100*max(derV)) and abs(min(derI))>abs(1/100*max(derI))))
        
__result =  float(sum(transistorActive))/len(transistorList)
""",
            'vector' : False
            }  


    # Parameters
    params = { 

    }
    # End parameters

    # Parameter order
    inOrder=list(params.keys())
    inOrder.sort()
    # End parameter order

    pe=PerformanceEvaluator(heads, analyses, measures, corners, debug=0)
    
    results,ancount=pe()   
    
    # Aggregate cost function
    ce=Aggregator(pe, definition, inOrder, debug=0)   


    # Vectorize parameters
    x=paramList(params, inOrder)

    # Evaluate aggregate function at vector x 
    cf=ce(x)

    # Print the results
    #print("")
    #print("cost=%e" % cf)
    #print(ce.formatParameters())
    #print(ce.formatResults(nMeasureName=15, nCornerName=15))
    #print("")

    # Print analysis count 
    # print("Analysis count: "+str(pe.analysisCount))

    # Cleanup intemediate files
    pe.finalize()
    return cf, results
