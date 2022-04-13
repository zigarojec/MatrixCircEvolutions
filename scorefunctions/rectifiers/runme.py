#!/usr/bin/env python

from pyopus.evaluator.performance import PerformanceEvaluator
from pyopus.evaluator.aggregate import formatParameters, Nbelow, Nabove, Rworst, Rexcluded, Aggregator
from pyopus.evaluator.auxfunc import paramList
import numpy as np

from pyopus.simulator.hspice import ipath

# ROBUST CIRCUIT EVOLUTION

def evaluate_rectifier(filename, **kwargs):
    """Evaluation script for pyopus. It evaluates performance of simple rectifier. 
    
    """    
    heads={
        'opus': {
            'simulator' : 'SpiceOpus',
            'moddefs': {
                #'def':     { 'file': '/spice/commonemitter/g_' + str(generation) + '_i_' + str(individual) + '_subckt.cir' },
                'def':     { 'file': str(filename) },
                'tb':      { 'file': 'topdc.cir'},
                'models':   {'file': 'models_for_start.inc'} # those files are !for now! to be accessible in main.py homefolder
                }, 
            'settings':{
                'debug': 0,
                'timeout' : 5, # Preizkusimo timeout 
                    },
            'params':{
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
            'params': {}, 
            'command': "op()",
        }, 
        'dc_sweep_tf': {
            'head': 'opus',
            'modules': [ 'def', 'tb', 'models' ], 
            'params': {}, 
            'saves': [ ],  
            'command': "dc(-10, 10, 'lin', 100, 'vin', 'dc')"
            },            
        }

    # This is a dictionary of Python variables that can be accessed in the 
    # expressions as var['name']
    #variables={
        #'saveInst': [ 'xmn2', 'xmn3', 'xmn1' ], 
        #'saveProp': [ 'vgs', 'vth', 'vds', 'vdsat' ]
    #}

            
    measures={
		'vout':{
			'analysis' : 'dc_sweep_tf',
			'corners' : [ 'nominal' ],
			'expression': 'v("vout")',
			'vector' : True
		 },        
		'scale' : {
			'analysis' : 'dc_sweep_tf',
			'corners' : [ 'nominal' ],
			'script': """__result = scale()""",
			'vector' : True,
		},
        'dcvout_rmse' : {
            'analysis' : 'dc_sweep_tf',
            'corners' : [ 'nominal' ],
            'script': """
#targets = np.abs(v("vin"))
targets = np.abs(scale())
outputs = v("vout")
#__result = np.sqrt((((outputs-outputs[0]) - (targets-targets[0])) ** 2).mean()) # Offset excluded!
__result = np.sqrt(((outputs - targets) ** 2).mean()) # Offset INCLUDED!
""",
            'vector' : False,
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
            'measure': 'dcvout_rmse',
            'norm': Nbelow(0.3, failure=10000.0), 	
            'reduce': Rworst()
        },         
    ]
    # End definition


    # Add definitions and analyses for nominal topology!
    if(len(kwargs)):
        if(kwargs["nominalTopology"] == True):
            # Append definition 
            definition.append(
                {
                    'measure': 'deviceActive',
                    'norm': Nabove(0.8, 1e-4, failure=10000.0), 	
                    'reduce': Rexcluded()#Rworst()
                }
            )
            # Append measure on nominal topology
            measures['deviceActive'] = {    # Count transisor that are "alive"  THIS SHOULD BE RUN ONLY ON NOMINAL TOPOLOGY NETLIST 
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
transistorList = ['xds_1','xds_2','xds_3', 'xds_4', 'xds_5', 'xds_6', 'xds_7', 'xds_8','xds_9','xds_10']  #
transistorActive = []
for t in transistorList:
    vpn = v(ipath('pint', [t, 'xcirc']), ipath('nint', [t, 'xcirc'])) #Vpn 
    #ipath() is universal pyopus function to access instance (node) e.g. cint in subcircuit (device) t that is part of another subcircuit xcirc. 
    ipn = i(ipath('vp', [t, 'xcirc']))    #Ipn
    x = scale()
    if(np.abs(vpn).mean() < 1e-12):     # Added against numerical noise of unconnected elements. Abs added!
        transistorActive.append(0.0)
    else:
        derV = m.dYdX(x, vpn)   # This works better than m.dYdX(vpn, x) !
        derI = m.dYdX(x, ipn)
        #transistorActive.append(int(abs(min(derV))>abs(1/100*max(derV)) and abs(min(derI))>abs(1/100*max(derI))))
        transistorActive.append(int(abs(min(derV))>abs(1e-4*max(derV))))
        
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
