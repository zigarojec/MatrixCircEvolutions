"""
buildingBlocksBank.py
Ziga Rojec, EDA group, Faculty of Electrical Engineering, University of Ljubljana, Slovenia. Jan 2018.

This script is a starting point in which you tell the topology synthesis process WHICH and HOW MANY electronic devices shall be used in the search algorithm and therefore in the final solution. 
You can also set  MIN and MAX values for certain electronic devices parameters in paramBounds. 
This script is in strong relation to method makeNetlist from reproduction.py (! NOTE when changing the software !). 

Note that including many ('Quantity') elements will result in a huge (sum of every 'Quantity' times 'NofPins' in buildBlocks) connection matrix. This can result in longer "matrix to netlist" conversion time and can also widen the algorithm search space, which can - on one hand - increase possibilities of reaching an optimum solution - but on the other hand - significantly increase the time to find a solution. 

My advice - before you trigger the run, make sure you know, what kind of a circuit are you looking for. 
"""
#NOTE: Set Quantity, do not touch others.
buildBlocks =  [
      {	#Simple resistor
	'SpiceElementType': 'r',	# How this element is encoded in Spice netlist
	'Element': 'Rs',
	'Quantity': 2, #<---NOTE
	'NofPins':  2,
	'Model': '',
	'ParamTypes': ['r'],
	  },
      {	#Simple capacitor
	'SpiceElementType': 'c',	
	'Element':'Cs',
	'Quantity': 0,#<---NOTE
	'NofPins':  2,
	'Model': '',
	'ParamTypes': ['c'],
	  },
      {	#Simple inductor
	'SpiceElementType': 'l',
	'Element':'Ls',
	'Quantity': 0,#<---NOTE
	'NofPins':  2,
	'Model': '',
	'ParamTypes': ['l'],
	  },
      {	#Zener diode
	'SpiceElementType': 'd',
	'Element':'ZDs',
	'Quantity': 0,#<---NOTE
	'NofPins':  2,
	'Model': 'zd4v7',
	'ParamTypes': [],
	  },
      {	#NPN BJ Transistor
	'SpiceElementType': 'q',
	'Element':'NPNs',
	'Quantity': 0,#<---NOTE
	'NofPins':  3,
	'Model': 'bc238b',
	'ParamTypes': [],     
	  },
      {	#Subcircuit with three parallel PNPs
	'SpiceElementType': 'x',
	'Element':'3PNPs',
	'Quantity': 2,#<---NOTE
	'NofPins':  3,
	'Model': 'par3pnp',
	'ParamTypes': [],     
	  },
      {	#PNP BJ Transistor
	'SpiceElementType': 'q',
	'Element':'PNPs',
	'Quantity': 1,#<---NOTE
	'NofPins':  3,
	'Model': 'BC308B',
	'ParamTypes': [],
	  },
      {	#NMos transistor
	'SpiceElementType': 'x',
	'Element':'NMOSs',
	'Quantity': 4,#<---NOTE
	'NofPins':  3,
	'Model': 'submodn',
	'ParamTypes': ['mos_w', 'mos_l'],
	  },
      {	#PMos transistor
	'SpiceElementType': 'x',
	'Element':'PMOSs',
	'Quantity': 6,#<---NOTE
	'NofPins':  3,
	'Model': 'submodp',
	'ParamTypes': ['mos_w', 'mos_l'],
	  },
      
      {	#Opamp
	'SpiceElementType': 'x',
	'Element':'OPAMPSs',
	'Quantity': 0,#<---NOTE
	'NofPins':  3,
	'Model': 'LM348T_plus',
	'ParamTypes': [],
	  },
      {	#PMosCurrSrc1stg 				1.
	'SpiceElementType': 'x',
	'Element':'PMosCurrSrc1stg',
	'Quantity': 0,#<---NOTE
	'NofPins':  3,
	'Model': 'PMosCurrSrc1stg',
	'ParamTypes': ['mos_w', 'mos_l'],     
	  },
      {	#CascPMosCurrSrc1stg				2.
	'SpiceElementType': 'x',
	'Element':'CascPMosCurrSrc1stg',
	'Quantity': 0,#<---NOTE
	'NofPins':  3,
	'Model': 'CascPMosCurrSrc1stg',
	'ParamTypes': ['mos_w', 'mos_l'],     
	  },      
      {	#NMosAmp1ResOnSrc				3.
	'SpiceElementType': 'x',
	'Element':'NMosAmp1ResOnSrc',
	'Quantity': 0,#<---NOTE
	'NofPins':  3,
	'Model': 'NMosAmp1ResOnSrc',
	'ParamTypes': ['mos_w', 'mos_l', 'r'],     
	  },   
      {	#BJTNPNCurrSink					4.
	'SpiceElementType': 'x',
	'Element':'BJTNPNCurrSink',
	'Quantity': 0,#<---NOTE
	'NofPins':  4,
	'Model': 'BJTNPNCurrSink',
	'ParamTypes': [],     
	  },
      {	#BJTPNPCurrSrc					5.
	'SpiceElementType': 'x',
	'Element':'BJTPNPCurrSrc',
	'Quantity': 0,#<---NOTE
	'NofPins':  4,
	'Model': 'BJTPNPCurrSrc',
	'ParamTypes': [],     
	  }, 
      {	#NMosCurrMirr					6.
	'SpiceElementType': 'x',
	'Element':'NMosCurrMirr',
	'Quantity': 0,#<---NOTE
	'NofPins':  4,
	'Model': 'NMosCurrMirr',
	'ParamTypes': ['mos_w', 'mos_l'],     
	  },
      {	#CascNMosCurrMirr				7.
	'SpiceElementType': 'x',
	'Element':'CascNMosCurrMirr',
	'Quantity': 0,#<---NOTE
	'NofPins':  4,
	'Model': 'CascNMosCurrMirr',
	'ParamTypes': ['mos_w', 'mos_l'],
	  },
      {	#PMosCurrSrc2Stg				8.
	'SpiceElementType': 'x',
	'Element':'PMosCurrSrc2Stg',
	'Quantity': 0,#<---NOTE
	'NofPins':  4,
	'Model': 'PMosCurrSrc2Stg',
	'ParamTypes': ['mos_w', 'mos_l', 'mos_w', 'mos_l'],     
	  },
      {	#CascPMosCurrSrc2Stg				9.
	'SpiceElementType': 'x',
	'Element':'CascPMosCurrSrc2Stg',
	'Quantity': 0,#<---NOTE
	'NofPins':  4,
	'Model': 'CascPMosCurrSrc2Stg',
	'ParamTypes': ['mos_w', 'mos_l', 'mos_w', 'mos_l'],     
	  },
      {	#PMosCascode					10.
	'SpiceElementType': 'x',
	'Element':'PMosCascode',
	'Quantity': 0,#<---NOTE
	'NofPins':  4,
	'Model': 'PMosCascode',
	'ParamTypes': ['mos_w', 'mos_l', 'mos_w', 'mos_l'],     
	  },      
      {	#NMosCascode					11.
	'SpiceElementType': 'x',
	'Element':'NMosCascode',
	'Quantity': 0,#<---NOTE
	'NofPins':  4,
	'Model': 'NMosCascode',
	'ParamTypes': ['mos_w', 'mos_l', 'mos_w', 'mos_l'],     
	  },  
]

paramBounds = {
  'r':{'min': 1e3,#Spremenjeno 6.9.2017
       'max': 1e7,
       'scale': 96,	#---not used--------scales (E):	96, 	48, 	24, 	12,	None for linear
	 },		#---------------------------	1%	2%	5%	10%
  'c':{'min': 1e-12,
       'max': 1e-6,
       'scale': 12,
	 },
  'l':{'min': 1e-9,
       'max': 1e-3,
       'scale': None,
	 },
  'mos_w':{'min': 1.8e-07,
       'max': 0.0001,
       'scale': None,
	 },
  'mos_l':{'min': 1.8e-07,
       'max': 4e-06,
       'scale': None,
	 },

  }
    
  #TODO: multiplication parameter (mos)
  #TODO: parameter values based on standard scale
	# add "standard scale" as option
	# problem for PSADE - optimizes continiousely
