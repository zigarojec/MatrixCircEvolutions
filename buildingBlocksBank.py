"""
buildingBlocksBank.py
Ziga Rojec, EDA group, Faculty of Electrical Engineering, University of Ljubljana, Slovenia. Jan 2018.

This script is a starting point in which you tell the topology synthesis process WHICH and HOW MANY electronic devices shall be used in the search algorithm and therefore in the final solution. 
You can also set  MIN and MAX values for certain electronic devices parameters in paramBounds. 
This script is in strong relation to method makeNetlist from reproduction.py (! NOTE when changing the software !). 

Note that including many ('Quantity') elements will result in a huge (sum of every 'Quantity' times 'NofPins' in buildBlocks) connection matrix. This can result in longer "matrix to netlist" conversion time and can also widen the algorithm search space, which can - on one hand - increase possibilities of reaching an optimum solution - but on the other hand - significantly increase the time to find a solution. 

'Model' property: Multiple models can be listed. In normal mode only the first one is used. In nofailure mode all models are available to be used while simulating the robustness of the circuit. 

My advice - before you trigger the run, make sure you know, what kind of a circuit are you looking for. 
"""
import numpy as np
# This array contains the set of circuit nodes, that are accessible to the outer world. 
outerConns = ['vout','gnd']

"""
    {	#Simple resistor
'SpiceElementType': 'r',	      # How this element is encoded (initiated) in Spice netlist  
'Element': 'Rs',               # How do we name it
'Quantity': 2, #<---NOTE       # How many devices
'NofPins':  2,                 # How many terminals
'Model': '',                   # Name of the model(s)
'ParamTypes': {'r':'r'},       # Dictionary of parameter names vs. parameter types
    },    
"""


#NOTE: Set Quantity, do not touch others.
buildBlocks =  [
      {	#Simple resistor
	'SpiceElementType': 'r',	# How this element is encoded in Spice netlist
	'Element': 'Rs',
	'Quantity': 10, #<---NOTE
	'NofPins':  2,
	'Model': '',
	'ParamTypes': {'r':'r'},
	  },
      {	#Simple capacitor
	'SpiceElementType': 'c',	
	'Element':'Cs',
	'Quantity': 0,#<---NOTE
	'NofPins':  2,
	'Model': '',
	'ParamTypes': {'c':'c'},
	  },
      {	#Simple inductor
	'SpiceElementType': 'l',
	'Element':'Ls',
	'Quantity': 0,#<---NOTE
	'NofPins':  2,
	'Model': '',
	'ParamTypes': {'l':'l'},
	  },
      {	#Zener diode
	'SpiceElementType': 'd',
	'Element':'ZDs',
	'Quantity': 0,#<---NOTE
	'NofPins':  2,
	'Model': 'zd4v7',
	'ParamTypes': {},
	  },
      {	#NPN BJ Transistor
	'SpiceElementType': 'q',
	'Element':'NPNs',
	'Quantity': 0,#<---NOTE
	'NofPins':  3,
	'Model': 'bc238b',
	'ParamTypes': {},     
	  },
      
      {	#NPN BJ Transistor MULTI MODEL ROBUST CHECK
	'SpiceElementType': 'x',   # Since it is a complex subckt it is an x, not a q... 
	'Element':'NPNs',
	'Quantity': 6,#<---NOTE	# Mind the analysis in runme.py if hardcoded to some elements. 
	'NofPins':  3,
	'Model': ['T2N2222_resil_nom', 'T2N2222_resil_himp', 'T2N2222_resil_sck'], #,  
	'ParamTypes': {},     
	  },
      
      {	#PNP BJ Transistor MULTI MODEL ROBUST CHECK TEST
	'SpiceElementType': 'x',   # Since it is a complex subckt it is an x, not a q... 
	'Element':'PNPs',
	'Quantity': 6,#<---NOTE # Mind the analysis in runme.py if hardcoded to some elements.
	'NofPins':  3,
	'Model': ['2N2907_resil_nom', '2N2907_resil_himp', '2N2907_resil_sck'], #,  
	'ParamTypes': {},     
	  },      
      {	#Zener diode ROBUST MODE
	'SpiceElementType': 'x',	 # Since it is a combined subckt it is an x, not a q... 
	'Element':'ZDs',
	'Quantity': 0,#<---NOTE
	'NofPins':  2,
	'Model': ['zd4v7_resil_nom', 'zd4v7_resil_sck'], #, 'zd4v7_resil_himp'
	'ParamTypes': {},
	  },

      {	#Voltage supply
	'SpiceElementType': 'x', # Check subcircuit in models
	'Element':'VDC',
	'Quantity': 2,#<---NOTE
	'NofPins':  2,
	'Model': 'vdc',
	'ParamTypes': {'v':'v'},     
	  },



      {	#Subcircuit with three parallel PNPs
	'SpiceElementType': 'x',
	'Element':'3PNPs',
	'Quantity': 0,#<---NOTE
	'NofPins':  3,
	'Model': 'par3pnp',
	'ParamTypes': {},     
	  },
      {	#PNP BJ Transistor
	'SpiceElementType': 'q',
	'Element':'PNPs',
	'Quantity': 0,#<---NOTE
	'NofPins':  3,
	'Model': 'BC308B',
	'ParamTypes': {},
	  },
      {	#NMos transistor
	'SpiceElementType': 'x',
	'Element':'NMOSs',
	'Quantity': 0,#<---NOTE
	'NofPins':  3,
	'Model': 'submodn',
	'ParamTypes': {'w':'mos_w', 'l':'mos_l'},
	  },
      {	#PMos transistor
	'SpiceElementType': 'x',
	'Element':'PMOSs',
	'Quantity': 0,#<---NOTE
	'NofPins':  3,
	'Model': 'submodp',
	'ParamTypes': {'w':'mos_w', 'l':'mos_l'},
	  },
      
      {	#Opamp
	'SpiceElementType': 'x',
	'Element':'OPAMPSs',
	'Quantity': 0,#<---NOTE
	'NofPins':  3,
	'Model': 'LM348T_plus',
	'ParamTypes': {},
	  },
      {	#PMosCurrSrc1stg 				1.
	'SpiceElementType': 'x',
	'Element':'PMosCurrSrc1stg',
	'Quantity': 0,#<---NOTE
	'NofPins':  3,
	'Model': 'PMosCurrSrc1stg',
	'ParamTypes': {'w':'mos_w', 'l':'mos_l'},    
	  },
      {	#CascPMosCurrSrc1stg				2.
	'SpiceElementType': 'x',
	'Element':'CascPMosCurrSrc1stg',
	'Quantity': 0,#<---NOTE
	'NofPins':  3,
	'Model': 'CascPMosCurrSrc1stg',
	'ParamTypes': {'w':'mos_w', 'l':'mos_l'},    
	  },      
      {	#NMosAmp1ResOnSrc				3.
	'SpiceElementType': 'x',
	'Element':'NMosAmp1ResOnSrc',
	'Quantity': 0,#<---NOTE
	'NofPins':  3,
	'Model': 'NMosAmp1ResOnSrc',
	'ParamTypes': {'w':'mos_w', 'l':'mos_l', 'r':'r'},     
	  },   
      {	#BJTNPNCurrSink					4.
	'SpiceElementType': 'x',
	'Element':'BJTNPNCurrSink',
	'Quantity': 0,#<---NOTE
	'NofPins':  4,
	'Model': 'BJTNPNCurrSink',
	'ParamTypes': {},     
	  },
      {	#BJTPNPCurrSrc					5.
	'SpiceElementType': 'x',
	'Element':'BJTPNPCurrSrc',
	'Quantity': 0,#<---NOTE
	'NofPins':  4,
	'Model': 'BJTPNPCurrSrc',
	'ParamTypes': {},     
	  }, 
      {	#NMosCurrMirr					6.
	'SpiceElementType': 'x',
	'Element':'NMosCurrMirr',
	'Quantity': 0,#<---NOTE
	'NofPins':  4,
	'Model': 'NMosCurrMirr',
	'ParamTypes': {'w':'mos_w', 'l':'mos_l'},     
	  },
      {	#CascNMosCurrMirr				7.
	'SpiceElementType': 'x',
	'Element':'CascNMosCurrMirr',
	'Quantity': 0,#<---NOTE
	'NofPins':  4,
	'Model': 'CascNMosCurrMirr',
	'ParamTypes': {'w':'mos_w', 'l':'mos_l'},
	  },
      {	#PMosCurrSrc2Stg				8.
	'SpiceElementType': 'x',
	'Element':'PMosCurrSrc2Stg',
	'Quantity': 0,#<---NOTE
	'NofPins':  4,
	'Model': 'PMosCurrSrc2Stg',
	'ParamTypes': {
        'w1':'mos_w', 
        'l1':'mos_l', 
        'w2':'mos_w', 
        'l2': 'mos_l'},     
	  },
      {	#CascPMosCurrSrc2Stg				9.
	'SpiceElementType': 'x',
	'Element':'CascPMosCurrSrc2Stg',
	'Quantity': 0,#<---NOTE
	'NofPins':  4,
	'Model': 'CascPMosCurrSrc2Stg',
	'ParamTypes': {
        'w1':'mos_w', 
        'l1':'mos_l', 
        'w2':'mos_w', 
        'l2': 'mos_l'},       
	  },
      {	#PMosCascode					10.
	'SpiceElementType': 'x',
	'Element':'PMosCascode',
	'Quantity': 0,#<---NOTE
	'NofPins':  4,
	'Model': 'PMosCascode',
	'ParamTypes': {
        'w1':'mos_w', 
        'l1':'mos_l', 
        'w2':'mos_w', 
        'l2': 'mos_l'},     
	  },      
      {	#NMosCascode					11.
	'SpiceElementType': 'x',
	'Element':'NMosCascode',
	'Quantity': 0,#<---NOTE
	'NofPins':  4,
	'Model': 'NMosCascode',
	'ParamTypes': {
        'w1':'mos_w', 
        'l1':'mos_l', 
        'w2':'mos_w', 
        'l2': 'mos_l'},     
	  },  
]
    
    
    

    

paramBounds = {
  'r':{'min': 1e1,
       'max': 1e4,
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
	'v':{
		'min': 0, #V
       	'max': 12,
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


NofOutConns = len(outerConns)		#number of outerConnections - such as Vin, Vout, GND. NOTE: When changing this, one has also zo update the makeNetlist method from reproduction.py . 

# Moved to buildingBlocksBank.

#Global values NOTE: strictly follow the sequence in buildBlocks if changing the software!

#Calculating various global variables connected to matrix representation of a circuit
ALLPINS = np.array([], dtype=int)	#number of pins for each device type in an array
for element in buildBlocks:
  ALLPINS = np.append(ALLPINS,  element['NofPins']*element['Quantity'])

BigMatrixSize = 0
for element in buildBlocks:
  BigMatrixSize = BigMatrixSize + element['NofPins']*element['Quantity']
BigMatrixSize = BigMatrixSize + NofOutConns

NofPARAMS = 0
for element in buildBlocks:
  NofPARAMS = NofPARAMS + len(element['ParamTypes'])*element['Quantity']
  
for element in buildBlocks:
  if element['NofPins'] == 0:
      raise ValueError('An electrical element in buildingBlocksBank cannot have 0 (zero) pins (electrical terminals).')


	  
Nof2poles = 0
Nof3poles = 0
Nof4poles = 0
for element in buildBlocks:
  if element['NofPins'] == 2:
    Nof2poles += element['Quantity']
  if element['NofPins'] == 3:
    Nof3poles += element['Quantity']
  if element['NofPins'] == 4:
    Nof4poles += element['Quantity']



# Model netlist scheme
def returnModelScheme():
    """
    This function builds and returns a Python dictionary of models available for each device. The ModelScheme is ised to iterate over the various models for the same device in robust scoreCircuit evaluation procedure. 
    
    Place this function in e.g. utils.py. 
    Call it within somewhere else (e.g. main.py). 
    """
    
    modelScheme = dict()
    for buildingBlockTypeNo, buildingBlockType in enumerate(buildBlocks):
        for buildingBlockNo in range(1, buildingBlockType['Quantity']+1):   # Count, not index.
            #print(buildingBlockTypeNo, buildingBlockNo, buildingBlockType )
    
            if isinstance(buildingBlockType['Model'], list):
                #print(buildingBlockType['Model'])
                elementName = buildingBlockType['Element'] + '_' + str(buildingBlockNo) #Warning. Underscore is used as delimiter in score functions and netlister. 
                modelScheme[elementName] = buildingBlockType['Model']
    
    return modelScheme
    
MODELSCHEME = returnModelScheme()   # Make modelScheme a global var. 
