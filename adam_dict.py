
outer_conns = 'outer_conns'
inner_netlist = 'inner_netlist'

primates = {
  'adam_circuit_name':{
    outer_conns: ['gnd', 'vin', 'vout'],
    inner_netlist: {
	'r1':[
	  ['vin',1], [10e3],
	],
	'r12':[
	  [1,'vout'], [10e3],
	],
	'r2':[
	  ['vout','gnd'], [10e3],
	],
    }
  },

  'filter2':{
    outer_conns: [],    
    inner_netlist: {
	'r1':[
	  ['vin',1], [10e3],
	],
	'dz1':[
	  [1,'vout'], 'zd4v7',
	],
	'r2':[
	  ['vout','gnd'], [10e3],
	],      
      
      }
  }
}


"""
Je bilo ukvarjanje z matrikami zabloda? 
Če vezje zapišemo v obliki Python slovarja, bi lahko mutacije izvajal skoraj neposredno na slovarju. To bi prineslo znatno pospešitev.

Edini pomislek, ki ga imamo je, da se z vidika gledanja na postopek kot genetski evolucijski postopek "lastnosti" osebkov ne bi dedovale neposredno s pomočjo križanja. 

Dokazati bi bilo potrebno, da križanja matrik nimajo znatnega vpliva na potek evolucije v primerjavi z mutacijami. 
"""

  
def checkNetlistDictionary(primates):
  """
  This is a test procedure. 
  Checks whether my netlist description dictionary complies with building blocks, stated in buildingBlocksBank.py
  Throws error, if primate_circuit netlist dictionary does not comply. 
  """
  string = ""
  string += "Checking circuit description against building blocks bank:" + "\n"
  for primate in primates:
    string += primate + "..." + "\n"
    try:
      if not (len(primates[primate][outer_conns]) >= 2):
	string += "\t" + primate + "not enough (<2) outer terminals" + "\n"
    except:
      string += "\t" + outer_conns + "field missing." + "\n"
      
    try:
      for elm in primates[primate][inner_netlist]:		# Go through every listed element of circuit.
	if elm[0] != 'r' and elm[0] != 'c' and elm[0] != 'l':	# If not rlc, then model name is wanted.
	  if isinstance(primates[primate][inner_netlist][elm][1], str):	# Is model name (string) listed or not?
	    string += "ok" + "\n"
	  else:
	    string += primate + "misses model name (in string) for element" + element + "\n"
	len(primates[primate][inner_netlist][elm])
      
    except:
      string += "\t" + inner_netlist + "field probably missing." + "\n"
      
  print(string)

checkNetlistDictionary(primates)






netlist = ""
adam = primates['adam_circuit_name']
NofOutConns = len(adam[outer_conns])


with open("adame.cir", "wt") as file:
  netlist += "\n"
  netlist += ".subckt "
  file.write(netlist)
  file.close()

primates['adam_circuit_name']


