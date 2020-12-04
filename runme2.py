"""
runme2.py
Ziga Rojec, EDA group, Faculty of Electrical Engineering, University of Ljubljana, Slovenia. Jan 2018.

You will find simulation and measurement setups for PyOpus here. With those functions you actually execute the SPICE simulations on a certain netlist. Those funstions are used within cost functions (scoreFunctions.py).

"""


from pyopus.evaluator.performance import PerformanceEvaluator
from pyopus.evaluator.aggregate import formatParameters
import numpy as np
#import pyopus.wxmplplot.plotitf as pyopl


#generation = 1
#individual = 1

def evaluateCmosVoltageRef(generation, individual):
	"""Evaluation script for pyopus. It evaluates performance of an voltage reference. """
	heads={
		'opus': {
			'simulator' : 'HSpice',
			'moddefs': {
				'def':     { 'file': 'g_' + str(generation) + '_i_' + str(individual) + '_subckt.cir' },
				'tb':      { 'file': 'topdc.cir'}#'testTopCircuit_hspice.cir' }, 
				}, 
			'settings':{
				'debug': 0
				  },
			'params':{
				'vdd' : 15.0,
				'temperature': 25,
				'rl': 10e6,
				'iref': 1e-2
				  },
			'options': {
				'method': 'trap', 
				#'noautoconv': True,		#prevent Spice to try to converge for every price
				  }
			}
		}
			
	analyses={
		'op': {
			'head': 'opus', 
			'modules': [ 'def', 'tb' ],  
			'params': {'rl':10e6
			}, 
			'command': "op()"
		}, 
		'vdd_sweep_t1': {
			'head': 'opus',
			'modules': [ 'def', 'tb' ], 
			'params': {'rl':10e6,
				   'temperature': -20,
			}, 
			'saves': [ ], 
			'command': "dc(1.8, 20, 'lin', 100, 'vdd', 'dc')" #if stop point changed, change also in reproduction cost function. 
			}, 
		'vdd_sweep_t2': {
			'head': 'opus',
			'modules': [ 'def', 'tb' ], 
			'params': {'rl':10e6,
				   'temperature': 25,
			}, 
			'saves': [ ], 
			'command': "dc(1.8, 20, 'lin', 100, 'vdd', 'dc')" #if stop point changed, change also in reproduction cost function. 
			}, 
		'vdd_sweep_t3': {
			'head': 'opus',
			'modules': [ 'def', 'tb' ], 
			'params': {'rl':10e6,
				   'temperature': 120,
			}, 
			'saves': [ ], 
			'command': "dc(1.8, 20, 'lin', 100, 'vdd', 'dc')" #if stop point changed, change also in reproduction cost function. 
			},
		'vdd_sweep_r1': {
			'head': 'opus',
			'modules': [ 'def', 'tb' ], 
			'params': {'rl':10e6,
				   'temperature': 25,
			}, 
			'saves': [ ], 
			'command': "dc(1.8, 20, 'lin', 100, 'vdd', 'dc')" #if stop point changed, change also in reproduction cost function. 
			}, 
		'vdd_sweep_r2': {
			'head': 'opus',
			'modules': [ 'def', 'tb' ], 
			'params': {'rl':10e4,
				   'temperature': 25,
			}, 
			'saves': [ ], 
			'command': "dc(1.8, 20, 'lin', 100, 'vdd', 'dc')" #if stop point changed, change also in reproduction cost function. 
			}, 
		'vdd_sweep_r3': {
			'head': 'opus',
			'modules': [ 'def', 'tb' ], 
			'params': {'rl':10e2,
				   'temperature': 25,
			}, 
			'saves': [ ], 
			'command': "dc(1.8, 20, 'lin', 100, 'vdd', 'dc')" #if stop point changed, change also in reproduction cost function. 
			}, 
		 }
			
	measures={
		'vout_vdd_temp1':{
			'analysis' : 'vdd_sweep_t1',
			'corners' : [ 'nominal' ],
			'expression': 'v("vout")',
			'vector' : True
		 },
		'vout_vdd_temp2':{
			'analysis' : 'vdd_sweep_t2',
			'corners' : [ 'nominal' ],
			'expression': 'v("vout")',
			'vector' : True
		 },
		'vout_vdd_temp3':{
			'analysis' : 'vdd_sweep_t3',
			'corners' : [ 'nominal' ],
			'expression': 'v("vout")',
			'vector' : True
		 },
		'vout_vdd_scale' : {
			'analysis' : 'vdd_sweep_t1',
			'corners' : [ 'nominal' ],
			'script': """__result = scale()""",
			'vector' : True,
		},
		'vout_vdd_res1':{
			'analysis' : 'vdd_sweep_r1',
			'corners' : [ 'nominal' ],
			'expression': 'v("vout")',
			'vector' : True
		 },
		'vout_vdd_res2':{
			'analysis' : 'vdd_sweep_r2',
			'corners' : [ 'nominal' ],
			'expression': 'v("vout")',
			'vector' : True
		 },
		'vout_vdd_res3':{
			'analysis' : 'vdd_sweep_r3',
			'corners' : [ 'nominal' ],
			'expression': 'v("vout")',
			'vector' : True
		 },
		'vout_vdd_res_scale' : {
			'analysis' : 'vdd_sweep_r1',
			'corners' : [ 'nominal' ],
			'script': """__result = scale()""",
			'vector' : True,
		},
		
		'power':{
			'analysis' : 'op',
			'corners' : [ 'nominal' ],
			'expression': '-v("vdd")*i("vdd")',
			'vector' : False
		 },
	}
		
	corners = {
		'nominal': {
			'params': {
				'temperature': 25, 
				'vcc': 15, 
				  }
			    }
		}		
	pe=PerformanceEvaluator(heads, analyses, measures, corners, debug=0)

	results,ancount=pe()
	pe.finalize()  
	
	return results


def evaluateVoltageRef(generation, individual):
	"""Evaluation script for pyopus. It evaluates performance of an voltage reference. """
	heads={
		'opus': {
			'simulator' : 'SpiceOpus',
			'moddefs': {
				'def':     { 'file': 'g_' + str(generation) + '_i_' + str(individual) + '_subckt.cir' },
				'tb':      { 'file': 'topdc.cir'}#'testTopCircuit_hspice.cir' }, 
				}, 
			'settings':{
				'debug': 0
				  },
			'params':{
				'vdd' : 15.0,
				'temperature': 25,
				'rl': 10e3,
				'iref': 1e-2
				  },
			'options': {
				'method': 'trap', 
				#'noautoconv': True,		#prevent Spice to try to converge for every price
				  }
			}
		}
			
	analyses={
		'op': {
			'head': 'opus', 
			'modules': [ 'def', 'tb' ],  
			'params': {'rl':10e3
			}, 
			'command': "op()"
		}, 
		'vdd_sweep': {
			'head': 'opus',
			'modules': [ 'def', 'tb' ], 
			'params': {'rl':10e3
			}, 
			'saves': [ ], 
			#'command': "dc(3, 20.0, 'lin', 100, 'vdd', 'dc')" #if stop point changed, change also in reproduction cost function. 
			'command': "dc(3, 20, 'lin', 100, 'vdd', 'dc')" #if stop point changed, change also in reproduction cost function. 
			}, 
		'rload_sweep': {
			'head': 'opus',
			'modules': [ 'def', 'tb' ],  
			'command': "dc(10, 100e3, 'dec', 10, 'rl[r]', 'dc')"#if stop point changed, change also in reproduction cost function. 
			},
		'temp_sweep': {
			'head': 'opus',
			'modules': [ 'def', 'tb' ], 
			'params': {'rl':10e3
			}, 
			'options': {
				'method': 'gear'
			},
			'command': "dc(-20, 120, 'lin', 100, '@@temp', 'dc')"#if stop point changed, change also in reproduction cost function. 
			},
		 }
			
		
	
	measures={
		'vout_vdd':{
			'analysis' : 'vdd_sweep',
			'corners' : [ 'nominal' ],
			'expression': 'v("vout")',
			'vector' : True
		 },
		'vout_vdd_scale' : {
			'analysis' : 'vdd_sweep',
			'corners' : [ 'nominal' ],
			'script': """__result = scale()""",
			'vector' : True,
		},
		'vout_rload':{
			'analysis' : 'rload_sweep',
			'corners' : [ 'nominal' ],
			'expression': 'v("vout")',
			'vector' : True
		 },
		'vout_rload_scale' : {
			'analysis' : 'rload_sweep',
			'corners' : [ 'nominal' ],
			'script': """__result = scale()""",
			'vector' : True,
		},
		'vout_temp':{
			'analysis' : 'temp_sweep',
			'corners' : [ 'nominal' ],
			'expression': 'v("vout")',
			'vector' : True
		 },
		'vout_temp_scale' : {
			'analysis' : 'temp_sweep',
			'corners' : [ 'nominal' ],
			'script': """__result = scale()""",
			'vector' : True,
		},
		'power':{
			'analysis' : 'op',
			'corners' : [ 'nominal' ],
			'expression': '-v("vdd")*i("vdd")',
			'vector' : False
		 },
	}
		
	corners = {
		'nominal': {
			'params': {
				'temperature': 25, 
				'vcc': 15, 
				  }
			    }
		}		
	pe=PerformanceEvaluator(heads, analyses, measures, corners, debug=0)

	results,ancount=pe()
	pe.finalize()  
	
	return results
      

def evaluatePassiveFilter_SUHAD(generation, individual):
	heads={
		'opus': {
			'simulator' : 'HSpice',
			'moddefs': {
				'def':     { 'file': 'g_' + str(generation) + '_i_' + str(individual) + '_subckt.cir' },
				'tb':      { 'file': 'testTopCircuit_hspice.cir' }, 
			}, 
			'settings':{
				'debug': 0
			},
			'params':{
				'vcc' : 15.0,
				'temperature': 25
			},
			'options': {
				'method': 'trap', 
				'noautoconv': True,		#prevent Spice to try to converge for every price
				'maxord': 2,
				'trapratio' : 10,
				'integdebug': True
			}
		}
	}

	analyses={
		'ac' : {
			'head':'opus',
			'modules': [ 'def', 'tb' ], 
			'command': "ac(1, 100e3, 'dec', 100)"	
		},
		#'tran':{
		#	'head':'opus',
		#	'modules': [ 'def', 'tb' ],
		#	#'command': "tran(1e-6, 200e-3, 150e-3)"
		#	'command': "tran(1e-6, 25e-3, 5e-3)"
		#},
	}

	measures={
		#'vout':{
		#	'analysis' : 'tran',
		#	'corners': ['nominal'],
		#	'expression': 'v("vout")',
		#	'vector' : True,
		# },
		#'t':{
		#	'analysis' : 'tran',
		#	'corners': ['nominal'],
		#	'expression': 'v("time")',
		#	'vector' : True,
		# },
		#'vinT':{
		#	'analysis' : 'tran',
		#	'corners': ['nominal'],
		#	'expression': 'v("vin")',
		#	'vector' : True,
		#},
		'x' : {
			'analysis' : 'ac',
			'corners' : [ 'nominal' ],
			'expression': 'scale()',
			'vector' : True
		},
		'y' : {
			'analysis' : 'ac',
			'corners' : [ 'nominal' ],
			#'expression': 'm.ACmag(m.ACtf(v("16"), v("1")), "db20")',
			'expression': 'm.ACmag(m.ACtf(v("vout"), v("vin")), "db20")',
			'vector': True
			},
		'ripple':{
			'analysis':None,
			'corners' : [ 'nominal' ],
			'script': """
x=result['x'][thisCorner]
y=result['y'][thisCorner]
i=m.IatXval(x,1e3)
y1=m.XatIrange(y,0,i)
rup=y1.max()

i=m.IatXval(x,0.75e3)
y1=m.XatIrange(y,0,i)
rdn=y1.min()

nom=y[0]

__result=np.maximum(rup-nom, nom-rdn)
"""
		},
		'gain':{
			'analysis': 'ac',
			'corners' : [ 'nominal' ],
			'expression': 'm.ACmag(m.ACtf(v("vout"), v("vin")), "db20")[0]'
		},
		'damping':{
			'analysis':None,
			'corners' : [ 'nominal' ],
			'script': """
x=result['x'][thisCorner]
y=result['y'][thisCorner]
i=m.IatXval(x,1.35e3)
y1=m.XatIrange(y,i,y.size-1)
__result=y[0] - y1.max()
"""
		},
	}
		
	corners = {
		'nominal': {
			'params': {
				'temperature': 25, 
				'vcc': 15, 
			}
		}
	}		
	pe=PerformanceEvaluator(heads, analyses, measures, corners, debug=0)

	results,ancount=pe()
	pe.finalize()
	return results


def evaluateActiveFilter_SUHAD(generation, individual):
	heads={
		'opus': {
			'simulator' : 'HSpice',
			'moddefs': {
				'def':     { 'file': 'g_' + str(generation) + '_i_' + str(individual) + '_subckt.cir' },
				'tb':      { 'file': 'testTopCircuit_hspice.cir' }, 
			}, 
			'settings':{
				'debug': 0
			},
			'params':{
				'vcc' : 15.0,
				'temperature': 25
			},
			'options': {
				'method': 'trap', 
				'noautoconv': True,		#prevent Spice to try to converge for every price
				'maxord': 2,
				'trapratio' : 10,
				'integdebug': True
			}
		}
	}

	analyses={
		'ac' : {
			'head':'opus',
			'modules': [ 'def', 'tb' ], 
			'command': "ac(1, 100e3, 'dec', 100)"	
		},
		'tran':{
			'head':'opus',
			'modules': [ 'def', 'tb' ],
			#'command': "tran(1e-6, 200e-3, 150e-3)"
			'command': "tran(1e-6, 25e-3, 5e-3)"
		},
	}

	measures={
		'vout':{
			'analysis' : 'tran',
			'corners': ['nominal'],
			'expression': 'v("vout")',
			'vector' : True,
		 },
		't':{
			'analysis' : 'tran',
			'corners': ['nominal'],
			'expression': 'scale()',
			'vector' : True,
		 },
		'vinT':{
			'analysis' : 'tran',
			'corners': ['nominal'],
			'expression': 'v("vin")',
			'vector' : True,
		 },
		'y' : {
			'analysis' : 'ac',
			'corners' : [ 'nominal' ],
			'expression': 'm.ACmag(m.ACtf(v("vout"), v("vin")), "db20")',
			'vector': True
			},		
		'x' : {
			'analysis' : 'ac',
			'corners' : [ 'nominal' ],
			'script': """__result = scale()""",
			'vector' : True,
		},
		'isOutVNonStationary':{		#stationary output detection (added 16.5.2016)
			'analysis': None,
			'corners' : ['nominal'],
			'script': """
voutMax = max(result['vout']['nominal'])
voutMin = min(result['vout']['nominal'])
if abs(voutMax-voutMin) < 0.5:
  __result = False
else:
  __result = True
""",
			'vector': False
		},
		
		
		'THD':{
			'analysis' : None,
			'corners': ['nominal'],
			'script':"""
t = result['t']['nominal']
vout = result['vout']['nominal']
sigLength = max(t)-min(t)

t_inter = np.linspace(min(t), max(t), 1024*8) #interpolation for fft
vout_inter = np.interp(t_inter, t, vout)

n = len(vout)           
k = np.arange(n)
T = n/len(vout)
frq = k/T/sigLength # two sides frequency range
freq = frq[range(n/2)]           # one side frequency range
Y = np.fft.fft(vout_inter)/n              # fft computing and normalization
Y = Y[range(n/2)]
absY = abs(Y)
sortedI = np.argsort(absY)
highH = absY[sortedI[-10:-1]]
baseH = absY[sortedI[-1]]
THD_ = np.sqrt(np.sum(np.power(highH,2)))/baseH*100
__result = THD_
"""
		 },
		'ripple':{
			'analysis':None,
			'corners' : [ 'nominal' ],
			'script': """
x=result['x'][thisCorner]
y=result['y'][thisCorner]
i=m.IatXval(x,1e3)
y1=m.XatIrange(y,0,i)
rup=y1.max()

i=m.IatXval(x,0.75e3)
y1=m.XatIrange(y,0,i)
rdn=y1.min()

nom=y[0]

__result=np.maximum(rup-nom, nom-rdn)
"""
		},
		'gain':{
			'analysis': 'ac',
			'corners' : [ 'nominal' ],
			'expression': 'm.ACmag(m.ACtf(v("vout"), v("vin")), "db20")[0]'
		},
		'damping':{
			'analysis':None,
			'corners' : [ 'nominal' ],
			'script': """
x=result['x'][thisCorner]
y=result['y'][thisCorner]
i=m.IatXval(x,1.35e3)
y1=m.XatIrange(y,i,y.size-1)
__result=y[0] - y1.max()
"""
		},
		'xb1' : {
			'analysis' : 'ac',
			'corners' : [ 'nominal' ],
			'expression': '[np.real(scale())[0] , 1e3]',
			'vector' : True,
		},
		'yb1' : {
			'analysis' : 'ac',
			'corners' : [ 'nominal' ],
			'expression': '[20.5, 20.5]',
			'vector' : True,
		},
		'xb2' : {
			'analysis' : 'ac',
			'corners' : [ 'nominal' ],
			'expression': '[np.real(scale())[0] , 750]',
			'vector' : True,
		},
		'yb2' : {
			'analysis' : 'ac',
			'corners' : [ 'nominal' ],
			'expression': '[19.5, 19.5]',
			'vector' : True,
		},
		'xb3' : {
			'analysis' : 'ac',
			'corners' : [ 'nominal' ],
			'expression': '[1.35e3, np.real(scale())[-1]]',
			'vector' : True,
		},
		'yb3' : {
			'analysis' : 'ac',
			'corners' : [ 'nominal' ],
			'expression': '[-40, -40]',
			'vector' : True,
		}
	}
		
	corners = {
		'nominal': {
			'params': {
				'temperature': 25, 
				'vcc': 15, 
			}
		}
	}		
	pe=PerformanceEvaluator(heads, analyses, measures, corners, debug=0)

	results,ancount=pe()
	pe.finalize()
	return results



def evaluatePassiveFilter(generation, individual):
	"""Evaluates passive filter performance."""
	#print "GEN-%s-IND-%s's evaluation started..." %(generation, individual)
	heads = {
		'opus': {
			'simulator': 'SpiceOpus', 
			'settings': {
				'debug': 0
			}, 
			'moddefs': {
				'def':     { 'file': 'g_' + str(generation) + '_i_' + str(individual) + '_subckt.cir' },
				'tb':      { 'file': 'testTopCircuit_pLP.cir' }, 
			}, 
			'options': {
				'method': 'trap', 
				'noautoconv': True		#prevent Spice to try to converge for every price
			}, 
			'params': {#no params
			}
		}
	}  
	analyses = {
		'op': {
			'head': 'opus', 
			'modules': [ 'def', 'tb' ], 
			'options': {
				'method': 'gear'
			}, 
			'params': {
			},
			# Save current through vdd. Save voltages at inp, inn, and out. 
			# Save vgs, vth, vds, and vdsat for mn2, mn3, and mn1. 
			'saves': [ 
			],  
			'command': "op()"
		}, 
		'ac': {
			'head': 'opus', 
			'modules': [ 'def', 'tb' ], 
			'options': {
				'method': 'gear'
			},
			'params': {
			},	
			'saves': [ ], 
			'command': "ac(1.0, 1e5, 'dec', 100)"
		}
	}
	measures = {
		'tf':{
			'analysis':'ac',
			'expression': "m.ACmag(m.ACtf(v('vout'), v('vin')), unit='abs')",
			'vector': True
		},
		'freq_scale':{
			'analysis':'ac',
			'expression': "np.real(v('frequency'))",
			'vector':True
		},
		'cutoff':{
			'analysis':'ac',
			'expression': "m.ACbandwidth(m.ACtf(v('vout'), v('vin')), v('frequency'), filter='lp', levelType='db', level=-60.0)"
		},
		'bw':{
			'analysis':'ac',
			'expression': "m.ACbandwidth(m.ACtf(v('vout'), v('vin')), v('frequency'), filter='lp', levelType='db', level=-3.0)"
		}
	}
	pe=PerformanceEvaluator(heads, analyses, measures, debug=0)

	results,ancount=pe()
	pe.finalize()
	#print "GEN-%s-IND-%s's evaluation finished." %(generation, individual)

	if results['freq_scale']['default']!=None:
	  import numpy as np
	  idealFilter = np.append(np.ones(len(results['freq_scale']['default'])-183), np.zeros(183)*(0))
	  results['idealFilter'] = idealFilter

	
	#plotting:
	"""
	plot = False
	if plot == True:
	  pyopl.init()
	  fig=pyopl.figure(windowTitle="Figure - single axes", figpx=(600,400), dpi=100)  
	  pyopl.lock(True)
	  if pyopl.alive(fig):
		  x=results['freq_scale']['default']
		  y=results['tf']['default']
		  y2=idealFilter
		  ax=fig.add_axes((0.12,0.12,0.76,0.76))
		  ax.semilogx(x, y, '-', color=(1,0,0))
		  ax.semilogx(x, y2, '-', color=(0,0,1))
		  pyopl.draw(fig)
	  pyopl.lock(False)
	  #pyopl.join()
	  pyopl.saveFigure(fig, "best.png")
	"""
	return results

def evaluateDiffStage(generation, individual):
	print "GEN-%s-IND-%s's evaluation started..." %(generation, individual)
	heads = {
		'opus': {
			'simulator': 'SpiceOpus', 
			'settings': {
				'debug': 0
			}, 
			'moddefs': {
				#'def':     { 'file': 'HOT_CIRCUIT_subckt.cir' }, 
				'def':     { 'file': 'G_' + str(generation) + '_I_' + str(individual) + '_subckt.cir' },
				'tb':      { 'file': 'testTopCircuit.cir' }, 
			}, 
			'options': {
				'method': 'trap', 
				'noautoconv': True		#prevent Spice to try to converge for every price
			}, 
			'params': {
				'acdif': 0.0, 
				'accom': 0.0
			}
		}
	}
	
	analyses = {
		'op': {
			'head': 'opus', 
			'modules': [ 'def', 'tb' ], 
			'options': {
				'method': 'gear'
			}, 
			'params': {
			},
			# Save current through vdd. Save voltages at inp, inn, and out. 
			# Save vgs, vth, vds, and vdsat for mn2, mn3, and mn1. 
			'saves': [ 
			],  
			'command': "op()"
		}, 
		'dc': {
			'head': 'opus', 
			'modules': [ 'def', 'tb' ], 
			'options': {
				'method': 'gear'
			},
			'params': {
			},	
			'saves': [ ], 
			'command': "dc(-1, 1, 'lin', 100, 'vd', 'dc')"
		},
		'acdif': {
			'head': 'opus', 
			'modules': [ 'def', 'tb' ], 
			'options': {
				'method': 'gear'
			},
			'params': {
				'acdif': 1.0,
				'accom': 0.0
			},	
			'saves': [ ], 
			'command': "ac(1.0, 1e9, 'dec', 100)"
		}, 
		'accom': {
			'head': 'opus', 
			'modules': [ 'def', 'tb' ], 
			'options': {
				'method': 'gear'
			},
			'params': {
				'acdif': 0.0,
				'accom': 1.0
			},	
			'saves': [ ], 
			'command': "ac(1.0, 1e9, 'dec', 100)"
		}
	}

	measures = {
		'isup': {
			'analysis': 'op', 
			'script': "__result=-i('vp')"
		}, 
		'isupN': {
			'analysis': 'op', 
			'script': "__result=-i('vn')"
		}, 

		'swing': {
			'analysis': 'dc', 
			'expression': "m.DCswingAtGain(v('vout1'), v('vin1', 'vin2'), 0.5, 'out')"
		},
		'swing2': {
			'analysis': 'dc', 
			'expression': "m.DCswingAtGain(v('vout2'), v('vin1', 'vin2'), 0.5, 'out')"
		},
		#'dcin': {
		#	'analysis': 'dc', 
		#	'expression': "v('vin1', 'vin2')", 
		#	'vector': True
		#},
		#'dcout': {
		#	'analysis': 'dc', 
		#	'expression': "v('vout1')", 
		#	'vector': True
		#},
		'difgain': {
			'analysis': 'acdif', 
			'expression': "m.ACgain(m.ACtf(v('vout1'), v('vin1', 'vin2')))"
		},
		'difgain2': {
			'analysis': 'acdif', 
			'expression': "m.ACgain(m.ACtf(v('vout2'), v('vin1', 'vin2')))"
		},
		'cmgain': {
			'analysis': 'accom', 
			'expression': "m.ACgain(m.ACtf(v('vout1'), v('com')))"
		},
		'cmgain2': {
			'analysis': 'accom', 
			'expression': "m.ACgain(m.ACtf(v('vout2'), v('com')))"
		},
		#'cmrr': {
		#	'analysis': None, 
		#	'expression': "result['difgain'][thisCorner]-result['cmgain'][thisCorner]", 
		#	'depends': [ 'difgain', 'cmgain' ]
		#},
	}

	pe=PerformanceEvaluator(heads, analyses, measures, debug=0)

	results,ancount=pe()
	
	#print "isup:", results['isup']['default']
	#print "swing:", results['swing']['default']
	#print "CMRR:", results['cmrr']['default']
	#print ancount

	#print("")
	#print(pe.formatResults(nMeasureName=10, nCornerName=15)

	#plotting:
	"""
	pyopl.init()
	fig=pyopl.figure(windowTitle="Figure - single axes", figpx=(600,400), dpi=100)  
	pyopl.lock(True)
	if pyopl.alive(fig):
		x=results['dcin']['default']
		y=results['dcout']['default']
		ax=fig.add_axes((0.12,0.12,0.76,0.76))
		ax.plot(x, y, '-', color=(1,0,0))
		pyopl.draw(fig)
	pyopl.lock(False)
	
	pyopl.join()
	"""
	pe.finalize()
	print "GEN-%s-IND-%s's evaluation finished." %(generation, individual)
	print
	return results
	

def evaluatePassiveFilter_2(generation, individual):
	heads={
		'opus': {
			'simulator' : 'HSpice',
			'moddefs': {
				'def':     { 'file': 'g_' + str(generation) + '_i_' + str(individual) + '_subckt.cir' },
				'tb':      { 'file': 'testTopCircuit_hspice.cir' }, 
			}, 
			'settings':{
				'debug': 0
			},
			'params':{
				'vcc' : 15.0,
				'temperature': 25
			},
			'options': {
				'method': 'trap', 
				'noautoconv': True,		#prevent Spice to try to converge for every price
				'maxord': 2,
				'trapratio' : 10,
				'integdebug': True
			}
		}
	}

	analyses={
		'ac' : {
			'head':'opus',
			'modules': [ 'def', 'tb' ], 
			'command': "ac(1, 100e3, 'dec', 100)"	
		},
		#'tran':{
		#	'head':'opus',
		#	'modules': [ 'def', 'tb' ],
		#	#'command': "tran(1e-6, 200e-3, 150e-3)"
		#	'command': "tran(1e-6, 25e-3, 5e-3)"
		#},
	}

	measures={
		#'vout':{
		#	'analysis' : 'tran',
		#	'corners': ['nominal'],
		#	'expression': 'v("vout")',
		#	'vector' : True,
		# },
		#'t':{
		#	'analysis' : 'tran',
		#	'corners': ['nominal'],
		#	'expression': 'v("time")',
		#	'vector' : True,
		# },
		#'vinT':{
		#	'analysis' : 'tran',
		#	'corners': ['nominal'],
		#	'expression': 'v("vin")',
		#	'vector' : True,
		#},
		'x' : {
			'analysis' : 'ac',
			'corners' : [ 'nominal' ],
			'expression': 'scale()',
			'vector' : True
		},
		'y' : {
			'analysis' : 'ac',
			'corners' : [ 'nominal' ],
			#'expression': 'm.ACmag(m.ACtf(v("16"), v("1")), "db20")',
			'expression': 'm.ACmag(m.ACtf(v("vout"), v("vin")), "db20")',
			'vector': True
			},
		'ripple':{
			'analysis':None,
			'corners' : [ 'nominal' ],
			'script': """
x=result['x'][thisCorner]
y=result['y'][thisCorner]
i=m.IatXval(x,1e3)
y1=m.XatIrange(y,0,i)
rup=y1.max()

i=m.IatXval(x,0.75e3)
y1=m.XatIrange(y,0,i)
rdn=y1.min()

nom=y[0]

__result=np.maximum(rup-nom, nom-rdn)
"""
		},
		'gain':{
			'analysis': 'ac',
			'corners' : [ 'nominal' ],
			'expression': 'm.ACmag(m.ACtf(v("vout"), v("vin")), "db20")[0]'
		},
		'damping':{
			'analysis':None,
			'corners' : [ 'nominal' ],
			'script': """
x=result['x'][thisCorner]
y=result['y'][thisCorner]
i=m.IatXval(x,1.35e3)
y1=m.XatIrange(y,i,y.size-1)
__result=y[0] - y1.max()
"""
		},
		'dumpingSlope':{
			'analysis':None,
			'corners':['nominal'],
			'script':"""
slopesDec = []
for i in range(0, len(result['x'][thisCorner])-100, 10):
  slope = result['y'][thisCorner][i] - result['y'][thisCorner][i+100]
  slopesDec.append(slope)
__result=max(slopesDec)
			  """
		 
		 },
		'bw':{
			'analysis': 'ac',
			'corners': ['nominal'],
			'expression': 'm.ACbandwidth(m.ACtf(v(\'vout\'), v(\'vin\')),scale(), filter="lp", levelType="db", level=-6.0)'
		},
	}
		
	corners = {
		'nominal': {
			'params': {
				'temperature': 25, 
				'vcc': 15, 
			}
		}
	}		
	pe=PerformanceEvaluator(heads, analyses, measures, corners, debug=0)

	results,ancount=pe()
	pe.finalize()
	return results


def evaluateActiveFilter_2(generation, individual):
	heads={
		'hspice1': {
			'simulator' : 'HSpice',
			'moddefs': {
				'def':     { 'file': 'g_' + str(generation) + '_i_' + str(individual) + '_subckt.cir' },
				'tb':      { 'file': 'testTopCirc_hspice_AF.cir' },
				'tb_rout': { 'file': 'testTopCirc_hspice_AF_outimped.cir'}
			}, 
			'settings':{
				'debug': 0
			},
			'params':{
				'vcc' : 15.0,
				'temperature': 25
			},
			'options': {
				'method': 'trap', 
				#'noautoconv': True,		#prevent Spice to try to converge for every price
				#'maxord': 2,
				#'trapratio' : 10,
				#'integdebug': True
			}
		}
	}

	analyses={
		'ac' : {
			'head':'hspice1',
			'modules': [ 'def', 'tb' ], 
			'command': "ac(1, 100e3, 'dec', 100)"	
		},
		'tran_Lf':{#-------------------tran_Lf (at Low Frequency)
			'head':'hspice1',
			'modules': [ 'def', 'tb' ],
			'params':{'frequ':100,
				  'rstopin':1e-6
			},
			#'command': "tran(1e-6, 125e-3, 105e-3)"
			'command': "tran("+str(1.0/100/160)+", 125e-3, 105e-3)"
		},
		#'tran_Hf':{#-------------------tran_Hf (at High Frequency - near cutoff)
		#	'head':'hspice1',
		#	'modules': [ 'def', 'tb' ],
		#	'params':{'frequ':800,
		#		  'rstopin':1e-6
		#	},
		#	'command': "tran(1e-6, 15e-3, 10e-3)"
		#},
		'ac_outImp':{#------------------Output impedance measurement (ac analysis)
			'head': 'hspice1', 
			'modules': [ 'def', 'tb_rout' ],  
			'command': "ac(1, 100e3, 'dec', 100)"
		},
	}

	measures={
		'vout_Lf':{	#-------------------tranLf
			'analysis' : 'tran_Lf',
			'corners': ['nominal'],
			'expression': 'v("vout")',
			'vector' : True,
		 },
		't_Lf':{
			'analysis' : 'tran_Lf',
			'corners': ['nominal'],
			'expression': 'scale()',
			'vector' : True,
		 },
		'vinT_Lf':{
			'analysis' : 'tran_Lf',
			'corners': ['nominal'],
			'expression': 'v("vin")',
			'vector' : True,
		 },
		#'vout_Hf':{	#-------------------tranHf
		#	'analysis' : 'tran_Hf',
		#	'corners': ['nominal'],
		#	'expression': 'v("vout")',
		#	'vector' : True,
		# },
		#'t_Hf':{
		#	'analysis' : 'tran_Hf',
		#	'corners': ['nominal'],
		#	'expression': 'scale()',
		#	'vector' : True,
		# },
		#'vinT_Hf':{
		#	'analysis' : 'tran_Hf',
		#	'corners': ['nominal'],
		#	'expression': 'v("vin")',
		#	'vector' : True,
		#},
		'y' : {		#-------------------AC analises and bode stuff
			'analysis' : 'ac', 
			'corners' : [ 'nominal' ],
			'expression': 'm.ACmag(m.ACtf(v("vout"), v("vin")), "db20")',
			'vector': True
			},		
		'x' : {
			'analysis' : 'ac',
			'corners' : [ 'nominal' ],
			'script': """__result = scale()""",
			'vector' : True,
		},
		
		'inimped' : {	#-------------------input impedance mesurement
			'analysis' : 'ac', 
			'corners' : [ 'nominal' ],
			'expression': 'max(m.ACmag(m.ACtf(v("vin"), i("vd")), "abs"))',
			'vector': False
			},	
		
		'outimped' : {	#-------------------output impedance mesurement
			'analysis' : 'ac_outImp', 
			'corners' : [ 'nominal' ],
			'expression': 'max(m.ACmag(m.ACtf(v("vout"), i("vload")), "abs"))',
			'vector': False
			},	
		
		'isOutVNonStationary':{		#stationary output detection (added 16.5.2016)
			'analysis': None,	#Added because lowest THD was if circuit was dead (or output was flat) of course.
			'corners' : ['nominal'],
			'script': """
voutMax = max(result['vout_Lf']['nominal'])
voutMin = min(result['vout_Lf']['nominal'])
if abs(voutMax-voutMin) < 1.5:
  __result = False
else:
  __result = True
""",
			'vector': False
		},
		
				#-------------------THD calculations
		'THD_Lf':{		#-------------------THD at low frequency
			'analysis' : None,
			'corners': ['nominal'],
			'script':"""
t = result['t_Lf']['nominal']
vout = result['vout_Lf']['nominal']
sigLength = max(t)-min(t)

t_inter = np.linspace(min(t), max(t), 1024*8) #interpolation for fft
vout_inter = np.interp(t_inter, t, vout)

n = len(vout)           
k = np.arange(n)
T = n/len(vout)
frq = k/T/sigLength # two sides frequency range
freq = frq[range(n/2)]           # one side frequency range
Y = np.fft.fft(vout_inter)/n              # fft computing and normalization
Y = Y[range(n/2)]
absY = abs(Y)
sortedI = np.argsort(absY)
highH = absY[sortedI[-10:-1]]
baseH = absY[sortedI[-1]]
THD_ = np.sqrt(np.sum(np.power(highH,2)))/baseH*100
__result = THD_
"""
		 },
		
#		'THD_Hf':{		#-------------------THD at high frequency
#			'analysis' : None,
#			'corners': ['nominal'],
#			'script':"""
#t = result['t_Hf']['nominal']
#vout = result['vout_Hf']['nominal']
#sigLength = max(t)-min(t)
#
#t_inter = np.linspace(min(t), max(t), 1024*8) #interpolation for fft
#vout_inter = np.interp(t_inter, t, vout)
#
#n = len(vout)           
#k = np.arange(n)
#T = n/len(vout)
#frq = k/T/sigLength # two sides frequency range
#freq = frq[range(n/2)]           # one side frequency range
#Y = np.fft.fft(vout_inter)/n              # fft computing and normalization
#Y = Y[range(n/2)]
#absY = abs(Y)
#sortedI = np.argsort(absY)
#highH = absY[sortedI[-10:-1]]
#baseH = absY[sortedI[-1]]
#THD_ = np.sqrt(np.sum(np.power(highH,2)))/baseH*100
#__result = THD_
#"""
#		 },
		
		
		
		'ripple':{
			'analysis':None,
			'corners' : [ 'nominal' ],
			'script': """
bw = result['bw'][thisCorner]
x=result['x'][thisCorner]
y=result['y'][thisCorner]
#i=m.IatXval(x,1e3)
i=m.IatXval(x,bw*1.2)
y1=m.XatIrange(y,0,i)
rup=y1.max()

#i=m.IatXval(x,0.75e3)
i=m.IatXval(x,bw*0.75)
y1=m.XatIrange(y,0,i)
rdn=y1.min()

nom=y[0]

__result=np.maximum(rup-nom, nom-rdn)
"""
		},
		'gain':{
			'analysis': 'ac',
			'corners' : [ 'nominal' ],
			'expression': 'm.ACmag(m.ACtf(v("vout"), v("vin")), "db20")[0]'
		},
		'damping':{
			'analysis':None,
			'corners' : [ 'nominal' ],
			'script': """
bw = result['bw'][thisCorner]
x=result['x'][thisCorner]
y=result['y'][thisCorner]

i0=m.IatXval(x,bw*0.75)
y0=m.XatIrange(y,0,i0)

i=m.IatXval(x,bw*1.35)
y1=m.XatIrange(y,i,y.size-1)
#__result=y[0] - y1.max()
__result=np.average(y0) - y1.max()
"""
		},
		'bw':{
			'analysis': 'ac',
			'corners': ['nominal'],
			'expression': 'm.ACbandwidth(m.ACtf(v(\'vout\'), v(\'vin\')),scale(), filter="lp", levelType="db", level=-6.0)'
		},

		'dampingSlopes':{
			'analysis':None,
			'corners':['nominal'],
			'vector' : True,
			'script':"""
slopesDec = []
for i in range(0, len(result['x'][thisCorner])-100, 10):
  slope = result['y'][thisCorner][i] - result['y'][thisCorner][i+100]
  slopesDec.append(slope)
__result= slopesDec
			  """
		 },
		
		'maxDampingSlopeFreq':{
			'analysis':None,
			'corners':['nominal'],
			'vector' : False,
			'script':"""
FreqI = np.argmax(result['dampingSlopes'][thisCorner])
dampingSlopeFreq = result['x'][thisCorner][FreqI*10]
__result= dampingSlopeFreq
			  """
		 },		
		

		'maxDampingSlope':{
			'analysis':None,
			'corners':['nominal'],
			'vector' : True,
			'script':"""
slopesDec = []
for i in range(0, len(result['x'][thisCorner])-100, 10):
  slope = result['y'][thisCorner][i] - result['y'][thisCorner][i+100]
  slopesDec.append(slope)
FreqI = np.argmax(slopesDec)
dampingSlopeFreq = result['x'][thisCorner][FreqI*10]  
__result= [ max(slopesDec), dampingSlopeFreq]
			  """
		 },

		'is_LP':{
			'analysis':None,
			'corners': ['nominal'],
			'vector' : False,
			'script':"""
gain0 = result['y'][thisCorner][0]
gainEnd = result['y'][thisCorner][-1]
__result=(np.power(10,float(gain0)/20)-np.power(10,float(gainEnd)/20))
""",		
		}	
	}
		
	corners = {
		'nominal': {
			'params': {
				'temperature': 25, 
				'vcc': 15, 
			}
		}
		
	}		
	pe=PerformanceEvaluator(heads, analyses, measures, corners, debug=0)

	results,ancount=pe()
	pe.finalize()
	return results


def evaluateActiveFilter_2_Hf(generation, individual, Hf):
	heads={
		'opus': {
			'simulator' : 'HSpice',
			'moddefs': {
				'def':     { 'file': 'g_' + str(generation) + '_i_' + str(individual) + '_subckt.cir' },
				'tb':      { 'file': 'testTopCirc_hspice_AF.cir' }, 
			}, 
			'settings':{
				'debug': 0
			},
			'params':{
				'vcc' : 15.0,
				'temperature': 25
			},
			'options': {
				'method': 'gear', 
				#'noautoconv': True,		#prevent Spice to try to converge for every price
				#'maxord': 2,
				#'trapratio' : 10,
				#'integdebug': True
			}
		}
	}
	
	
	frequHf = Hf*0.8
	start = 10e-3
	analyses={
		'tran_Hf':{#-------------------tran_Hf (at High Frequency - near cutoff)
			'head':'opus',
			'modules': [ 'def', 'tb' ],
			'params':{'frequ':frequHf,
				  'rstopin':1e-6
			},
			#		step			stop (start+t0x4 periode)	start (fixed)
			'command': "tran("+str(1/frequHf/160) +","+ str(start+1/frequHf*1) + ", " + str(start) + ")"
		}
	}

	measures={
		'vout_Hf':{	#-------------------tranHf
			'analysis' : 'tran_Hf',
			'corners': ['nominal'],
			'expression': 'v("vout")',
			'vector' : True,
		 },
		't_Hf':{
			'analysis' : 'tran_Hf',
			'corners': ['nominal'],
			'expression': 'scale()',
			'vector' : True,
		 },
		'vinT_Hf':{
			'analysis' : 'tran_Hf',
			'corners': ['nominal'],
			'expression': 'v("vin")',
			'vector' : True,
		 },
		
		'isOutVNonStationary':{		#stationary output detection (added 16.5.2016)
			'analysis': None,	#Added because lowest THD was if circuit was dead (or output was flat) of course.
			'corners' : ['nominal'],
			'script': """
voutMax = max(result['vout_Lf']['nominal'])
voutMin = min(result['vout_Lf']['nominal'])
if abs(voutMax-voutMin) < 1.5:
  __result = False
else:
  __result = True
""",
			'vector': False
		},
		
				#-------------------THD calculations		
		'THD_Hf':{		#-------------------THD at high frequency
			'analysis' : None,
			'corners': ['nominal'],
			'script':"""
t = result['t_Hf']['nominal']
vout = result['vout_Hf']['nominal']
sigLength = max(t)-min(t)

t_inter = np.linspace(min(t), max(t), 1024*8) #interpolation for fft
vout_inter = np.interp(t_inter, t, vout)

n = len(vout)           
k = np.arange(n)
T = n/len(vout)
frq = k/T/sigLength 		# two sides frequency range
freq = frq[range(n/2)]        	# one side frequency range
Y = np.fft.fft(vout_inter)/n	# fft computing and normalization
Y = Y[range(n/2)]
absY = abs(Y)
sortedI = np.argsort(absY)
highH = absY[sortedI[-10:-1]]
baseH = absY[sortedI[-1]]
THD_ = np.sqrt(np.sum(np.power(highH,2)))/baseH*100
__result = THD_
"""
		 },			
	}
		
	corners = {
		'nominal': {
			'params': {
				'temperature': 25, 
				'vcc': 15, 
			}
		}
		
	}	
	#print analyses
	pe=PerformanceEvaluator(heads, analyses, measures, corners, debug=0)

	results,ancount=pe()
	pe.finalize()
	return results
      
      


def evaluateActiveF_2_PB(generation, individual):
  """
  To write. 
  """
  
  heads={
	'hspice1': {
	    'simulator' : 'HSpice',
	    'moddefs': {
		    'def':     { 'file': 'g_' + str(generation) + '_i_' + str(individual) + '_subckt.cir' },
		    'tb':      { 'file': 'testTopCirc_hspice_AF.cir' },
		    'tb_rout': { 'file': 'testTopCirc_hspice_AF_outimped.cir'}
	    }, 
	    'settings':{
		    'debug': 0
	    },
	    'params':{
		    'vcc' : 15.0,
		    'temperature': 25
	    },
	    'options': {
		    'method': 'trap', 
		    #'noautoconv': True,		#prevent Spice to try to converge for every price
		    #'maxord': 2,
		    #'trapratio' : 10,
		    #'integdebug': True
	    }
	  }
	}  

  analyses={
	'ac' : {
		'head':'hspice1',
		'modules': [ 'def', 'tb' ], 
		'command': "ac(1, 100e4, 'dec', 100)"	
	},
	'ac_outImp':{#------------------Output impedance measurement (ac analysis)
		'head': 'hspice1', 
		'modules': [ 'def', 'tb_rout' ],  
		'command': "ac(1, 100e3, 'dec', 100)"
	},
  }
  
  measures={
	'y' : {		#-------------------AC analises and bode stuff
		'analysis' : 'ac', 
		'corners' : [ 'nominal' ],
		'expression': 'm.ACmag(m.ACtf(v("vout"), v("vin")), "db20")',
		'vector': True
		},		
	'x' : {
		'analysis' : 'ac',
		'corners' : [ 'nominal' ],
		'script': """__result = scale()""",
		'vector' : True,
	},
    
	'bw_start':{
		'analysis': 'ac',
		'corners': ['nominal'],
		'expression': 'm.ACbandwidth(m.ACtf(v(\'vout\'), v(\'vin\')),scale(), filter="hp", levelType="db", level=-3.0)'
	}, 
	'bw_stop':{
		'analysis': 'ac',
		'corners': ['nominal'],
		'expression': 'm.ACbandwidth(m.ACtf(v(\'vout\'), v(\'vin\')),scale(), filter="lp", levelType="db", level=-3.0)'
	},  


	'gain':{				
		'analysis':None,		
		'corners' : [ 'nominal' ],
		'script': """
		  
x=result['x'][thisCorner]			
y=result['y'][thisCorner]			
						
bw_start = result['bw_start'][thisCorner]       #measure start and stop freq
if result['bw_start'][thisCorner] == None:      #if freqs not found,  consider whole vector
  bw_start = 1					#find the middle freq between bw_start*1.1 and bw_stop*0.9

bw_stop =  result['bw_stop'][thisCorner]
if result['bw_stop'][thisCorner] == None:
  bw_stop = max(x)

i_start=m.IatXval(x,bw_start*1.1)
i_stop=m.IatXval(x,bw_stop*0.9)

i = m.floor((i_stop-i_start)/2)
i = m.floor(i_start+i)

__result=y[i[0]]
"""
	}, 


	'ripple':{
		'analysis':None,
		'corners' : [ 'nominal' ],
		'script': """
		  
x=result['x'][thisCorner]
y=result['y'][thisCorner]                         #measure start and stop freq
                                                  #if freqs not found,  consider whole vector
bw_start = result['bw_start'][thisCorner]         #check response between bw_start*1.1 and bw_stop*0.9
if result['bw_start'][thisCorner] == None:
  bw_start = 1

bw_stop =  result['bw_stop'][thisCorner]
if result['bw_stop'][thisCorner] == None:
  bw_stop = max(x)

i_start=m.floor(m.IatXval(x,bw_start*1.1))
i_stop=m.floor(m.IatXval(x,bw_stop*0.9))

y1=(m.XatIrange(y,i_start, i_stop)).min()
y2=(m.XatIrange(y,i_start, i_stop)).max()

__result=abs(y1-y2)
"""
	},   
	
	'damping_L':{				#damping at Lower stop band
		'analysis':None,
		'corners' : [ 'nominal' ],
		'script': """
x=result['x'][thisCorner]
y=result['y'][thisCorner]		  

bw_start = result['bw_start'][thisCorner]
if result['bw_start'][thisCorner] == None:
  bw_start = 1

i_start=m.floor(m.IatXval(x,bw_start))	# 100 elements before start freq WARNING! If resolution changed, change!
if (i_start-100)>0:
  i_start-=100
else:
  istart = 0

y1=(m.XatIrange(y,0, i_start)).max()

__result= y1
"""	  
	},
	'damping_H':{				#damping at Higher stop band
		'analysis':None,
		'corners' : [ 'nominal' ],
		'script': """
x=result['x'][thisCorner]
y=result['y'][thisCorner]	

bw_stop = result['bw_stop'][thisCorner]
if result['bw_stop'][thisCorner] == None:
  bw_stop = max(x)

i_stop=m.floor(m.IatXval(x,bw_stop))	## 100 elements after stop freq WARNING! If resolution changed, change!

if (i_stop+100) > len(x):
  i_stop = len(x)-1
else:
  i_stop+=100

y2=(m.XatIrange(y,i_stop, len(y)-1)).max()

__result= y2
"""	  
	},	
	
	
  }
	
  corners = {
	  'nominal': {
		  'params': {
			  'temperature': 25, 
			  'vcc': 15, 
		  }
	  }
	  
  }	
  
  pe=PerformanceEvaluator(heads, analyses, measures, corners, debug=0)

  results,ancount=pe()
  pe.finalize()
  return results



def evaluateActiveF_2_HP(generation, individual):
  """
  Evaluation script for high-pass filter. More to come.
  """
  
  heads={
	'hspice1': {
	    'simulator' : 'HSpice',
	    'moddefs': {
		    'def':     { 'file': 'g_' + str(generation) + '_i_' + str(individual) + '_subckt.cir' },
		    'tb':      { 'file': 'testTopCirc_hspice_AF.cir' },
		    'tb_rout': { 'file': 'testTopCirc_hspice_AF_outimped.cir'}
	    }, 
	    'settings':{
		    'debug': 0
	    },
	    'params':{
		    'vcc' : 15.0,
		    'temperature': 25
	    },
	    'options': {
		    'method': 'trap', 
		    #'noautoconv': True,		#prevent Spice to try to converge for every price
		    #'maxord': 2,
		    #'trapratio' : 10,
		    #'integdebug': True
	    }
	  }
	}    

  analyses={
	'ac' : {
		'head':'hspice1',
		'modules': [ 'def', 'tb' ], 
		'command': "ac(1, 100e4, 'dec', 100)"	
	},
	'ac_outImp':{#------------------Output impedance measurement (ac analysis)
		'head': 'hspice1', 
		'modules': [ 'def', 'tb_rout' ],  
		'command': "ac(1, 100e4, 'dec', 100)"
	},
  }

  measures={
	'y' : {		#-------------------AC analises and bode stuff
		'analysis' : 'ac', 
		'corners' : [ 'nominal' ],
		'expression': 'm.ACmag(m.ACtf(v("vout"), v("vin")), "db20")',
		'vector': True
		},		
	'x' : {
		'analysis' : 'ac',
		'corners' : [ 'nominal' ],
		'script': """__result = scale()""",
		'vector' : True,
	},
	'bw_start':{
		'analysis': 'ac',
		'corners': ['nominal'],
		'expression': 'm.ACbandwidth(m.ACtf(v(\'vout\'), v(\'vin\')),scale(), filter="hp", levelType="db", level=-3.0)'
	}, 
	'gain':{				
		'analysis':None,		
		'corners' : [ 'nominal' ],
		'script': """
x=result['x'][thisCorner]			
y=result['y'][thisCorner]			
						
bw_start = result['bw_start'][thisCorner]       #measure start and stop freq
if result['bw_start'][thisCorner] == None:      #if freqs not found,  consider whole vector
  bw_start = 1					#find the middle freq between bw_start*1.1 and bw_stop*0.9

bw_stop =  max(x)

i_start=m.IatXval(x,bw_start*1.1)
i_stop=m.IatXval(x,bw_stop)

i = m.floor((i_stop-i_start)/2)
i = m.floor(i_start+i)

__result=y[i[0]]
"""
	}, 	
	
	
	'ripple':{
		'analysis':None,
		'corners' : [ 'nominal' ],
		'script': """
		  
x=result['x'][thisCorner]
y=result['y'][thisCorner]                         #measure start freq
                                                  #if freq not found, consider whole vector
bw_start = result['bw_start'][thisCorner]         #bw_stop is at the end
if result['bw_start'][thisCorner] == None:
  bw_start = 1

bw_stop =  max(x)

i_start=m.floor(m.IatXval(x,bw_start*2))
i_stop=m.floor(m.IatXval(x,bw_stop))

y1=(m.XatIrange(y,i_start, i_stop)).min()
#y2=(m.XatIrange(y,i_start, i_stop)).max()
y2=(m.XatIrange(y,1, i_stop)).max() #through the whole range

__result=abs(y1-y2)
"""
	}, 
	'damping_L':{				#damping at Lower stop band
		'analysis':None,
		'corners' : [ 'nominal' ],
		'script': """
x=result['x'][thisCorner]
y=result['y'][thisCorner]		  

bw_start = result['bw_start'][thisCorner]
if result['bw_start'][thisCorner] == None:
  bw_start = 1

i_start=m.floor(m.IatXval(x,bw_start))	# 100 elements before start freq WARNING! If resolution changed, change!
if (i_start-20)>0:
  i_start-=20
else:
  istart = 0

y1=(m.XatIrange(y,0, i_start)).max()

__result= y1
"""	  
	}
  }	
  corners = {
	  'nominal': {
		  'params': {
			  'temperature': 25, 
			  'vcc': 15, 
		  }
	  }
	  
  }	
  
  pe=PerformanceEvaluator(heads, analyses, measures, corners, debug=0)

  results,ancount=pe()
  pe.finalize()
  return results	
      
def evaluateCmosVoltageRef(generation, individual):
	"""Evaluation script for an analog voltage reference. It evaluates performance of an voltage reference. 
	"""
	heads={
		'opus': {
			'simulator' : 'HSpice',
			'moddefs': {
				'def':     { 'file': 'g_' + str(generation) + '_i_' + str(individual) + '_subckt.cir' },
				'tb':      { 'file': 'topdc.cir'} ,
				'tb_psrr': { 'file': 'topdc_psrr.cir'}
				}, 
			'settings':{
				'debug': 0
				  },
			'params':{
				'vdd' : 15.0,
				'temperature': 25,
				'rl': 10e6,
				'iref': 1e-2
				  },
			'options': {
				'method': 'gear', #'method': 'trap' gets stuck often
				'noautoconv': True,		#prevent Spice to try to converge for every price
				  }
			}
		}

	analyses={
		'op': {	#------------------------------------- for power calculation
			'head': 'opus', 
			'modules': [ 'def', 'tb' ],  
			'params': {'rl':10e6
			}, 
			'command': "op()"
		}, 
		'vdd_sweep_t1': {#-------------------------------------vdd sweeps at temperatures and fixed rload
			'head': 'opus',
			'modules': [ 'def', 'tb' ], 
			'params': {'rl':10e6,
				   'temperature': -20,
			}, 
			'saves': [ ], 
			'command': "dc(1.8, 20, 'lin', 100, 'vdd', 'dc')" #if stop point changed, change also in reproduction cost function. 
			}, 
		'vdd_sweep_t2': {
			'head': 'opus',
			'modules': [ 'def', 'tb' ], 
			'params': {'rl':10e6,
				   'temperature': 25,
			}, 
			'saves': [ ], 
			'command': "dc(1.8, 20, 'lin', 100, 'vdd', 'dc')" #if stop point changed, change also in reproduction cost function. 
			}, 
		'vdd_sweep_t3': {
			'head': 'opus',
			'modules': [ 'def', 'tb' ], 
			'params': {'rl':10e6,
				   'temperature': 120,
			}, 
			'saves': [ ], 
			'command': "dc(1.8, 20, 'lin', 100, 'vdd', 'dc')" #if stop point changed, change also in reproduction cost function. 
			},
		'vdd_sweep_r1': {#-------------------------------------vdd sweeps at rloads and fixed temperature (obsolete)
			'head': 'opus',
			'modules': [ 'def', 'tb' ], 
			'params': {'rl':10e6,
				   'temperature': 25,
			}, 
			'saves': [ ], 
			'command': "dc(1.8, 20, 'lin', 100, 'vdd', 'dc')" #if stop point changed, change also in reproduction cost function. 
			}, 
		'vdd_sweep_r2': {
			'head': 'opus',
			'modules': [ 'def', 'tb' ], 
			'params': {'rl':10e4,
				   'temperature': 25,
			}, 
			'saves': [ ], 
			'command': "dc(1.8, 20, 'lin', 100, 'vdd', 'dc')" #if stop point changed, change also in reproduction cost function. 
			}, 
		'vdd_sweep_r3': {
			'head': 'opus',
			'modules': [ 'def', 'tb' ], 
			'params': {'rl':10e2,
				   'temperature': 25,
			}, 
			'saves': [ ], 
			'command': "dc(1.8, 20, 'lin', 100, 'vdd', 'dc')" #if stop point changed, change also in reproduction cost function. 
			}, 
		
		'vdd_psrr_t2':{#-------------------------------------vdd PSRR
			'head': 'opus',
			'modules': [ 'def', 'tb_psrr' ], 
			'params': {'rl':10e6,
				   'temperature': 25,
			}, 
			'saves': [ ], 
			'command': "tran(10e-5, 20e-3, 10e-3)"			
		  },
		 }
			
	measures={
		'vout_vdd_temp1':{
			'analysis' : 'vdd_sweep_t1',
			'corners' : [ 'nominal' ],
			'expression': 'v("vout")',
			'vector' : True
		 },
		'vout_vdd_temp2':{
			'analysis' : 'vdd_sweep_t2',
			'corners' : [ 'nominal' ],
			'expression': 'v("vout")',
			'vector' : True
		 },
		'vout_vdd_temp3':{
			'analysis' : 'vdd_sweep_t3',
			'corners' : [ 'nominal' ],
			'expression': 'v("vout")',
			'vector' : True
		 },
		'vout_vdd_scale' : {
			'analysis' : 'vdd_sweep_t1',
			'corners' : [ 'nominal' ],
			'script': """__result = scale()""",
			'vector' : True,
		},
		'vout_vdd_res1':{
			'analysis' : 'vdd_sweep_r1',
			'corners' : [ 'nominal' ],
			'expression': 'v("vout")',
			'vector' : True
		 },
		'vout_vdd_res2':{
			'analysis' : 'vdd_sweep_r2',
			'corners' : [ 'nominal' ],
			'expression': 'v("vout")',
			'vector' : True
		 },
		'vout_vdd_res3':{
			'analysis' : 'vdd_sweep_r3',
			'corners' : [ 'nominal' ],
			'expression': 'v("vout")',
			'vector' : True
		 },
		'vout_vdd_res_scale' : {
			'analysis' : 'vdd_sweep_r1',
			'corners' : [ 'nominal' ],
			'script': """__result = scale()""",
			'vector' : True,
		},
		
		'power':{
			'analysis' : 'op',
			'corners' : [ 'nominal' ],
			'expression': '-v("vdd")*i("vdd")',
			'vector' : False
		 },
		't' : {
			'analysis' : 'vdd_psrr_t2',
			'corners' : [ 'nominal' ],
			'expression': 'scale()',
			'vector' : True,
		},
		'vin_psrr':{
			'analysis' : 'vdd_psrr_t2',
			'corners': ['nominal'],
			'expression': 'v("vdd")',
			'vector' : True,
		 },		
		'vout_psrr':{
			'analysis' : 'vdd_psrr_t2',
			'corners': ['nominal'],
			'expression': 'v("vout")',
			'vector' : True,
		 },
		'psrr':{
			'analysis' : None,
			'corners' : [ 'nominal' ],
			'vector' : False,
			'script': """
vout = result['vout_psrr']['nominal']
vin = result['vin_psrr']['nominal']
psrr = 10*np.log10(pow(max(vin)-min(vin),2)/pow(max(vout)-min(vout),2))
__result = psrr
"""
		},
		#ADD a measurement of temp sensitivity
		'tempSens':{
			'analysis' : None,
			'corners' : [ 'nominal' ],
			'vector' : False,		
			'script': """
voutt1 = result['vout_vdd_temp1']['nominal']
voutt2 = result['vout_vdd_temp2']['nominal']
voutt3 = result['vout_vdd_temp3']['nominal']
vvdd = result['vout_vdd_scale']['nominal']
i0=int(m.floor(min(m.IatXval(vvdd,11.0))))

voutsens = (voutt3[i0] - voutt1[i0])/(140) #-20...+120 degrees
__result = voutsens

"""
		},
		
		'VatLowVdd':{
			'analysis' : None,
			'corners' : [ 'nominal' ],
			'vector' : False,
			'script': """
vout = result['vout_vdd_temp2']['nominal']
vvdd = result['vout_vdd_scale']['nominal']
i0=int(m.floor(min(m.IatXval(vvdd,2.5, slope='rising'))))
VatVdd = vout[i0]
__result = VatVdd
"""
		},		
	}
		
	corners = {
		'nominal': {
			'params': {
				'temperature': 25, 
				'vcc': 15, 
				  }
			    }
		}		
	pe=PerformanceEvaluator(heads, analyses, measures, corners, debug=0)

	results,ancount=pe()
	pe.finalize()  
	
	return results