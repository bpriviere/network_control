import numpy as np 

param = {}

# time
param['t0'] = 0.
param['tf'] = 10.
param['dt'] = 0.01
param['T'] = np.arange( param.get('t0'), param.get('tf'), param.get('dt'))
param['nt'] = len(param.get('T'))

# number of agents
param['na'] = 4
param['ni'] = param.get('na') # agents plus virtual leader and basis 

# number of dimensions
param['nd'] = 2
param['dof'] = 2*param.get('nd') # degrees of freedom per agent
param['n'] = param.get('dof')*param.get('ni')
param['m'] = param.get('dof')*param.get('ni')

# network
param['R_comm'] = np.inf
param['lambda_a'] = 0.8

# dynamics
param['k1'] = 10
param['k2'] = 10

# desired formation
param['case_xb'] = None # basis for scaling and rotation invariance
param['case_xl'] = None # leader for translation invariance
param['case_xd'] = 'n_node_regular_polygon' # formation

# init
param['plim'] = 2
param['vlim'] = 1 
param['min_dist'] = 0.1

# plotting
param['fn_plots'] = 'plots.pdf'
param['start_marker'] = 's'
param['stop_marker'] = 'x'
param['AgentColor']  = 'b'
param['LeaderColor'] = 'g'
param['BasisColor']  = 'k'