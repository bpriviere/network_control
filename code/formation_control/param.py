import numpy as np 

param = {}

# time
param['t0'] = 0.
param['tf'] = 10.
param['dt'] = 0.05
param['T'] = np.arange( param.get('t0'), param.get('tf'), param.get('dt'))
param['nt'] = len(param.get('T'))

# number of agents
param['na'] = 5
param['ni'] = param.get('na') # agents plus virtual leader and basis 

# number of dimensions
param['nd'] = 2
param['dof'] = 2*param.get('nd') # degrees of freedom per agent
param['n'] = param.get('dof')*param.get('ni')
param['m'] = param.get('dof')*param.get('ni')

# network
param['R_comm'] = 1
param['lambda_a'] = 1

# dynamics
k = 20
param['k1'] = k
param['k2'] = 0.3*k
param['k_c'] = 1e-18
param['R_safe'] = 0.001

# task assignment
param['ta_on'] = True

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