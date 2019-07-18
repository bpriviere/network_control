import numpy as np 

param = {
	
	# time
	't0' : 0., 
	'tf' : 20,
	'dt' : 0.1,

	# agents
	'na' : 2,
	'nb' : 1,

	# dimension
	'nd' : 2, 

	# dynamics
	'model' : 'reynolds', 
	'kx' : 1, 
	'kv' : 1,
	'R_des' : 1, 
	'gamma': 3, # relative degree ( do not change)
	
	# controller
	'controller' : 'fdbk', # [ 'empty', 'fdbk', 'clf', 'scp-clf']
	'mpc_horizon' : 10, # in timesteps	
	'control_max' : 1, 
	'p_u' : 0, 
	'p_v' : 2000, 
	'k_fdbk' : 100, 

	# desired trajectory of swarm centroid
	# 0: x_d(t) = (1,1) (point)
	# 1: x_d(t) = (t,t) (line at 45 degrees)
	# 2: x_d(t) = circle 
	'case_xd' : 2,

	# initialization
	'plim' : 1,
	'vlim' : 1,
	'min_dist' : 1, 

	# plotting
	'ControlAgentColor' : 'r', 
	'FreeAgentColor' : 'b',
	'DesiredTrajectoryColor' : 'g',
	'CentroidColor': 'k',
	'start_marker' : 's', 
	'stop_marker' : 'x', 
	'plot_buffer': 0.2, 
	'gif_on': 0, 
	'nframes_gif' : 50, 
	'fn_gif' : 'gif.gif', 
	'fn_plots' : 'plots.pdf', 
}

# time vector
T = np.arange( param.get('t0'), param.get('tf'), param.get('dt'))
param['T'] = np.reshape(T, (len(T),-1))
# total number of agents
param['ni'] = param.get('na') + param.get('nb')
# state dimension
param['n'] = 2 * param.get('nd') * param.get('ni')
# input dimension
param['m'] = param.get('nd') * param.get('nb')

param['u_prev'] = np.zeros((param.get('m'),1))
