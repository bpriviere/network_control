
from param import param
import dynamics
import utilities as util
import numpy as np 
import task_assignment as ta

def get_dynamics_u(x,t):

	# get
	A = dynamics.get_A(x)
	pi_pv = util.permute_to_pv()

	# formation task assignment 
	pi_ta_nd = np.kron(ta.get_centralized_ta(x,t), \
		np.eye(param.get('nd')))
	pi_ta_dof = np.kron(ta.get_centralized_ta(x,t), \
		np.eye(param.get('dof')))
	T = np.dot(np.dot(
		pi_ta_nd.T, dynamics.get_T(x,t)), pi_ta_nd)	
	Tv = np.dot(np.dot(
		pi_ta_dof.T, dynamics.get_Tv(x,t)), pi_ta_dof)
 
	# leader/basis	
	x_l = dynamics.get_xl(x,t)
	x_b = dynamics.get_xb(x,t)
	my_1 = np.kron( np.ones((param.get('ni'),1)),\
		np.eye(param.get('dof')))

	# transform state into global coordinate frame
	z = x - np.dot(my_1,x_l) - \
		np.dot(np.dot(Tv, my_1),x_b)	

	# Laplacian
	L = dynamics.get_L(x)
	L = np.dot(
		np.kron(L, np.eye(param.get('nd'))), T)

	I = np.eye(param.get('ni')*param.get('nd')) 

	U = np.dot(np.dot(
		np.hstack((-param.get('k1')*(L+I), -param.get('k2')*(L+I))),
		pi_pv), z)

	# collision avoidance
	# not implemented 
	# for i in range(param.get('ni')):
	# 	p_i = np.dot(util.get_p_i(i),x)
	# 	for j in range(param.get('ni')):
	# 		p_j = np.dot(util.get_p_i(j),x)
	# 		dist = np.linalg.norm( p_i - p_j)

	# 		idx = i*param.get('nd') + np.arange(0, param.get('nd'))

	# 		U[idx] = U[idx] - \
	# 			0*(param.get('k_c')*(p_j-p_i)/(dist*(dist - param.get('R_safe'))))

	return U
