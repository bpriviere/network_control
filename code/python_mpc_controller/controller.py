

import numpy as np 
import dynamics
import utilities as util 
from param import param
from cvxopt import matrix, solvers

def get_u(x,t):
	# input

	if param.get('controller') is 'empty':
		u = np.zeros( [param.get('m'),1])

	elif param.get('controller') is 'fdbk': 
		u = get_fdbk_controller( x,t)

	elif param.get('controller') is 'clf':
		u = get_clf_controller( x,t)

	elif param.get('controller') is 'mpc-clf':
		u = get_mpc_clf_controller( x,t)

	return u

def get_fdbk_controller( x,t):
	k 			   = np.where(param.get('T') == t)[0][0]
	my_1 		   = util.get_my_1()
	eta 		   = dynamics.get_eta( x,t)
	dvdota_dvb 	   = dynamics.get_dvdota_dvb(x,t)
	xtilde_dot	   = dynamics.get_xtildedot( x,t)
	dvdota_dxtilde = dynamics.get_dvdota_dxtilde( x,t)

	A = np.matmul( np.transpose( my_1), \
		np.matmul( dvdota_dxtilde, xtilde_dot)) - \
		util.list_to_vec( param.get('ad')[k,:])
	B = np.matmul( np.transpose( my_1), dvdota_dvb); 
	K = param.get('k_fdbk')*np.kron( np.eye(param.get('nd')), \
		np.ones((1,param.get('gamma'))));
	u = np.matmul(np.linalg.pinv(B), - A - np.matmul(K,eta));
	return u

def get_clf_controller( x,t):

	LgV = dynamics.get_LgV( x,t) 
	LfV = dynamics.get_LfV( x,t)
	V   = dynamics.get_V( x,t)
	lambda_v = util.get_stabilization_rate()

	# if LfV + lambda_v*V < 0:
	# 	u = np.zeros( (param.get('m'),1))
	# else:
	# 	u = np.matmul( np.linalg.pinv( LgV), -LfV - lambda_v*V)

	# solve using cvx solver
	# min xTQx + pTx 
	# st. G x < h
	#     A x = b
	Q = matrix( np.eye( param.get('m')))
	p = matrix( np.zeros( (param.get('m'),1)))
	G = matrix( LgV)
	h = matrix( -LfV - lambda_v*V)
	A = matrix( np.zeros( (param.get('m'), param.get('m'))))
	b = matrix( np.zeros( (param.get('m'), 1)))

	sol = solvers.qp(Q, p, G, h)
	u = sol['x']
	return util.list_to_vec([x for x in u])

def get_mpc_clf_controller( x,t):

	LgV_t = []
	LfV_t = []
	V_t = []
	x_t 
	for i_t in range(param.get('mpc_horizon')):
		LgV_t.append( dynamics.get_LgV( x,t) )

	LfV = dynamics.get_LfV( x,t)
	V   = dynamics.get_V( x,t)
	lambda_v = util.get_stabilization_rate()

	if LfV + lambda_v*V < 0:
		u = np.zeros( (param.get('m'),1))
	else:
		u = np.matmul( np.linalg.pinv( LgV), -LfV - lambda_v*V)

	# solve using cvx solver
	# min xTQx + pTx 
	# st. G x < h
	#     A x = b
	# Q = matrix( np.eye( param.get('m')))
	# p = matrix( np.zeros( (param.get('m'),1)))
	# G = matrix( LgV)
	# h = matrix( -LfV - lambda_v*V)
	# A = matrix( np.zeros( (param.get('m'), param.get('m'))))
	# b = matrix( np.zeros( (param.get('m'), 1)))

	# sol = solvers.qp(Q, p, G, h)
	# u = sol['x']
	return util.list_to_vec([x for x in u])
