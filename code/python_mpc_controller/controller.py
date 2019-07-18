

import numpy as np 
import dynamics
import utilities as util 
from param import param
import cvxpy as cp
import plotter

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
	B = np.matmul( np.transpose( my_1), dvdota_dvb)
	K = param.get('k_fdbk')*np.kron( np.eye(param.get('nd')), \
		np.ones((1,param.get('gamma'))))
	u = np.matmul(np.linalg.pinv(B), - A - np.matmul(K,eta))

	u = np.clip( u, -param.get('control_max'), param.get('control_max'))
	return u

def get_clf_controller( x,t):

	LgV = dynamics.get_LgV( x,t) 
	LfV = dynamics.get_LfV( x,t)
	V   = dynamics.get_V( x,t)
	lambda_v = util.get_stabilization_rate()

	# cvxpy
	u = cp.Variable(param.get('m'))
	delta_v = cp.Variable()
	constraints = [ 
		LfV + LgV*u + lambda_v * V <= delta_v,
		cp.abs( u)  <= param.get('control_max'),
		delta_v >= 0
	]
	obj = cp.Minimize( \
		cp.sum_squares(u) + \
		param.get('p_v') * delta_v \
		)
	
	prob = cp.Problem(obj, constraints)
	prob.solve(verbose = True, solver = cp.ECOS)
	
	u = util.list_to_vec([x for x in u.value])
	param['u_prev'] = u
	return u

def get_scp_clf_controller():

	xbar, ubar = get_scp_initial_trajectory()
	lambda_v = util.get_stabilization_rate()*param.get('dt')

	T = param.get('T')
	p_v = param.get('p_v')
	control_max = param.get('control_max')

	u  = cp.Variable( param.get('m'), len(T) )
	tx = cp.Variable( param.get('n'), len(T) )
	tV = cp.Variable( len(T))
	delta_v = cp.Variable( len(T)-1 )

	print('Building Problem')
	constraints = []
	obj_fn = cp.sum_squares( u)
	for k in range( len(T)):
		if k < len(T)-1:
			F_k, B_k, d_k = dynamics.get_linear_dynamics( xbar[k], ubar[k], T[k])

			constraints.append( 
				tV[k+1] <= (1-lambda_v) * tV[k] + delta_v[k])

			constraints.append(
				tx[:,k] == F_k * tx[:,k] + B_k * u[:,k] + d_k)

			constraints.append(
				delta_v[k] >= 0)

			obj_fn += p_v * delta_v[k]

		R_k, w_k = dynamics.get_linear_lyapunov( xbar[k], ubar[k], T[k])
		constraints.append( 
			tV[k] == R_k * tx[:,k] + w_k)
		constraints.append( 
			cp.abs(u) <= control_max)

	obj = cp.Minimize( obj_fn)
	prob = cp.Problem( obj, constraints)

	prob.solve(verbose = True, solver = cp.ECOS)

	X = np.transpose(np.squeeze(np.asarray([y for y in tx.value])))
	U = np.asarray([y for y in u.value])

	plotter.plot_SS( np.squeeze(xbar), T)
	plotter.plot_SS( X, T)

	return


def get_mpc_clf_controller( x,t):
	pass

def get_scp_initial_trajectory():
	# use feedback linearizing solution

	X = []
	U = []
	x_curr = param.get('x0')
	print(param.get('controller'))
	for t in param.get('T'):
		u_curr = get_u( x_curr, t) 
		x_next = x_curr + dynamics.get_dxdt( x_curr, u_curr, t) * param.get('dt')

		X.append(x_curr)
		U.append(u_curr)
		x_curr = x_next

	X = np.asarray(X)
	U = np.asarray(U)
	return X,U
