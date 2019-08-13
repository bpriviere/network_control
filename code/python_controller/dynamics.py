

import autograd.numpy as np 
from autograd import grad, elementwise_grad, jacobian
import utilities as util
from param import param


def compute_jacobians():
	param['partial_eta_x'] = jacobian(get_eta, 0)
	param['partial_eta_t'] = jacobian(get_eta, 1)
	param['partial_g_x'] = jacobian(get_g, 0)
	param['partial_f_x'] = jacobian(get_f, 0)

def get_dxdt(x,u,t):
    f = get_f(x)
    g = get_g(x)
    dxdt = f + np.dot(g,u);
    return dxdt

def get_f(x):
	# drift dynamics
	if param.get('model') is 'reynolds':

		# free agents
		for i in range(param.get('na')):
			v_i = np.dot( util.get_v_i(i), x)
			a_i = reynolds( x, i)
			try:
				f = np.vstack( (f, v_i, reynolds( x, i)))
			except:
				f = np.vstack( (v_i, reynolds( x, i))) 

		# control agents
		for i in range(param.get('nb')):
			v_i = np.dot( util.get_v_i(i + param.get('na')), x)
			a_i = np.zeros([param.get('nd'),1])
			f = np.vstack( (f, v_i, a_i))
	return f

def get_g(x):
	# control matrix
	g = np.zeros( [param.get('n'), param.get('m')])

	for i in range(param.get('nb')):
		row_idx = 2*param.get('nd')*i + \
			2*param.get('na')*param.get('nd') + param.get('nd')
		col_idx = i*param.get('nd') 
		g[row_idx:row_idx + param.get('nd'), col_idx:col_idx + param.get('nd')] = \
			np.eye(param.get('nd')) 
	return g

def get_LgV( x,t):
	dVdeta = get_dVdeta(x,t)
	detadx = get_detadx(x,t)
	g 	   = get_g(x)
	return np.dot( dVdeta, np.dot( detadx, g))

def get_LfV( x,t):
	dVdeta = get_dVdeta(x,t)
	detadx = get_detadx(x,t)
	detadt = get_detadt(x,t)
	f 	   = get_f(x)
	return np.dot( dVdeta, np.dot( detadx, f) + detadt)

def get_V( x,t):
	# Lyapunov (not velocity)
	P = util.get_lyapunov_metric()
	eta = get_eta( x,t)
	return np.dot( np.transpose(eta), np.dot(P, eta))

def get_dVdeta(x,t):
	return 2 * np.dot( np.transpose(get_eta(x,t)), \
		util.get_lyapunov_metric())

def get_detadx(x,t):
	detadx = np.squeeze( param.get('partial_eta_x')( x,t))
	return detadx

def get_detadt( x,t):
	detadt = param.get('partial_eta_t')( x,t).reshape( \
		(param.get('gamma')*param.get('nd'),1))
	return detadt

def get_eta( x,t):
	# tracking output dynamics
	my_1 = util.get_my_1()
	p_a = np.dot( util.get_p_a(), x)
	v_a = np.dot( util.get_v_a(), x)
	vdot_a = get_vdot_a( x)

	y   = np.dot( np.transpose(my_1), p_a) \
		- util.get_pd_t(t)
	yp  = np.dot( np.transpose(my_1), v_a) \
		- util.get_vd_t(t)
	ypp = np.dot( np.transpose(my_1), vdot_a) \
		- util.get_ad_t(t)

	eta = np.vstack( (y, yp, ypp))
	eta = np.dot( util.permute_eta_rows(), eta)

	return eta

def get_vdot_a( x):
	
	for i in range(param.get('na')):
		try: 
			vdot_a = np.vstack( (vdot_a, reynolds( x,i)))
		except:
			vdot_a = reynolds( x,i)
	return vdot_a

def reynolds( x, i):
	A   = util.get_A(x)
	p_i = np.dot( util.get_p_i(i), x)
	v_i = np.dot( util.get_v_i(i), x)
	a_i = np.zeros((param.get('nd'),1))
	for j in range(param.get('ni')):
		if i is not j:
			p_j = np.dot( util.get_p_i(j), x)
			v_j = np.dot( util.get_v_i(j), x)

			r_ij = p_j - p_i
			dist = np.linalg.norm(r_ij)

			a_i = a_i + A[i,j]*( \
                param.get('kv')*(v_j - v_i) + \
                param.get('kx')*r_ij*(1 - param.get('R_des')/dist) 
                )
	return a_i 

def get_linear_dynamics( xbar, ubar, t):

	def prod(dgdx, u):
		s = 0
		for i in range(len(u)):
			s += dgdx[:,i,:]*u[i]
		return np.squeeze(s)

	dfdx = get_dfdx(xbar,t)
	dgdx = get_dgdx(xbar)
	g = get_g(xbar)
	f = get_f(xbar)

	# continuous time 
	F = dfdx + prod(dgdx,ubar)
	B = g
	d = f + np.dot( g, ubar) \
		- np.dot( dfdx, xbar) \
		- np.dot( prod( dgdx, ubar), xbar) \
		- np.dot( g, ubar)

	# discrete time
	F_k = np.eye(param.get('n')) + F*param.get('dt') 
	B_k = B*param.get('dt')
	d_k = d*param.get('dt')

	return F_k, B_k, d_k

def get_linear_lyapunov( xbar, ubar, t):

	dVdeta = get_dVdeta( xbar, t)
	detadx = get_detadx( xbar, t)
	V = get_V( xbar, t)

	R_k = np.dot( dVdeta, detadx)
	w_k = V - np.dot( dVdeta, np.dot( detadx, xbar))
	return R_k, w_k

def get_dgdx(x):
	dgdx = np.squeeze( param.get('partial_g_x')( x))
	return dgdx

def get_dfdx(x,t):
	dfdx = np.squeeze( param.get('partial_f_x')( x))
	return dfdx
