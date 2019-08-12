
from param import param
from autograd import jacobian
import utilities as util
import autograd.numpy as np


def compute_jacobians():
	param['partial_g_x'] = jacobian(get_g)
	param['partial_f_x'] = jacobian(get_f)

def get_dxdt(x,u,t):

	f = get_f(x)
	g = get_g(x)
	return f + np.dot(g,u)

def get_f( x):

	# drift dynamics

	# free agents
	for i in range(param.get('na')):
		v_i = util.to_vec(np.dot( util.get_v_i(i), x))
		a_i = util.to_vec(reynolds( x, i))

		try:
			f = np.vstack( (f, v_i, a_i))
		except:
			f = np.vstack( (v_i, a_i)) 

	# control agents
	for i in range(param.get('nb')):
		v_i = util.to_vec(np.dot( util.get_v_i(i + param.get('na')), x))
		a_i = util.to_vec(np.zeros( param.get('nd')))
		f = np.vstack( (f, v_i, a_i))
	return f

def get_g( x):
	# control matrix
	g = np.zeros( [param.get('n'), param.get('m')])

	for i in range(param.get('nb')):
		row_idx = 2*param.get('nd')*i + \
			2*param.get('na')*param.get('nd') + param.get('nd')
		col_idx = i*param.get('nd') 
		g[row_idx:row_idx + param.get('nd'), col_idx:col_idx + param.get('nd')] = \
			np.eye(param.get('nd')) 
	return g

def get_dfdx(x):
	dfdx = np.squeeze( param.get('partial_f_x')( x))
	return dfdx

def get_dgdx(x):
	dgdx = np.squeeze( param.get('partial_g_x')( x))
	return dgdx	

def get_linear_dynamics(xbar,ubar,t):

	def prod(dgdx, u):
		s = 0
		for i in range(len(u)):
			s += dgdx[:,i,:]*u[i]
		return np.squeeze(s)

	dfdx = get_dfdx(xbar)
	dgdx = get_dgdx(xbar)
	g = get_g(xbar)
	f = get_f(xbar)

	# continuous time 
	F = dfdx + prod(dgdx,ubar)
	B = g
	d = f + util.to_vec(np.dot( g, ubar)) \
		- util.to_vec(np.dot( dfdx, xbar)) \
		- util.to_vec(np.dot( prod( dgdx, ubar), xbar)) \
		- util.to_vec(np.dot( g, ubar))

	# discrete time
	F_k = np.eye(param.get('n')) + F*param.get('dt') 
	B_k = B*param.get('dt')
	d_k = d*param.get('dt')

	return F_k, B_k, d_k

def reynolds( x, i):
	A   = util.get_A(x)
	p_i = util.to_vec(np.dot( util.get_p_i(i), x))
	v_i = util.to_vec(np.dot( util.get_v_i(i), x))
	a_i = util.to_vec(np.zeros(param.get('nd')))
	for j in range(param.get('ni')):
		if i is not j:
			p_j = util.to_vec(np.dot( util.get_p_i(j), x))
			v_j = util.to_vec(np.dot( util.get_v_i(j), x))

			r_ij = p_j - p_i
			dist = np.linalg.norm(r_ij)

			a_i = a_i + A[i,j]*( \
                param.get('kv')*(v_j - v_i) + \
                param.get('kx')*r_ij*(1 - param.get('R_des')/dist) 
                )
	return a_i 