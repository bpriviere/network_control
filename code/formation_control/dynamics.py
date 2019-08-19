
from param import param
import utilities as util
import task_assignment as ta 
import numpy as np
from scipy.linalg import block_diag as block_diag


def get_dxdt(x,u,t):
	# full dynamics
	f = get_f(x,t)
	g = get_g(x,t)
	return f + np.dot(g,u)


def get_f(x,t):
	# single integrator dynamics
	pi_pv = util.permute_to_pv()
	I = np.eye(param.get('nd')*param.get('ni'))
	F = np.vstack((
		np.hstack((0.*I,1.*I)),
		np.hstack((0.*I,0.*I))))

	return np.dot(np.dot(np.dot(
		pi_pv.T, F), pi_pv), x)


def get_g(x,t):
	# control matrix
	g = np.zeros((2*param.get('ni'), param.get('ni')))
	for i in range(param.get('ni')):
		g[2*i+1,i] = 1
	return np.kron(g,np.eye(param.get('nd')))


def get_A(x):
	# Adjacency matrix
	
	P = np.dot(util.get_p(),x)	
	X = np.dot( 
		np.ones((param.get('ni'),1)),
		util.to_vec(P[np.mod(np.arange(0,len(P)),2)==False]).T)
	Y = np.dot( 
		np.ones((param.get('ni'),1)),
		util.to_vec(P[np.mod(np.arange(0,len(P)),2)==True]).T)
	D = np.sqrt( 
		np.power( X - X.T + 1e-9,2) + 
		np.power( Y - Y.T + 1e-9,2))
	A = np.exp( -param.get('lambda_a')*D) * \
		1.0 * (D < param.get('R_comm'))
	A = A / np.linalg.norm(A, ord=2, axis=1)
	return A


def get_L(x):
	# Laplacian Matrix
	A = get_A(x)
	D = np.zeros( np.shape(A))
	for i in range( np.shape(A)[0]):
		D[i,i] = sum( A[i,:])
	return D - A


def get_T(x,t):
	# add phase shift for all agents
	T = param.get('radius_d')[0]*util.get_R(param.get('phase_d')[0])
	for i in range(1,param.get('na')): 
		r_i = param.get('radius_d')[i]
		phi_i = param.get('phase_d')[i]
		T = block_diag(T, r_i*util.get_R(phi_i))
	return T

def get_Tv(x,t):
	# add phase shift for all agents
	T = param.get('radius_d')[0]*util.get_R(param.get('phase_d')[0])
	T = block_diag(T, 0*np.eye(param.get('nd')))
	for i in range(1,param.get('na')): 
		r_i = param.get('radius_d')[i] 
		phi_i = param.get('phase_d')[i] 
		T = block_diag(T, r_i*util.get_R(phi_i))
		T = block_diag(T, 0*np.eye(param.get('nd')))
	return T


def get_xl(x,t):
	# leader state
	xl = np.zeros((param.get('dof')))
	return util.to_vec(xl)


def get_xb(x,t):
	# basis state
	xb = np.zeros((param.get('dof')))
	xb[0] = 1
	return util.to_vec(xb)


def get_phase(p_b):
	return np.arctan2(p_b[1],p_b[0])


# def get_f(x,t):

# 	pi_pv = util.permute_to_pv()
# 	pi_ta_nd = np.kron(ta.get_centralized_ta(x,t), \
# 		np.eye(param.get('nd')))	
# 	pi_ta_dof = np.kron(ta.get_centralized_ta(x,t), \
# 		np.eye(param.get('dof')))

# 	T = np.dot(np.dot(
# 		pi_ta_nd.T, get_T(x,t)), pi_ta_nd)

# 	Tv = np.dot(np.dot(
# 		pi_ta_dof.T, get_Tv(x,t)), pi_ta_dof)

# 	L = np.dot(
# 		np.kron( get_L(x), np.eye(param.get('nd'))),
# 		T
# 		)

# 	xl = get_xl(x,t)
# 	xb = get_xb(x,t)
# 	my_1 = np.kron( np.ones((param.get('ni'),1)),\
# 		np.eye(param.get('dof')))

# 	z = x - np.dot(my_1,xl) - \
# 		np.dot(np.dot(Tv, my_1),xb)

# 	I = np.eye(param.get('ni')*param.get('nd'))
# 	F = np.vstack((
# 		np.hstack(( 
# 			0.*I,
# 			1.*I)),
# 		np.hstack(( 
# 			-param.get('k1')*(L+I), 
# 			-param.get('k2')*(L+I))) ))

# 	f = np.dot(np.dot(np.dot( 
# 		pi_pv.T, F), pi_pv), z)

# 	return f	