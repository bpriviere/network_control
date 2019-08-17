
from param import param
import numpy as np 
from scipy.linalg import block_diag as block_diag


def augment(x):
	return np.kron( x, np.eye(param.get('nd')))


def to_vec(x):
	return np.array(np.reshape(x, (len(x), 1)))


def get_p():
	# get matrix to extract stacked positions of all agents
	pi = np.zeros( (param.get('ni'), 2*param.get('ni') ))
	for i in range(param.get('ni')):
		pi[i,2*i] = 1
	return np.kron( pi, np.eye(param.get('nd')))


def get_p_a():
	# get matrix to extract stacked positions of all agents
	pi = np.zeros( (param.get('na'), 2*param.get('ni') ))
	for i in range(param.get('na')):
		pi[i,2*i] = 1
	return np.kron( pi, np.eye(param.get('nd')))


def get_v_i(i):
	# get row vector to extract velocity of a single agent
	pi = np.zeros( (1, 2*param.get('ni')))
	pi[0, 2*i + 1] = 1
	return np.kron( pi, np.eye(param.get('nd')))


def get_p_i(i):
	# get row vector to extract position of a single agent
	pi = np.zeros( (1, 2*param.get('ni')))
	pi[0, 2*i] = 1
	return np.kron( pi, np.eye(param.get('nd')))


def get_p_i_idx(i):
	return i*param.get('dof') + np.arange( 0, param.get('nd'))


def get_v_i_idx(i):
	return i*param.get('dof') + param.get('nd') + np.arange( 0, param.get('nd'))


def get_pv_i_idx(i):
	return i*param.get('dof') + np.arange(0, param.get('dof'))


def permute_to_pv():
	pi = np.zeros((2*param.get('ni'),2*param.get('ni')))
	for i in range(param.get('ni')):
		p1_idx = 2*i
		v1_idx = p1_idx + 1
		p2_idx = i
		v2_idx = i + param.get('ni')
		pi[p2_idx, p1_idx] = 1 
		pi[v2_idx, v1_idx] = 1 

	return np.kron( pi, np.eye(param.get('nd')))


def permute_to_pivi():
	pi = permute_to_pv()
	return pi.T


def get_R(phase):
	# rotation matrix, 2D
	return np.vstack(( 
		np.hstack(( np.cos(phase), -np.sin(phase) )),
		np.hstack(( np.sin(phase),  np.cos(phase) )) ))


def get_xd():

	if param.get('case_xd') is "n_node_regular_polygon":
		param['phase_d'] = np.linspace( 2*np.pi/param.get('ni'), 2*np.pi, param.get('ni'))
		param['radius_d'] = np.ones((param.get('ni'),1))

	elif param.get('case_xd') is "3_node_triangle":

		if param.get('ni') is not 3:
			print('Error: Number of Agents Mismatch')
			exit()

		param['phase_d'] = [
			np.pi*2/3,
			np.pi*4/3,
			np.pi*6/3,
		]

		param['radius_d'] = [
			1,
			1,
			1
		]
	print('Desired Phase: ' + str( param.get('phase_d')))


def get_x0():
	min_dist = -np.inf
	
	for _ in range(100):		
		x0 = np.zeros(param.get('n'))
		for i in range(param.get('na')):
			p_i = param.get('plim')*np.random.rand(param.get('nd')) \
				- param.get('plim')/2.
			v_i = param.get('vlim')*np.random.rand(param.get('nd')) \
				- param.get('vlim')/2.
			x0[get_p_i_idx(i)] = p_i
			x0[get_v_i_idx(i)] = v_i
		min_dist = get_min_dist( x0)

		if min_dist < param.get('min_dist'):
			# x0[get_pv_i_idx(param.get('na'))] = get_xb0()
			# x0[get_pv_i_idx(param.get('na') + 1)] = get_xl0()
			param['x0'] = to_vec(x0)
			return

	print('Error: Incompatible Initial Conditions')
	exit()


def get_xb0():
	x = np.zeros((param.get('dof')))
	x[get_p_i_idx(0)] = np.asarray([0,0])
	x[get_v_i_idx(0)] = np.zeros((param.get('nd')))
	return x


def get_xl0():
	x = np.zeros((param.get('dof')))
	return x	


def get_min_dist( x):
	min_dist = np.inf 
	for i in range(param.get('ni')):
		pose_i = np.dot(get_p_i(i), x)
		for j in range(param.get('ni')):
			if i is not j:
				pose_j = np.dot(get_p_i(j), x)
				dist = np.linalg.norm( pose_i - pose_j)
				if dist < min_dist:
					min_dist = dist
	return min_dist