
import numpy as np 
import cvxpy as cp
from param import param
# import cvxpy.atoms.affine.sum as cpatom
# import cvxpy.atoms.affine.binary_operators as cpbinary

def get_centralized_ta(x,t):

	# return np.eye(param.get('ni'))

	C = get_cij(x,t)
	Pi = cp.Variable(param.get('ni'),param.get('ni'))
	constraints = []

	for i in range(param.get('ni')):
		constraints.append( 
			sum(Pi[:,i]) == 1)
		constraints.append( 
			sum(Pi[i,:]) == 1)

	obj = cp.Minimize(cp.sum(cp.multiply(C,Pi))) 
	prob = cp.Problem(obj,constraints)
	prob.solve(verbose = False)
	Pi = Pi.value
	# print(Pi)
	# print(Pi.shape)
	return Pi

def get_cij(x,t):
	return np.ones((param.get('ni'),param.get('ni')))