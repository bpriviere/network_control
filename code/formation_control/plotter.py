
import numpy as np 
import matplotlib.pyplot as plt 
import utilities as util
import dynamics
from matplotlib.animation import FuncAnimation
from matplotlib.backends.backend_pdf import PdfPages
from param import param
import os 

def show():
	plt.show()


def save_figs():
	fn = os.path.join( os.getcwd(), param.get('fn_plots'))

	pp = PdfPages(fn)
	for i in plt.get_fignums():
		pp.savefig(plt.figure(i))
		plt.close(plt.figure(i))
	pp.close()


def plot_SS(X,T, title = None, fig = None, ax = None):

	if fig is None:
		fig, ax = plt.subplots()
		plt.axis('equal')

	if title is not None:
		plt.title(title)		

	# plot agents
	for i in range(param.get('na')):
		P_i = np.dot( \
			util.get_p_i(i), X).transpose()
		
		# plot trajectory
		ax.plot( P_i[:,0], P_i[:,1], color = param.get('AgentColor'))
		ax.scatter( P_i[0,0], P_i[0,1], color = param.get('AgentColor'), 
			marker = param.get('start_marker'))
		ax.scatter( P_i[-1,0], P_i[-1,1], color = param.get('AgentColor'), 
			marker = param.get('stop_marker'))
	

def plot_ss(x,t, title = None, fig = None, ax = None):

	if fig is None:
		fig, ax = plt.subplots()
		plt.axis('equal')

	if title is not None:
		plt.title(title)		

	# plot agents
	for i in range(param.get('na')):
		p_i = np.dot( \
			util.get_p_i(i), x)
		
		# plot trajectory
		ax.scatter( p_i[0], p_i[1], color = param.get('AgentColor'), 
			marker = param.get('stop_marker'))		

	p_i = dynamics.get_xl(x,t)[0:param.get('nd')]
	ax.scatter( p_i[0], p_i[1], color = param.get('LeaderColor'), 
		marker = param.get('stop_marker'))
	p_i = dynamics.get_xb(x,t)[0:param.get('nd')]
	ax.scatter( p_i[0], p_i[1], color = param.get('BasisColor'), 
		marker = param.get('stop_marker'))


def plot_VsTx( delta, T, title = 'VsTx'):

	fig, ax = plt.subplots()
	plt.imshow( np.abs(delta)) 
	plt.colorbar()  


def plot(T,X, title = None, fig = None, ax = None):
	if fig is None or ax is None:
		fig, ax = plt.subplots()
	plt.plot( T,X)
	if title is not None:
		plt.title(title)


def make_fig():
	fig, ax = plt.subplots()
	plt.axis('equal')
	return fig, ax
	