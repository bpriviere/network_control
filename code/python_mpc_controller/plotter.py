

import numpy as np 
import matplotlib.pyplot as plt 
import utilities as util
from param import param
from matplotlib.animation import FuncAnimation
from matplotlib.backends.backend_pdf import PdfPages
import os 

def plot_V(V,T):
	fig, ax = plt.subplots()
	Vdot = np.gradient(V, param.get('dt'), axis = 0)

	ax.plot(T,np.squeeze(V),label = 'V')
	ax.plot(T,np.squeeze(Vdot),label = 'Vdot')
	plt.title('Lyapunov Convergence')
	plt.legend()
	ax.grid(True)

def plot_U(U,T):
	fig, ax = plt.subplots()

	for i in range(param.get('m')):
		ax.plot(T,U[:,i],label = 'U' + str(i))
	plt.title('Control Input')
	plt.legend()


def plot_test( V, LgV, LfV, U, T):
	# currently testing lyapunov and derivatives
	
	Vdot_n = np.gradient(V, param.get('dt'), axis = 0)
	Vdot_a = []
	for t in range(len(T)):
		Vdot_a.append( LfV[t] + np.matmul( LgV[t], U[t]))
	Vdot_a = np.asarray(Vdot_a)

	fig, ax = plt.subplots()
	plt.plot( T, np.squeeze(V), label = 'V')
	plt.plot( T, np.squeeze(Vdot_n), label = 'Vdot n')
	plt.plot( T, np.squeeze(Vdot_a), label = 'Vdot a')
	plt.legend()
	plt.title('Testing')

def plot_test2( pe, ve, ae, T):

	fig, ax = plt.subplots()
	plt.plot( T, pe[:,0], label = 'pe_x')
	plt.plot( T, pe[:,1], label = 'pe_y')
	plt.plot( T, ve[:,0], label = 've_x')
	plt.plot( T, ve[:,1], label = 've_y')
	plt.plot( T, ae[:,0], label = 'ae_x')
	plt.plot( T, ae[:,1], label = 'ae_y')
	plt.legend()
	plt.title('Errors')



def plot_SS(X,T):

	fig, ax = plt.subplots()
	plt.axis('equal')

	# plot agents
	for i in range(param.get('ni')):
		P_i = util.get_P(X,i)
		if i < param.get('na'):
			color = param.get('FreeAgentColor')
		else:
			color = param.get('ControlAgentColor')
		# plot trajectory
		ax.plot( P_i[:,0], P_i[:,1], color = color)
		ax.scatter( P_i[0,0], P_i[0,1], color = color, 
			marker = param.get('start_marker'))
		ax.scatter( P_i[-1,0], P_i[-1,1], color = color, 
			marker = param.get('stop_marker'))

	# plot desired
	ax.plot( param.get('pd')[:,0], param.get('pd')[:,1], 
		color = param.get('DesiredTrajectoryColor'))

	# plot centroid
	P_bar = util.get_P_bar( X,T)
	ax.plot( P_bar[:,0], P_bar[:,1], color = param.get('CentroidColor'))
	ax.plot( P_bar[0,0], P_bar[0,1], color = param.get('CentroidColor'),
		marker = param.get('start_marker'))
	ax.plot( P_bar[-1,0], P_bar[-1,1], color = param.get('CentroidColor'),
		marker = param.get('stop_marker'))
	
	plt.title('State Space')

def plot_ss(x,k,scat):

	scat_x = []
	scat_y = []
	scat_c = []
	
	# plot agents
	for i in range(param.get('ni')):
		p_i = util.get_p(x,i)
		scat_x.append(p_i[0])
		scat_y.append(p_i[1])
		if i < param.get('na'):
			scat_c.append(param.get('FreeAgentColor'))
		else:
			scat_c.append(param.get('ControlAgentColor'))
	
	# plot centroid
	my_1 = util.get_my_1()
	p_a = util.get_p_a(x)
	p_c = np.matmul( np.transpose( my_1), p_a)
	scat_x.append( p_c[0])
	scat_y.append( p_c[1])
	scat_c.append( param.get('CentroidColor'))

	# plot desired
	xd = param.get('pd')[k,:]
	scat_x.append( xd[0])
	scat_y.append( xd[1])
	scat_c.append( param.get('DesiredTrajectoryColor'))

	return plt.scatter(np.asarray(scat_x), np.asarray(scat_y), 
		color = scat_c)

def make_gif(X,T):

	fig, ax = plt.subplots()
	xmin,xmax,ymin,ymax = util.get_plot_lim(X,T)

	def animate(i):
		plt.cla()
		ax.set_xlim(xmin,xmax)
		ax.set_ylim(ymin,ymax)
		plt.title('State Space at t = ' + str(np.round(T[i],2)))
		return plot_ss(X[i],i,scat),

	scat = plt.scatter([],[])
	step = np.round(len(T)/param.get('nframes_gif'))
	anim = FuncAnimation(fig, animate, np.arange(0,len(T),step,int))
	anim.save( param.get('fn_gif'))


def show():
	plt.show()

def save_figs():
	fn = os.path.join( os.getcwd(), param.get('fn_plots'))

	pp = PdfPages(fn)
	for i in plt.get_fignums():
		pp.savefig(plt.figure(i))
		plt.close(plt.figure(i))
	pp.close()