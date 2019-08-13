

import autograd.numpy as np 
import matplotlib.pyplot as plt 
import utilities as util
import dynamics
from param import param
from matplotlib.animation import FuncAnimation
from matplotlib.backends.backend_pdf import PdfPages
import os 

def plot_cost( C, title = 'Objective Value'):
	fig, ax = plt.subplots()

	prop_cycle = plt.rcParams['axes.prop_cycle']
	colors = prop_cycle.by_key()['color']
	
	for k,u in enumerate(param.get('controllers')):
		ax.axhline( C[k], label = u, c = colors[k])		

	plt.title(title)
	plt.xlabel('Iteration')
	plt.ylabel('Objective')
	plt.legend()

def plot_cost_iterations( C, title = 'Iteration Cost'):
	fig, ax = plt.subplots()
	ax.plot( C, '-s')
	plt.title(title)

def debug_scp_iteration_plot( tx_next, u_next, xbar, ubar, x0, T, i_iter):

	unl = u_next
	x_curr = x0
	
	Xnl = []
	Vnl_nlx = []
	Vnl_lx = []
	tV_nlx = []
	tV_lx = []

	for k,t in enumerate(T):
		x_next = x_curr + dynamics.get_dxdt( x_curr, unl[:,k], t) * param.get('dt')
		R_k, w_k = dynamics.get_linear_lyapunov( xbar[:,k], ubar[:,k], t)
		Vnl_nlx.append( dynamics.get_V( x_curr, t))
		Vnl_lx.append( dynamics.get_V( tx_next[:,k],t))
		tV_nlx.append( np.matmul( R_k, x_curr) + w_k )
		tV_lx.append( np.matmul( R_k, tx_next[:,k]) + w_k)		
		Xnl.append( x_curr)
		x_curr = x_next

	Xnl = np.asarray(Xnl)
	Vnl_nlx = np.asarray(Vnl_nlx)
	Vnl_lx = np.asarray(Vnl_lx)
	tV_nlx = np.asarray(tV_nlx)
	tV_lx = np.asarray(tV_lx)

	plot_scp_iteration_state( Xnl, np.transpose(tx_next,(1,0,2)), \
		np.transpose(xbar,(1,0,2)), T, title = str(param.get('controller')) + ' State' + \
		'\nIteration: ' + str(i_iter) + '\nTime: ' + str(T[0]))

	plot_scp_iteration_lyapunov( np.squeeze(Vnl_nlx), np.squeeze(Vnl_lx), np.squeeze( tV_nlx), \
		np.squeeze( tV_lx), T, title = str(param.get('controller')) + ' Lyapunov' + \
		'\nIteration: ' + str(i_iter) + '\nTime: ' + str(T[0]))

def plot_scp_iteration_lyapunov( Vnl_nlx, Vnl_lx, tV_nlx, tV_lx, T, title = None):
	fig, ax = plt.subplots()

	ax.plot( T, Vnl_nlx, label = 'NLV NLX', color = 'b')
	ax.plot( T, Vnl_lx, label = 'NLV LX', color = 'r')
	ax.plot( T, tV_nlx, label = 'LV NLX', color = 'c')
	ax.plot( T, tV_lx, label = 'LV LX', color = 'g')

	plt.legend()
	plt.title(title)


def plot_scp_iteration_state( Xnl, tX, tX_prev, T, title = None):
	fig, ax = plt.subplots() 
	plt.axis('equal')

	color_na_nl = 'b'
	color_nb_nl = 'g'

	color_na_txt = 'r'
	color_nb_txt = 'c'

	color_na_txtm1 = 'm'
	color_nb_txtm1 = 'y'

	for i in range(param.get('ni')):

		P_i = np.squeeze(np.dot( \
			util.get_p_i(i), Xnl)).transpose()

		if i < param.get('na'):
			color = color_na_nl
		else:
			color = color_nb_nl
		# plot trajectory
		ax.plot( P_i[:,0], P_i[:,1], color = color)
		ax.scatter( P_i[0,0], P_i[0,1], color = color, 
			marker = param.get('start_marker'))
		ax.scatter( P_i[-1,0], P_i[-1,1], color = color, 
			marker = param.get('stop_marker'))
	ax.plot( np.nan, np.nan, color = color_na_nl, label = 'Nonlinear X Free')
	ax.plot( np.nan, np.nan, color = color_nb_nl, label = 'Nonlinear X Control')		

	for i in range(param.get('ni')):
		
		P_i = np.squeeze(np.dot( \
			util.get_p_i(i), tX)).transpose()

		if i < param.get('na'):
			color = color_na_txt
		else:
			color = color_nb_txt
		# plot trajectory
		ax.plot( P_i[:,0], P_i[:,1], color = color)
		ax.scatter( P_i[0,0], P_i[0,1], color = color, 
			marker = param.get('start_marker'))
		ax.scatter( P_i[-1,0], P_i[-1,1], color = color, 
			marker = param.get('stop_marker'))
	ax.plot( np.nan, np.nan, color = color_na_txt, label = 'Linear X Free')
	ax.plot( np.nan, np.nan, color = color_nb_txt, label = 'Linear X Control')		

	for i in range(param.get('ni')):

		P_i = np.squeeze(np.dot( \
			util.get_p_i(i), tX_prev)).transpose()

		if i < param.get('na'):
			color = color_na_txtm1
		else:
			color = color_nb_txtm1
		# plot trajectory
		ax.plot( P_i[:,0], P_i[:,1], color = color)		
		ax.scatter( P_i[0,0], P_i[0,1], color = color, 
			marker = param.get('start_marker'))
		ax.scatter( P_i[-1,0], P_i[-1,1], color = color, 
			marker = param.get('stop_marker'))
	ax.plot( np.nan, np.nan, color = color_na_txtm1, label = 'Prev X Free')
	ax.plot( np.nan, np.nan, color = color_nb_txtm1, label = 'Prev X Control')


	plt.legend()
	plt.title( title)


def plot_V(V,T,title = 'Lyapunov Convergence'):
	fig, ax = plt.subplots()
	Vdot = np.gradient(V, param.get('dt'), axis = 0)

	ax.plot(T,np.squeeze(V),label = 'V')
	ax.plot(T,np.squeeze(Vdot),label = 'Vdot')
	plt.title(title)
	plt.legend()
	ax.grid(True)

def plot_U(U,T, title = 'Control Input'):
	fig, ax = plt.subplots()

	for i in range(param.get('m')):
		ax.plot(T,U[:,i],label = 'U' + str(i))
	plt.title(title)
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

def plot_test2( pe, ve, ae, T, title = 'Errors'):

	fig, ax = plt.subplots()
	plt.plot( T, pe[:,0], label = 'pe_x')
	plt.plot( T, pe[:,1], label = 'pe_y')
	plt.plot( T, ve[:,0], label = 've_x')
	plt.plot( T, ve[:,1], label = 've_y')
	plt.plot( T, ae[:,0], label = 'ae_x')
	plt.plot( T, ae[:,1], label = 'ae_y')
	plt.legend()
	plt.title(title)



def plot_SS(X,T, title = 'State Space'):

	fig, ax = plt.subplots()
	plt.axis('equal')

	# plot agents
	for i in range(param.get('ni')):
		P_i = np.dot( \
			util.get_p_i(i), np.transpose(X)).transpose()
		
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
	
	plt.title(title)

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