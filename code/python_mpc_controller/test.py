from param import param
import numpy as np 
import dynamics
import controller
import plotter
import utilities as util
import subprocess, os, timeit 


def fn2():

	pdf_path = os.path.join( os.getcwd(), param.get('fn_plots'))

	# remove if exists 
	if os.path.exists(pdf_path):
		os.remove(pdf_path)

	np.random.seed(88)

	util.get_x0()
	util.get_xd()

	# trjaectory to linearize about
	fX = []
	fU = []
	fV = []
	# scp trajectory 
	tX = []
	tU = []
	tV = []
	# integrated trajectory with scp control
	scpV = []
	# linearize about trajectory
	bV = []

	# feedback linearizing baseline 
	x_curr = param.get('x0')
	for k,t in enumerate(param.get('T')):
		u_curr = controller.get_fdbk_controller( x_curr,t)
		x_next = x_curr + param.get('dt')*dynamics.get_dxdt( x_curr, u_curr, t)
		v = dynamics.get_V( x_curr, t)

		fX.append(x_curr)
		fU.append(u_curr)
		fV.append(v)
		x_curr = x_next
	fX = np.squeeze(np.asarray(fX))
	fU = np.asarray(fU)

	# scp
	scpX,scpU,bX,bU = controller.get_scp_clf_controller()
	
	# integrate
	tx_curr = param.get('x0')
	x_curr  = param.get('x0')
	for k,t in enumerate( param.get('T')):
		
		ub = util.list_to_vec(bU[ k,:])
		xb = util.list_to_vec(bX[ k,:])

		u_curr = np.transpose(util.list_to_vec(scpU[ k,:]))

		F_k, B_k, d_k = dynamics.get_linear_dynamics( xb, 
			ub, t)

		R_k, w_k = dynamics.get_linear_lyapunov( xb, 
			ub, t)

		tx_next = np.matmul( F_k, tx_curr) + np.matmul( B_k, u_curr) + d_k
		x_next = x_curr + param.get('dt')*dynamics.get_dxdt( x_curr, u_curr, t)

		tv = np.matmul( R_k, tx_curr) + w_k
		scpv = dynamics.get_V( x_curr, t)

		tX.append(tx_curr)
		tV.append(tv)
		scpV.append(scpv)
		bV.append( dynamics.get_V( xb,t))
		
		tx_curr = tx_next
		x_curr = x_next

	tX = np.squeeze(np.asarray(tX))
	bV = np.asarray(bV)
	
	plotter.plot_SS(fX, param.get('T'), title = 'Fdbk Linearize SS')
	plotter.plot_V(fV, param.get('T'), title = 'Fdbk Linearize V')
	plotter.plot_U(fU, param.get('T'), title = 'Fdbk Linearize U')

	plotter.plot_SS(bX, param.get('T'), title = 'What SCP is linearizing about: SS')
	plotter.plot_V(bV, param.get('T'), title = 'What SCP is linearizing about: V')
	# plotter.plot_U(bU, param.get('T'), title = 'What SCP is linearizing about: U')

	plotter.plot_SS(tX, param.get('T'), title = 'What SCP thinks its doing: SS')
	plotter.plot_V(tV, param.get('T'), title = 'What SCP thinks its doing: V')
	# plotter.plot_U(tU, param.get('T'), title = 'What SCP thinks its doing: U')

	plotter.plot_SS(scpX, param.get('T'), title = 'What SCP is actually doing: SS')
	plotter.plot_V(scpV, param.get('T'), title = 'What SCP is actually doing: V')
	# plotter.plot_U(scpU, param.get('T'), title = 'What SCP is actually doing: U')

	plotter.save_figs()
	subprocess.call(["xdg-open", pdf_path])

def fn1():


	pdf_path = os.path.join( os.getcwd(), param.get('fn_plots'))

	# remove if exists 
	if os.path.exists(pdf_path):
		os.remove(pdf_path)

	np.random.seed(88)

	util.get_x0()
	util.get_xd()

	# baseline trajectory
	bX = []
	bU = np.zeros(( len(param.get('T')),param.get('m') ))

	# true perturbed trajectory
	X = []
	Xdot = []
	V = []
	Vdot = []	

	# linearized perturbed trajectory
	tX = []
	tXdot = []
	tV = []
	tVdot = []

	# collect baseline 
	x_curr = param.get('x0')
	for k,t in enumerate(param.get('T')):
		u_curr = util.list_to_vec(bU[ k,:])
		x_next = x_curr + param.get('dt')*dynamics.get_dxdt( x_curr, u_curr, t)

		bX.append(x_curr)
		x_curr = x_next
	bX = np.squeeze(np.asarray( bX))


	# now run two trajectories and look at divergence
	eps_x0 = 0.0*np.random.uniform( size = ( param.get('n') , 1))
	eps_u = 0.*np.random.uniform( size = ( param.get('nt') , param.get('m')))
	x_curr = param.get('x0') + eps_x0
	tx_curr = param.get('x0') + eps_x0
	scpU = bU + eps_u 
	for k,t in enumerate(param.get('T')):

		# current control
		u_curr = np.reshape( np.transpose( scpU[ k,:]), (param.get('m'),1)) 

		# base
		ub = util.list_to_vec(bU[ k,:])
		xb = util.list_to_vec(bX[ k,:])

		# true
		dxdt = dynamics.get_dxdt( x_curr, u_curr, t)
		x_next = x_curr + dxdt * param.get('dt') 
		v = dynamics.get_V( x_curr, t)
		vdot = dynamics.get_LfV( x_curr, t) + np.matmul( dynamics.get_LgV( x_curr, t), u_curr) 

		# approximate 
		F_k, B_k, d_k = dynamics.get_linear_dynamics( xb, \
			ub, t)
		R_k, w_k = dynamics.get_linear_lyapunov( xb, ub, t)

		tx_next = np.matmul( F_k, tx_curr) + \
			np.matmul( B_k, u_curr) + d_k
		tv = np.matmul(R_k, tx_curr) + w_k

		X.append(x_curr)
		Xdot.append(dxdt)
		V.append(v)
		Vdot.append( vdot)

		tX.append(tx_curr)
		tV.append(tv)

		x_curr = x_next
		tx_curr = tx_next

		print('Timestep:' + str(k) + '/' + str(len(param.get('T'))))

	X = np.squeeze(np.asarray(X))
	V = np.squeeze(np.asarray(V))
	
	tX = np.squeeze(np.asarray(tX))
	tV = np.asarray(tV) 
	tXdot = np.gradient(tX, param.get('dt'), axis = 0)
	tVdot = np.gradient(tV, param.get('dt'), axis = 0)

	import matplotlib.pyplot as plt 

	fig = plt.figure()

	for i in range(np.shape(X)[1]):
		plt.plot( param.get('T'), np.abs( X[:,i] - tX[:,i]), label = 'index' + str(i))
		plt.legend()
	plt.yscale('log')
	plt.title( 'State Error')

	fig = plt.figure()
	plt.plot( param.get('T'), np.abs( V - tV) , label = 'tV')
	plt.plot( param.get('T'), np.abs( tlV - V), label = 'tlV')
	plt.title('Lyapunov')
	plt.legend()
	plt.yscale('log')


	plotter.plot_SS( bX, param.get('T'), title = 'Unperturbed State Space')
	plotter.plot_SS( X, param.get('T'), title = 'Nonlinear Perturbed State Space')
	plotter.plot_SS( tX, param.get('T'), title = 'Linear Perturbed State Space')		

	plotter.save_figs()
	subprocess.call(["xdg-open", pdf_path])

def fn3():

	pdf_path = os.path.join( os.getcwd(), param.get('fn_plots'))

	# remove if exists 
	if os.path.exists(pdf_path):
		os.remove(pdf_path)

	np.random.seed(88)

	util.get_x0()
	util.get_xd()

	# baseline trajectory
	bX = []
	bU = np.zeros(( len(param.get('T')),param.get('m') ))

	# true perturbed trajectory
	X = []
	Xdot = []
	V = []
	Vdot = []	

	# linearized perturbed trajectory
	tX = []
	tXdot = []
	tV = []
	tVdot = []

	# collect baseline 
	x_curr = param.get('x0')
	for k,t in enumerate(param.get('T')):
		u_curr = util.list_to_vec(bU[ k,:])
		x_next = x_curr + param.get('dt')*dynamics.get_dxdt( x_curr, u_curr, t)

		bX.append(x_curr)
		x_curr = x_next
	bX = np.squeeze(np.asarray( bX))


	# now run two trajectories and look at divergence
	eps_x0 = 0.0*np.random.uniform( size = ( param.get('n') , 1))
	eps_u = 0.0*np.random.uniform( size = ( param.get('nt') , param.get('m')))
	x_curr = param.get('x0') + eps_x0
	tx_curr = param.get('x0') + eps_x0
	scpU = bU + eps_u 
	for k,t in enumerate(param.get('T')):

		# current control
		u_curr = np.reshape( np.transpose( scpU[ k,:]), (param.get('m'),1)) 

		# base
		ub = util.list_to_vec(bU[ k,:])
		xb = util.list_to_vec(bX[ k,:])

		# true
		dxdt = dynamics.get_dxdt( x_curr, u_curr, t)
		x_next = x_curr + dxdt * param.get('dt') 
		v = dynamics.get_V( x_curr, t)
		vdot = dynamics.get_LfV( x_curr, t) + np.matmul( dynamics.get_LgV( x_curr, t), u_curr) 

		# approximate 
		F_k, B_k, d_k = dynamics.get_linear_dynamics( xb, \
			ub, t)
		R_k, w_k = dynamics.get_linear_lyapunov( xb, ub, t)

		tx_next = np.matmul( F_k, tx_curr) + \
			np.matmul( B_k, u_curr) + d_k
		tv = np.matmul(R_k, tx_curr) + w_k

		X.append(x_curr)
		Xdot.append(dxdt)
		V.append(v)
		Vdot.append( vdot)

		tX.append(tx_curr)
		tV.append(tv)

		x_curr = x_next
		tx_curr = tx_next

		print('Timestep:' + str(k) + '/' + str(len(param.get('T'))))

	X = np.squeeze(np.asarray(X))
	V = np.reshape(np.asarray(V), (len(V), 1))
	Xdot = np.squeeze(np.asarray(Xdot))
	
	tX = np.squeeze(np.asarray(tX))
	tV = np.reshape(np.asarray(tV), (len(tV), 1)) 
	tXdot = np.gradient(tX, param.get('dt'), axis = 0)
	tVdot = np.gradient(tV, param.get('dt'), axis = 0)

	print( np.shape(X))
	print( np.shape(V))
	print( np.shape(Xdot))
	print( np.shape(tX))
	print( np.shape(tV))
	print( np.shape(tXdot))
	print( np.shape(tVdot))

	epa1 = np.linalg.norm( X[:,0:2] - tX[:,0:2], axis = 1)
	eva1 = np.linalg.norm( X[:,2:4] - tX[:,2:4], axis = 1)
	epa2 = np.linalg.norm( X[:,4:6] - tX[:,4:6], axis = 1)
	eva2 = np.linalg.norm( X[:,6:8] - tX[:,6:8], axis = 1)
	epb = np.linalg.norm( X[:,8:10] - tX[:,8:10], axis = 1)
	evb = np.linalg.norm( X[:,10:12] - tX[:,10:12], axis = 1)

	eXdot = np.linalg.norm( Xdot - tXdot, axis = 1)
	eV = np.linalg.norm( V - tV, axis = 1)
	eVdot = np.linalg.norm( Vdot - tVdot, axis = 1)

	import matplotlib.pyplot as plt 

	fig, ax = plt.subplots()
	plt.plot( epa1, eXdot, label = 'eXdot')
	plt.plot( epa1, eV, label = 'eV')
	ax.set_xscale('log')
	ax.set_yscale('log')
	ax.set_xlabel('epa')
	plt.legend()

	fig, ax = plt.subplots()
	plt.plot( eva1, eXdot, label = 'eXdot')
	plt.plot( eva1, eV, label = 'eV')
	ax.set_xscale('log')
	ax.set_yscale('log')
	ax.set_xlabel('eva')
	plt.legend()

	fig, ax = plt.subplots()
	plt.plot( epb, eXdot, label = 'eXdot')
	plt.plot( epb, eV, label = 'eV')
	# ax.set_xscale('log')
	ax.set_yscale('log')
	ax.set_xlabel('epb')
	plt.legend()

	fig, ax = plt.subplots()
	plt.plot( evb, eXdot, label = 'eXdot')
	plt.plot( evb, eV, label = 'eV')
	# ax.set_xscale('log')
	ax.set_yscale('log')
	ax.set_xlabel('evb')
	plt.legend()	


	# plt.plot( eX, eVdot, label = 'eVdot')

	plotter.save_figs()
	subprocess.call(["xdg-open", pdf_path])



fn3()




