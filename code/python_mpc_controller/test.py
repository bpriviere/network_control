from param import param
import numpy as np 
import dynamics
import controller
import plotter
import utilities as util
import subprocess, os, timeit 


pdf_path = os.path.join( os.getcwd(), param.get('fn_plots'))

# remove if exists 
if os.path.exists(pdf_path):
	os.remove(pdf_path)

np.random.seed(88)

util.get_x0()
util.get_xd()
#param['controller'] = 'scp'
print('Sim...')
controller.get_scp_clf_controller()
plotter.save_figs()
subprocess.call(["xdg-open", pdf_path])

