import matplotlib.pyplot as plt
import pandas as pd
import sys
import numpy as np

def plot(DIR):
	"""
	Function to read the global minumum outputs and plot the result
	
	Parameters
	-------------
	
	filename : 'str'
		Path to the files. Taken automatically.
	
	    
	"""
	#DIR = '/Users/will_rocha_starplan/Downloads/Data_StarPlan_Project/Fitting/WLT0/WL17_5_8/Workspace/Processing/Interp_proc/'
	df = pd.read_csv(DIR+'Interp_proc/Best_comb.csv',sep=',', header=1)
	n_genes = df.shape[1] - 3 # number of genes
		
	data = pd.read_csv(DIR+'Interp_proc/Best_comb.csv',sep=',', usecols=['name'], nrows=n_genes)#pd.read_csv((DIR+'Interp_proc/Best_comb.csv', delimiter=",", low_memory=True, usecols=[2], nrows=n_genes)
	spn = DIR+'Store_interp/'+data+'.dat'
	list = spn.T.values.tolist()[0]
		
	if sys.version_info[0] == 3:
		from ENIIGMA.GA import create3
		create3.create_file3f(list)
	else:
		import create
		create.create_file2f(list)
	
	header = []
	for h in range(n_genes):
		header.append('w'+str(h+1))
		
	data = pd.read_csv(DIR+'Interp_proc/Best_comb.csv', sep=',', low_memory=True, usecols=header, nrows=1)
	cmin = data.T.values.tolist()
		
	t1 = pd.read_csv(DIR+'Interp_proc/output_file_final.txt',sep='\s+', header=None)
	Ysp = pd.read_csv(DIR+'Interp_proc/output_file_final.txt',sep='\s+', header=None, usecols=range(1,t1.shape[1],2))
	
	vv = 'Workspace/Processing/'
	
	NDIR = DIR[:len(DIR) - len(vv)]
	
	l1,f1, ef1= np.loadtxt(NDIR+'/New_tau_GA.txt',dtype=float, usecols=(0,1,2)).T
		
	yff = 0.
	crange = range(n_genes)
	ysprange = range(1,t1.shape[1],2)
	for i,j in zip(crange,ysprange):
		#print i,j
		yff += cmin[i]*Ysp[j]
		np.savetxt('Component_'+str(j)+'.comp', np.transpose([t1[0], cmin[i]*Ysp[j]]))
	
	chi_square = np.sum(((f1 - yff)/(ef1))**2)
	fchi = (1. / (len(l1) - 1. - n_genes)) * chi_square
	rmse2 = np.sqrt(np.mean((f1 - yff) ** 2))
	
	#AIC - Akaike criterion
	p = n_genes
	N = len(l1)
	ff = 2*p +((2*p*(p+1))/(N-p-1))
	AIC = chi_square + ff
	
	print(' ')
	print(' ')
	print('Final score values [Best fit]')
	print('----------------------------------')
	print('Reduced chi-square:', fchi)
	print('Chi-square:', chi_square)
	print('Akaike Information Criterion:', AIC)
	print('Root-mean-square error (RMSE):', rmse2)
	print('----------------------------------')
	
	np.savetxt('Final_plot.txt', np.transpose([t1[0], yff]))
	plt.plot(t1[0], yff, color='limegreen', label='Model')
	plt.legend()
