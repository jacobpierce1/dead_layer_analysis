import numdifftools as ndt

# good function for nelder-mead minimize
# calculates the R^2 value for a fit, then multiplies by a big constant
def fitgdpeaks(x, histo):
	# calculate model
	model = modelbortels1(gdchannelbinlefts, *x)

	# calculate residuals and R^2
	ss_res = np.sum((histo - model)**2)
	ss_tot = np.sum((histo - np.mean(histo))**2)
	r2 = 1.0-ss_res/ss_tot

	# scale to help with fitting
	return (1.0 - r2)*100000

# good function for nelder-mead minimize
# calculates the R^2 value for a fit, then multiplies by a big constant
def fitcmpeaks(x, histo):
	# calculate model
	model = modelbortels2(cmchannelbinlefts, *x)

	# calculate residuals and R^2
	ss_res = np.sum((histo - model)**2)
	ss_tot = np.sum((histo - np.mean(histo))**2)
	r2 = 1.0-ss_res/ss_tot

	toreturn =  (1.0 - r2)*100000

	# scale to help with fitting
	return toreturn

# Bortels model for a single peak
def modelbortels1(x, *params):
	a = params[0]
	p1 = params[1]
	n = params[2]
	b = params[3]
	c = params[4]
	s = params[5]

	return a/2.0 * (((1.0 - n)/b)*np.exp((x - p1)/b + (s**2)/(2.0*b**2))*erfc(1.0/np.sqrt(2.0) * ((x - p1)/s + s/b)) + (n/c)*np.exp((x - p1)/c + (s**2)/(2.0*c**2))*erfc(1.0/np.sqrt(2.0) * ((x - p1)/s + s/c)))

# Bortels model for two peaks
def modelbortels2(x, *params):
	a2 = params[0]
	a3 = params[1]
	p2 = params[2]
	p3 = params[3]
	n = params[4]
	b = params[5]
	c = params[6]
	s = params[7]

	return a2/2.0 * (((1.0 - n)/b)*np.exp((x - p2)/b + (s**2)/(2.0*b**2))*erfc(1.0/np.sqrt(2.0) * ((x - p2)/s + s/b)) + (n/c)*np.exp((x - p2)/c + (s**2)/(2.0*c**2))*erfc(1.0/np.sqrt(2.0) * ((x - p2)/s + s/c))) + a3/2.0 * (((1.0 - n)/b)*np.exp((x - p3)/b + (s**2)/(2.0*b**2))*erfc(1.0/np.sqrt(2.0) * ((x - p3)/s + s/b)) + (n/c)*np.exp((x - p3)/c + (s**2)/(2.0*c**2))*erfc(1.0/np.sqrt(2.0) * ((x - p3)/s + s/c)))

# bighisto is a list containing the counts per channel
# i want it normalized so that sum over i of bighisto_i * binwidth_i = 1
# that's why it's normalized and then divided by the binwidth in the args argument

# Cm peak fitting
cm1centerguess = 28400
cm2centerguess = 28575
initguess=[0.23,0.77,cm1centerguess,cm2centerguess,0.1,45,800,30]

# we're going to use this function to do the initial minimization
# args passes arguments to the objective function
nelderans = minimize(fitcmpeaks,initguess,tol=0.00000000001,args=bighisto/np.sum(bighisto)/cmchannelbinwidth, options={'maxiter':15000}, method='Nelder-Mead')
popt = nelderans.x

# this is how to estimate the errors of the minimization
Hfun = ndt.Hessian(fitcmpeaks, full_output=True)
hessian_ndt, info = Hfun(nelderans['x'],histo=bighisto/np.sum(bighisto)/cmchannelbinwidth)
# sometimes this matrix has negative values? i stuck an absolute value in there
perr = np.sqrt(np.abs(np.diag(np.linalg.inv(hessian_ndt))))

# Gd peak fitting
gdcenterguess = 15600
initguess=[1.0,gdcenterguess,0.2,35,500,30]

# we're going to use this function to do the initial minimization
nelderans = minimize(fitgdpeaks,initguess,args=bighisto/np.sum(bighisto)/gdchannelbinwidth,tol=0.000000001, options={'maxiter':15000},method='Nelder-Mead')
popt = nelderans.x

# this is how to estimate the errors of the minimization
Hfun = ndt.Hessian(fitgdpeaks, full_output=True)
hessian_ndt, info = Hfun(nelderans['x'],histo=bighisto/np.sum(bighisto)/gdchannelbinwidth)
perr = np.sqrt(np.abs(np.diag(np.linalg.inv(hessian_ndt))))