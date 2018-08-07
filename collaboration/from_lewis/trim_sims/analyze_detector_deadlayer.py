# analysis of TRIM MC for alphas transmitted through a 50 nm Si detector dead layer
import numpy as np
from matplotlib import pyplot as plt
from scipy.optimize import curve_fit
from scipy import exp
import os

datadirec = "C:/Users/LV/Desktop/uchicago/anl/source_deadlayer/detector_deadlayer_srim/"
savefilenameGd = 'detector_deadlayer_fittingsGd.csv'
savefilenameCm = 'detector_deadlayer_fittingsCm.csv'

def main():
	# energies in keV
	alphaEnergies = np.asarray(list(range(3100,3301)))
	savefileGd = open(savefilenameGd,'w')
	print('alpha energy, mean (keV), mean error, S. D. (keV), S. D. error, norm, norm error',file=savefileGd)

	means = []
	meanserrs = []
	sds = []
	sdserrs = []
	for alphaEnergy in alphaEnergies: 
		# perform the fitting analysis and return the fits and errors
		params, errors = analysis(alphaEnergy)

		means.append(params[0])
		meanserrs.append(errors[0])
		sds.append(np.absolute(params[1]))
		sdserrs.append(errors[1])

		#print('alpha energy: %f\nmean: %f +- %f\nS. D.: %f +- %f\nnorm: %f +- %f\n ' % (alphaEnergy, params[0],errors[0],params[1],errors[1],params[2],errors[2]))
		print('%f, %f, %f, %f, %f, %f, %f' %(alphaEnergy, params[0],errors[0],params[1],errors[1],params[2],errors[2]), file=savefileGd)

	savefileGd.close()
	popt,pcov = curve_fit(linear, alphaEnergies,means,sigma=meanserrs,p0=[1.0,1.0])
	perr = np.sqrt(np.diag(pcov))
	print('\nGd means\nm: %f +- %f\nb: %f +- %f'%(popt[0],perr[0],popt[1],perr[1]))
	plt.errorbar(alphaEnergies,means,fmt='b.',yerr=meanserrs)
	plt.plot(alphaEnergies,linear(alphaEnergies,*popt),'k')
	plt.show()

	avgsd = np.average(np.asarray(sds), weights=1.0/(np.asarray(sdserrs)*np.asarray(sdserrs)))
	avgsderr = np.sqrt(1.0/np.sum(1.0/(np.asarray(sdserrs)*np.asarray(sdserrs))))
	print('\nGd S. D.:\nS. D.: %f +- %f'%(avgsd, avgsderr))

	plt.errorbar(alphaEnergies,sds,fmt='b.',yerr=sdserrs)
	plt.hlines(avgsd,alphaEnergies[0],alphaEnergies[-1],'k')
	plt.hlines(avgsd+avgsderr,alphaEnergies[0],alphaEnergies[-1],'k',linestyles='dashed')
	plt.hlines(avgsd-avgsderr,alphaEnergies[0],alphaEnergies[-1],'k',linestyles='dashed')
	plt.show()
	print('\ndone with Gd!')

	# energies in keV
	alphaEnergies = np.asarray(list(range(5500,5901)))
	savefileCm = open(savefilenameCm,'w')
	print('alpha energy, mean (keV), mean error, S. D. (keV), S. D. error, norm, norm error',file=savefileCm)

	means = []
	meanserrs = []
	sds = []
	sdserrs = []
	for alphaEnergy in alphaEnergies: 
		# perform the fitting analysis and return the fits and errors
		params, errors = analysis(alphaEnergy)

		means.append(params[0])
		meanserrs.append(errors[0])
		sds.append(np.absolute(params[1]))
		sdserrs.append(errors[1])

		#print('alpha energy: %f\nmean: %f +- %f\nS. D.: %f +- %f\nnorm: %f +- %f\n ' % (alphaEnergy, params[0],errors[0],params[1],errors[1],params[2],errors[2]))
		print('%f, %f, %f, %f, %f, %f, %f' %(alphaEnergy, params[0],errors[0],params[1],errors[1],params[2],errors[2]), file=savefileCm)

	savefileCm.close()
	popt,pcov = curve_fit(linear, alphaEnergies,means,sigma=meanserrs,p0=[1,1])
	perr = np.sqrt(np.diag(pcov))
	print('\nCm means\nm: %f +- %f\nb: %f +- %f'%(popt[0],perr[0],popt[1],perr[1]))
	plt.errorbar(alphaEnergies,means,fmt='r.',yerr=meanserrs)
	plt.plot(alphaEnergies,linear(alphaEnergies,*popt),'k')
	plt.show()

	avgsd = np.average(np.asarray(sds), weights=1.0/(np.asarray(sdserrs)*np.asarray(sdserrs)))
	avgsderr = np.sqrt(1.0/np.sum(1.0/(np.asarray(sdserrs)*np.asarray(sdserrs))))
	print('\nCm S. D.:\nS. D.: %f +- %f'%(avgsd, avgsderr))

	plt.errorbar(alphaEnergies,sds,fmt='b.',yerr=sdserrs)
	plt.hlines(avgsd,alphaEnergies[0],alphaEnergies[-1],'k')
	plt.hlines(avgsd+avgsderr,alphaEnergies[0],alphaEnergies[-1],'k',linestyles='dashed')
	plt.hlines(avgsd-avgsderr,alphaEnergies[0],alphaEnergies[-1],'k',linestyles='dashed')
	plt.show()
	print('\ndone with Cm!\n')

#Gd means (for energy in keV)
# y = ms + b
#m: -0.001602 +- 0.000025
#b: 14.329482 +- 0.079438

#Gd S. D. (in keV)
#weighted mean, since it seems roughly constant
#S. D.: 1.915266 +- 0.001106

#Cm means (for energy in keV)
# y = ms + b
#m: -0.000679 +- 0.000008
#b: 10.434548 +- 0.045485

#Cm S. D. (in keV)
#weighted mean, since it seems roughly constant
#S. D.: 1.841217 +- 0.000740

def analysis(alphaEnergy):
	# find number of events in the file (12 header lines)
	numEvents = sum((1 for i in open(datadirec+'detector_deadlayer_%i_transmit.txt' % alphaEnergy, 'rb'))) - 12

	# open file with transmitted ion data
	trimdata = open(datadirec+'detector_deadlayer_%i_transmit.txt' % alphaEnergy)

	# ignore header info
	for i in range(0,12):
		trimdata.readline()

	# import data (only need energy information)
	energies = [0]*numEvents
	for i in range(0,numEvents):
		columns = trimdata.readline().strip().split()

		# take difference between initial energy and transmitted energy
		energies[i] = alphaEnergy*1E3 - float(columns[3])

	# calculate guesses for the gaussian fitting
	guessmean = sum(energies)/len(energies)
	guessvar = np.sqrt(sum((np.asarray(energies)-guessmean)**2)/len(energies))

	# bin the data into 10 eV bins
	ebins = np.linspace(0,20000,2001)
	hist = np.histogram(np.asarray(energies),bins=ebins)

	# define the centers of the bins
	ebinscenters = ebins-5.0
	ebinscenters = np.delete(ebinscenters,0,0)

	# fit the data
	popt,pcov = curve_fit(gaussian,ebinscenters,hist[0],p0=[guessmean,guessvar,1.0])

	# compute the errors in the fit
	perr = np.sqrt(np.diag(pcov))

	trimdata.close()

	#plt.plot(ebinscenters,hist[0],'b.')
	#plt.plot(ebinscenters, gaussian(ebinscenters,*popt),'r')
	#plt.show()

	return (popt*1E-3,perr*1E-3)

# gaussian function
def gaussian(x,u,s,n):
	return n*exp(-((x-u)**2)/(2*s**2))

# linear function
def linear(x,m,b):
	return m * x + b

main()
