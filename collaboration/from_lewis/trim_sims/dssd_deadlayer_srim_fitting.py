# analysis of TRIM MC for alphas transmitted through a 50 nm Si detector dead layer
import numpy as np
from matplotlib import pyplot as plt
from scipy.optimize import curve_fit
from scipy import exp
import os
from scipy.special import erfc, gamma

datadirec = "C:/Users/LV/Desktop/uchicago/anl/source_deadlayer/dssd_deadlayer/"

def main():
	alphaEnergy = 3177.0
	thicknesses = np.array(list(range(22,51)))
	thicknesses = thicknesses * thicknesses
	means = []
	meanserrs = []
	sds = []
	sdserrs = []
	for thickness in thicknesses: 
		# perform the fitting analysis and return the fits and errors
		params, errors = analysis(alphaEnergy, thickness)

		means.append(params[1])
		meanserrs.append(errors[1])
		sds.append(np.absolute(params[2]))
		sdserrs.append(errors[2])

		#print('alpha energy: %f\nmean: %f +- %f\nS. D.: %f +- %f\nnorm: %f +- %f\n ' % (alphaEnergy, params[0],errors[0],params[1],errors[1],params[2],errors[2]))
		#print('%i, %i, %f, %f, %f, %f, %f, %f' %(alphaEnergy, thickness, params[0],errors[0],params[1],errors[1],params[2],errors[2]))
	
	popt,pcov = curve_fit(linear, thicknesses, means,sigma=meanserrs,p0=[1.0,1.0])
	perr = np.sqrt(np.diag(pcov))
	print('\nGd means\nm: %f +- %f\nb: %f +- %f'%(popt[0],perr[0],popt[1],perr[1]))
	plt.errorbar(thicknesses,means,fmt='b.',yerr=meanserrs)
	plt.plot(thicknesses,linear(thicknesses,*popt),'k')
	plt.xlabel('Dead layer thickness (angstrom)')
	plt.ylabel('Mean energy loss (kev)')
	plt.title('Detector dead layer energy loss')
	plt.show()

	popt,pcov = curve_fit(linear, np.sqrt(thicknesses), sds,sigma=sdserrs,p0=[1.0,1.0])
	perr = np.sqrt(np.diag(pcov))
	print('\nGd sds\nm: %f +- %f\nb: %f +- %f'%(popt[0],perr[0],popt[1],perr[1]))
	plt.errorbar(thicknesses,sds,fmt='b.',yerr=sdserrs)
	plt.plot(thicknesses,linear(np.sqrt(thicknesses),*popt),'k')
	plt.xlabel('Dead layer thickness (angstrom)')
	plt.ylabel('Energy loss SD (kev)')
	plt.title('Detector dead layer energy loss')
	plt.show()


def analysis(alphaEnergy, thickness):

	# find number of events in the file (12 header lines)
	numEvents = sum((1 for i in open(datadirec+'%ikev_%iA_Si_transmit.txt' % (alphaEnergy, thickness), 'rb'))) - 12

	# open file with transmitted ion data
	trimdata = open(datadirec+'%ikev_%iA_Si_transmit.txt' % (alphaEnergy, thickness))

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
	guessvar = np.sqrt(np.sqrt(sum((np.asarray(energies)-guessmean)**2)/len(energies)))*10.0

	# bin the data into 100 eV bins
	ebins = np.linspace(0,50000,501)
	hist = np.histogram(np.asarray(energies),bins=ebins)

	# define the centers of the bins
	ebinscenters = ebins-50.0
	ebinscenters = np.delete(ebinscenters,0,0)

	# fit the data
	popt,pcov = curve_fit(gennorm,ebinscenters,hist[0],p0=[200.0,guessmean,guessvar])

	# compute the errors in the fit
	perr = np.sqrt(np.diag(pcov))

	trimdata.close()

	# plt.plot(ebinscenters,hist[0],'b.')
	# plt.plot(ebinscenters, gennorm(ebinscenters,*popt),'r')
	# plt.show()

	return (popt*1E-3,perr*1E-3)

# gaussian function
def gaussian(x,u,s,n):
	return n*exp(-((x-u)**2)/(2*s**2))

# linear function
def linear(x,m,b):
	return m * x + b

# generalized normal distribution that allows for different kurtosis
# it's normalized
# i've fixed the b parameter to 2.4, which has negative excess kurtosis
def gennorm(x, *params):
	n = params[0]
	u = params[1]
	a = params[2]
	b = 2.4

	return n *(b)/(2*a*gamma(1.0/b)) * np.exp(-1.0 * ((np.abs(x - u))/(a))**b)*np.heaviside(x,1)

main()
