# this should work right off the bat, just try running it

# i have defined some functions specifically for gd or for cm and these functions
# are not completely symmetric, i. e. the cm ones might return more/different
# things than the gd ones.

# believe that this script should be self-contained, no dependencies on any
# other functions I have written.

# think that you probably don't need all these packages imported, but I didn't
# go through and check very carefully.

# any deadlayer parameter/variable refers to the source deadlayer, not the
# detector deadlayer. The detector deadlayer was taken to be 50 nm of Si,
# and the distribution uses the results from that model.

# also note that these convolutions do not take into account the small energy
# dependence of the losses. It is ignored, and the energy for generating the
# loss distributions is takent to be the mean incoming energy.

# the output plots show the evolution of the distribution as it is convolved
# with the different loss mechanisms.

import numpy as np
import math
from matplotlib import pyplot as plt
from scipy import exp
from scipy.optimize import curve_fit, minimize
from scipy.special import erfc, gamma
from scipy.stats import burr, norm
from scipy.integrate import quad
from scipy.interpolate import splev
import sys
import time

showPlots = True # you need this one
verbose = True # not sure if you need this, might have taken it out of these functions

def main():

	# these return the things you want, namely a convolution and its bins in keV
	# and also a convolution for each peak in eV and their bins.

	# everything is normalized and the cm peaks are already the correct proportion.
	makegdlineshape(1000, 0, 200, 3000, 0, verbose=True)
	print('')
	makecmlineshape(1000, 0, 200, 3000, 0, verbose=True)

########## LINESHAPE FUNCTIONS ##########
# find how much energy to shift the convolution over to match the data
def findgdshift(conv, convbins, gdenergypopt):
	maxconv = convbins[np.argmax(conv)]
	maxdata = np.linspace(3000,3300,301)[np.argmax(modelbortels1(np.linspace(3000,3300,301),*gdenergypopt))]

	return maxdata - maxconv

# make the Gd lineshape
def makegdlineshape(sourcethickness, deadlayerthickness, distributionlengthscale,sdnoise, angle, verbose=False):
	start = time.time()	
	alphaType = 'gd'

	energy = get_decayenergy(alphaType)
	sourcedist, ebins = make_gdsourcedistconv_nonhomog_angle(sourcethickness, deadlayerthickness, distributionlengthscale, angle, verbose=verbose)
	sourcemeanenergy = np.average(ebins, weights=sourcedist)

	print('source mean energy loss %f keV'%((energy-sourcemeanenergy)/1000.0))

	detectordeadlayerdist = make_detectordeadlayerdist(alphaType)
	conv1 = np.convolve(sourcedist, detectordeadlayerdist, mode='same')
	meanenergyloss = (-0.001602 * sourcemeanenergy + 14.329482 * 1000)
	newmeanenergy = sourcemeanenergy -  meanenergyloss# subtract off mean of detector dead layer loss

	print('detector deadlayer mean energy loss %f keV'%((meanenergyloss)/1000.0))

	# find the proper ebins shift (up to rounding)
	wheremeanshouldbe = newmeanenergy
	wheremeanis = np.average(ebins, weights=conv1)
	shift = int(round(wheremeanis - wheremeanshouldbe))
	conv1 = np.delete(conv1, np.s_[0:shift], 0)
	conv1 = np.append(conv1, np.zeros(shift))

	convmeanenergy = np.average(ebins, weights=conv1)

	# print(newmeanenergy)
	# print(convmeanenergy)

	# plt.figure()
	# plt.plot(ebins, sourcedist, 'k', label='source')
	# plt.plot(ebins, conv1, 'r', label='detector deadlayer')
	# plt.legend()
	# plt.title('sourcedist s:%0.2f d:%0.3f'%(sourcethickness, deadlayerthickness))
	# plt.show()

	detectornoniondist, meanenergyloss = make_detectornoniondist(alphaType, newmeanenergy)
	conv2 = np.convolve(conv1, detectornoniondist, mode='same')
	newmeanenergy = newmeanenergy - meanenergyloss

	print('non-ionizing mean energy loss %f keV'%((meanenergyloss)/1000.0))

	# find the proper ebins shift (up to rounding)
	wheremeanshouldbe = newmeanenergy
	wheremeanis = np.average(ebins, weights=conv2)
	shift = int(round(wheremeanis - wheremeanshouldbe))
	conv2 = np.delete(conv2, np.s_[0:shift], 0)
	conv2 = np.append(conv2, np.zeros(shift))

	convmeanenergy = np.average(ebins, weights=conv2)

	# print(newmeanenergy)
	# print(convmeanenergy)

	# plt.figure()
	# plt.plot(ebins, sourcedist, 'k', label='source')
	# plt.plot(ebins, conv1, 'r', label='detector deadlayer')
	# plt.plot(ebins, conv2, 'b', label='non-ionizing')
	# plt.legend()
	# plt.title('sourcedist s:%0.2f d:%0.3f'%(sourcethickness, deadlayerthickness))
	# plt.show()

	fanodist = make_fanoandnoise(newmeanenergy, sdnoise)
	conv3 = np.convolve(conv2, fanodist, mode='same')

	# find the proper ebins shift (up to rounding)
	wheremeanshouldbe = newmeanenergy
	wheremeanis = np.average(ebins, weights=conv3)
	shift = int(round(wheremeanis - wheremeanshouldbe))
	conv3 = np.delete(conv3, np.s_[0:shift], 0)
	conv3 = np.append(conv3, np.zeros(shift))

	convmeanenergy = np.average(ebins, weights=conv3)
	# print(newmeanenergy)
	# print(convmeanenergy)

	# plt.figure()
	# plt.plot(ebins, sourcedist, 'k', label='source')
	# plt.plot(ebins, conv1, 'r', label='detector deadlayer')
	# plt.plot(ebins, conv2, 'b', label='non-ionizing')
	# plt.plot(ebins, conv3, 'g', label='final distribution')
	# #plt.plot(kevbins*1000.0, convkev, 'g.', label='final distribution (kev)')
	# plt.legend()
	# #plt.xlim([3100000,3200000])
	# plt.title('s:%0.2f d:%0.3f'%(sourcethickness, deadlayerthickness))
	# plt.show()

	indices = np.where(ebins%1000==0)
	kevbins = [ebins[i] for i in indices][0]/1000.0
	convkev = [conv3[i] for i in indices[0]]
	#print(np.sum(convkev)*1000.0)

	# cut off excess points
	cutbegin = np.where(kevbins==3000)[0][0]
	cutend = np.where(kevbins==3201)[0][0]
	kevbins = np.array(kevbins[cutbegin:cutend])
	convkev = np.array(convkev[cutbegin:cutend])

	end = time.time()
	print('elapsed time: %i:%0.2f'%(int((end-start)/60),(end-start)%60.0))

	if showPlots:
		plt.figure()
		plt.plot(ebins, sourcedist, 'k', label='source')
		plt.plot(ebins, conv1, 'r', label='detector deadlayer')
		plt.plot(ebins, conv2, 'b', label='non-ionizing')
		plt.plot(ebins, conv3, 'g', label='final distribution')
		plt.plot(kevbins*1000.0, convkev, 'g.', label='final distribution (kev)')
		plt.legend()
		plt.xlim([3100000,3200000])
		plt.title('a:%i s:%0.2f d:%0.2f l:%0.2f'%(angle, sourcethickness, deadlayerthickness, distributionlengthscale))
		plt.show()

	# cut off excess points
	# cutbegin = np.where(ebins==3100000)[0][0]
	# cutend = np.where(ebins==3200001)[0][0]
	# ebins = ebins[cutbegin:cutend]
	# conv3 = conv3[cutbegin:cutend]

	# the width of the bins of this convolution is 1 keV, in units of keV.
	return convkev, kevbins, conv3, ebins

# make the Cm lineshape
def makecmlineshape(sourcethickness, deadlayerthickness, distributionlengthscale,sdnoise, angle, verbose=False):
	start = time.time()	

	sourcedistall, ebins, sourcedistpeak1, sourcedistpeak2 = make_cmsourcedistconv_nonhomog_angle(sourcethickness, deadlayerthickness, distributionlengthscale, angle, verbose)
	
	print('\npeak 1')
	alphaType = 'cm1'
	energy = get_decayenergy(alphaType)
	sourcedist = sourcedistpeak1

	sourcemeanenergy = np.average(ebins, weights=sourcedist)

	print('source mean energy loss %f keV'%((energy-sourcemeanenergy)/1000.0))

	detectordeadlayerdist = make_detectordeadlayerdist(alphaType)
	conv1 = np.convolve(sourcedist, detectordeadlayerdist, mode='same')
	meanenergyloss = (-0.001602 * sourcemeanenergy + 14.329482 * 1000)
	newmeanenergy = sourcemeanenergy -  meanenergyloss# subtract off mean of detector dead layer loss

	print('detector deadlayer mean energy loss %f keV'%((meanenergyloss)/1000.0))

	# find the proper ebins shift (up to rounding)
	wheremeanshouldbe = newmeanenergy
	wheremeanis = np.average(ebins, weights=conv1)
	shift = int(round(wheremeanis - wheremeanshouldbe))
	conv1 = np.delete(conv1, np.s_[0:shift], 0)
	conv1 = np.append(conv1, np.zeros(shift))

	convmeanenergy = np.average(ebins, weights=conv1)

	# print(newmeanenergy)
	# print(convmeanenergy)

	# plt.figure()
	# plt.plot(ebins, sourcedist, 'k', label='source')
	# plt.plot(ebins, conv1, 'r', label='detector deadlayer')
	# plt.legend()
	# plt.title('sourcedist s:%0.2f d:%0.3f'%(sourcethickness, deadlayerthickness))
	# plt.show()

	detectornoniondist, meanenergyloss = make_detectornoniondist(alphaType, newmeanenergy)
	conv2 = np.convolve(conv1, detectornoniondist, mode='same')
	newmeanenergy = newmeanenergy - meanenergyloss

	print('non-ionizing mean energy loss %f keV'%((meanenergyloss)/1000.0))

	# find the proper ebins shift (up to rounding)
	wheremeanshouldbe = newmeanenergy
	wheremeanis = np.average(ebins, weights=conv2)
	shift = int(round(wheremeanis - wheremeanshouldbe))
	conv2 = np.delete(conv2, np.s_[0:shift], 0)
	conv2 = np.append(conv2, np.zeros(shift))

	convmeanenergy = np.average(ebins, weights=conv2)

	# print(newmeanenergy)
	# print(convmeanenergy)

	# plt.figure()
	# plt.plot(ebins, sourcedist, 'k', label='source')
	# plt.plot(ebins, conv1, 'r', label='detector deadlayer')
	# plt.plot(ebins, conv2, 'b', label='non-ionizing')
	# plt.legend()
	# plt.title('sourcedist s:%0.2f d:%0.3f'%(sourcethickness, deadlayerthickness))
	# plt.show()

	fanodist = make_fanoandnoise(newmeanenergy, sdnoise)
	conv3 = np.convolve(conv2, fanodist, mode='same')	

	# find the proper ebins shift (up to rounding)
	wheremeanshouldbe = newmeanenergy
	wheremeanis = np.average(ebins, weights=conv3)
	shift = int(round(wheremeanis - wheremeanshouldbe))
	conv3 = np.delete(conv3, np.s_[0:shift], 0)
	conv3 = np.append(conv3, np.zeros(shift))

	convmeanenergy = np.average(ebins, weights=conv3)

	# print(newmeanenergy)
	# print(convmeanenergy)

	peak1conv1 = conv1
	peak1conv2 = conv2
	peak1conv3 = conv3
	peak1mean = convmeanenergy

	print('\npeak 2')
	alphaType = 'cm2'
	energy = get_decayenergy(alphaType)
	sourcedist = sourcedistpeak2

	sourcemeanenergy = np.average(ebins, weights=sourcedist)

	print('source mean energy loss %f keV'%((energy-sourcemeanenergy)/1000.0))

	detectordeadlayerdist = make_detectordeadlayerdist(alphaType)
	conv1 = np.convolve(sourcedist, detectordeadlayerdist, mode='same')
	meanenergyloss = (-0.001602 * sourcemeanenergy + 14.329482 * 1000)
	newmeanenergy = sourcemeanenergy -  meanenergyloss# subtract off mean of detector dead layer loss

	print('detector deadlayer mean energy loss %f keV'%((meanenergyloss)/1000.0))

	# find the proper ebins shift (up to rounding)
	wheremeanshouldbe = newmeanenergy
	wheremeanis = np.average(ebins, weights=conv1)
	shift = int(round(wheremeanis - wheremeanshouldbe))
	conv1 = np.delete(conv1, np.s_[0:shift], 0)
	conv1 = np.append(conv1, np.zeros(shift))

	convmeanenergy = np.average(ebins, weights=conv1)

	# print(newmeanenergy)
	# print(convmeanenergy)

	# plt.figure()
	# plt.plot(ebins, sourcedist, 'k', label='source')
	# plt.plot(ebins, conv1, 'r', label='detector deadlayer')
	# plt.legend()
	# plt.title('sourcedist s:%0.2f d:%0.3f'%(sourcethickness, deadlayerthickness))
	# plt.show()

	detectornoniondist, meanenergyloss = make_detectornoniondist(alphaType, newmeanenergy)
	conv2 = np.convolve(conv1, detectornoniondist, mode='same')
	newmeanenergy = newmeanenergy - meanenergyloss

	print('non-ionizing mean energy loss %f keV'%((meanenergyloss)/1000.0))

	# find the proper ebins shift (up to rounding)
	wheremeanshouldbe = newmeanenergy
	wheremeanis = np.average(ebins, weights=conv2)
	shift = int(round(wheremeanis - wheremeanshouldbe))
	conv2 = np.delete(conv2, np.s_[0:shift], 0)
	conv2 = np.append(conv2, np.zeros(shift))

	convmeanenergy = np.average(ebins, weights=conv2)

	# print(newmeanenergy)
	# print(convmeanenergy)

	# plt.figure()
	# plt.plot(ebins, sourcedist, 'k', label='source')
	# plt.plot(ebins, conv1, 'r', label='detector deadlayer')
	# plt.plot(ebins, conv2, 'b', label='non-ionizing')
	# plt.legend()
	# plt.title('sourcedist s:%0.2f d:%0.3f'%(sourcethickness, deadlayerthickness))
	# plt.show()

	fanodist = make_fanoandnoise(newmeanenergy, sdnoise)
	conv3 = np.convolve(conv2, fanodist, mode='same')	

	# find the proper ebins shift (up to rounding)
	wheremeanshouldbe = newmeanenergy
	wheremeanis = np.average(ebins, weights=conv3)
	shift = int(round(wheremeanis - wheremeanshouldbe))
	conv3 = np.delete(conv3, np.s_[0:shift], 0)
	conv3 = np.append(conv3, np.zeros(shift))

	convmeanenergy = np.average(ebins, weights=conv3)

	# print(newmeanenergy)
	# print(convmeanenergy)

	peak2conv1 = conv1
	peak2conv2 = conv2
	peak2conv3 = conv3
	peak2mean = convmeanenergy

	finalconv = peak1conv3 * 0.2310 + peak2conv3 * 0.7690

	indices = np.where(ebins%1000==0)
	kevbins = [ebins[i] for i in indices][0]/1000.0
	convkev = [finalconv[i] for i in indices[0]]
	# print(np.sum(convkev)*1000.0)

	cutbegin = np.where(kevbins==5650)[0][0]
	cutend = np.where(kevbins==5851)[0][0]
	kevbins = np.array(kevbins[cutbegin:cutend])
	convkev = np.array(convkev[cutbegin:cutend])

	end = time.time()
	print('elapsed time: %i:%0.2f'%(int((end-start)/60),(end-start)%60.0))

	if showPlots:
		plt.figure()
		plt.plot(ebins, sourcedistall, 'k', label='source')
		plt.plot(ebins, peak1conv1 * 0.2310 + peak2conv1 * 0.7690, 'r', label='detector deadlayer')
		plt.plot(ebins, peak1conv2 * 0.2310 + peak2conv2 * 0.7690, 'b', label='non-ionizing')
		plt.plot(ebins, finalconv, 'g', label='final distribution')
		plt.plot(kevbins*1000.0, convkev, 'g.', label='final distribution (kev)')
		plt.legend()
		plt.xlim([5650000,5850000])
		plt.title('s:%0.2f d:%0.2f l:%0.2f'%(sourcethickness, deadlayerthickness, distributionlengthscale))
		plt.show()

	# cutbegin = np.where(ebins==5650000)[0][0]
	# cutend = np.where(ebins==5850001)[0][0]
	# ebins = ebins[cutbegin:cutend]
	# peak1conv3 = peak1conv3[cutbegin:cutend]
	# peak2conv3 = peak2conv3[cutbegin:cutend]

	# the width of the bins of this convolution is 1 keV, in units of keV.
	return convkev, kevbins, peak1conv3 * 0.2310, peak2conv3 * 0.7690, ebins

# get the decay energy in eV
def get_decayenergy(alphaType):
	alphaEnergy = 3182.690 * 1000.0
	if alphaType == 'cm1':
		alphaEnergy = 5762.64 * 1000.0
	elif alphaType == 'cm2':
		alphaEnergy = 5804.77 * 1000.0

	return alphaEnergy

# makes the detector dead layer distribution
def make_detectordeadlayerdist(alphaType):
	# these parameters were in keV, multiply by 1000 to get eV
	if alphaType == 'gd':
		#mean = -0.001602 * energy + 14.329482*1000
		mean = 0
		SD = 1.915266*1000
	elif alphaType == 'cm1' or alphaType == 'cm2':
		#mean = -0.000679 * energy + 10.434548*1000
		mean = 0
		SD = 1.841217*1000

	minbin = int(mean-5*SD)
	maxbin = int(mean+5*SD)
	numbins = maxbin - minbin  + 1
	lossbins = np.linspace(minbin, maxbin, numbins)

	losses = norm.pdf(lossbins, loc=mean, scale=SD)

	return losses

# incEnergy must be in eV
def make_detectornoniondist(alphaType, incEnergy):
	incEnergy = incEnergy / 1000.0 # parameters fitted for an energy in keV

	# distribution was modeled for only energies around the peaks, so need to know
	# which peak to take parameters from
	c = 0
	d = 0
	loc = 0
	scale = 0
	if alphaType == 'gd':
		c = 2.602
		d = 1.289
		loc = 0.042 * incEnergy + 3093.0
		scale = 0.224 * incEnergy + 2678.0
	elif alphaType == 'cm1' or alphaType == 'cm2':
		c = 2.560
		d = 1.450
		loc = 0.146 * incEnergy + 2762.0
		scale = 0.116 * incEnergy + 2959.0

	lossbins = np.linspace(100000, 0, 100001)

	# parameters output energy loss in eV
	losses = burr.pdf(lossbins,c,d,loc,scale)

	mean = float(burr.stats(c, d, loc=loc, scale=scale, moments='m'))	
	# return the loss and some additional information
	return (losses, mean)

# makes the electron-hole pair statistics and the electronics noise distributions
def make_fanoandnoise(incEnergy, sdnoise):
	# in eV
	mean = 0.0
	SD = np.sqrt((2.35*np.sqrt(3.67 * incEnergy * 0.13))**2 + sdnoise**2)

	minbin = int(mean-5*SD)
	maxbin = int(mean+5*SD)
	numbins = maxbin - minbin  + 1
	lossbins = np.linspace(minbin, maxbin, numbins)

	losses = norm.pdf(lossbins, loc=mean, scale=SD)

	return losses

# perform a convolution to generate the source distribution
# but make it non-homogeneous
# and also angle-dependent
def make_gdsourcedistconv_nonhomog_angle(sourcethickness, deadlayerthickness, distributionlengthscale, angle, verbose=False):
	start = time.time()
	ebins = np.linspace(89999,-10000,100000)
	decayenergy = 3182.690 * 1000.0 # in eV
	theta = angle * np.pi / 180.0

	biglosses = ebins * 0.0

	# the thicknesses to create the distributions for
	# exclude what would be zero thickness if no dead layer exists
	# these are 1000 evenly-spaced intervals.
	numsamples = 1000
	thicknesses = np.linspace(sourcethickness + deadlayerthickness, deadlayerthickness, numsamples, endpoint=False)

	weights = []
	for thickness in thicknesses:
		# to account for the rotation of the source
		effectivethickness = thickness/np.cos(theta)

		# Gd parameters
		# mean m:   24.86063
		# SD   a:   35.03233
		# SD   c:    1.5188
		# SD   m:  133.97887
		# SD   b: -183.19972

		mean = 24.677440 * effectivethickness
		SD = linpowsqrt(effectivethickness, *[35.03233, 1.5188, 133.97887, -183.19972])

		# distribute the radioactive particles
		weight = np.exp(-1.0*(thickness-deadlayerthickness)/(distributionlengthscale))
		weights.append(weight)
		
		losses = gennorm(ebins, *[mean, SD])
		biglosses = biglosses + losses * weight

	# normalize to 1
	biglosses = biglosses/np.sum(biglosses)
		
	ebins = np.insert(ebins, 0, np.linspace(189999,90000,100000))
	ebins = np.append(ebins, np.linspace(-10001,-110000,100000))
	ebins = decayenergy - ebins

	# plt.figure()
	# plt.plot(thicknesses, weights)
	# plt.title('non-homog source weighting')
	# plt.xlabel('Thickness (angstrom)')
	# plt.ylabel('weight')
	# plt.show(block=False)

	biglosses = np.insert(biglosses, 0, np.zeros(100000))
	biglosses = np.append(biglosses, np.zeros(100000))

	sourcemean = np.average(ebins, weights=biglosses)
	if verbose:
		print('source mean loss: %0.2f'%((3182690 - sourcemean)/1000.0))
	end = time.time()
	print('elapsed time: %i:%0.2f'%(int((end-start)/60),(end-start)%60.0))
	
	return biglosses, ebins

# perform a convolution to generate the source distribution
# but make it non-homogeneous
# and also angle-dependent
#
# for Cm, make the two peaks independently, then multiply by the correct
# weight and add them together
def make_cmsourcedistconv_nonhomog_angle(sourcethickness, deadlayerthickness, distributionlengthscale, angle, verbose=False):
	start = time.time()
	ebins = np.linspace(89999,-10000,100000)	
	theta = angle * np.pi / 180.0

	# do Cm 1 first
	decayenergy = 5762.64 * 1000.0 # in eV

	biglosses1 = ebins * 0.0

	# the thicknesses to create the distributions for
	# exclude what would be zero thickness if no dead layer exists
	# these are 1000 evenly-spaced intervals.
	numsamples = 1000
	thicknesses = np.linspace(sourcethickness + deadlayerthickness, deadlayerthickness, numsamples, endpoint=False)

	weights = []
	for thickness in thicknesses:
		# to account for the rotation of the source
		effectivethickness = thickness/np.cos(theta)

		# Cm 1 parameters
		# mean m:   16.72004
		# SD   a:   19.31771
		# SD   c:    1.71423
		# SD   m:  137.32133
		# SD   b: -354.10649

		mean = 16.72004 * effectivethickness
		SD = linpowsqrt(effectivethickness, *[19.31771, 1.71423, 137.32133, -354.10649])

		# distribute the radioactive particles
		weight = np.exp(-1.0*(thickness-deadlayerthickness)/(distributionlengthscale))
		weights.append(weight)
		
		losses = gennorm(ebins, *[mean, SD])
		biglosses1 = biglosses1 + losses * weight

	# normalize to 1
	biglosses1 = biglosses1/np.sum(biglosses1)
		
	ebins1 = decayenergy - ebins

	# now do Cm 2
	decayenergy = 5804.77 * 1000.0 # in eV

	biglosses2 = ebins * 0.0

	# the thicknesses to create the distributions for
	# exclude what would be zero thickness if no dead layer exists
	# these are 1000 evenly-spaced intervals.
	numsamples = 1000
	thicknesses = np.linspace(sourcethickness + deadlayerthickness, deadlayerthickness, numsamples, endpoint=False)

	weights = []
	for thickness in thicknesses:
		# to account for the rotation of the source
		effectivethickness = thickness/np.cos(theta)

		# Cm 1 parameters
		# mean m:   16.48520
		# SD   a:   17.80547
		# SD   c:    1.76044
		# SD   m:  136.01682
		# SD   b: -336.08447

		mean = 16.48520 * effectivethickness
		SD = linpowsqrt(effectivethickness, *[17.80547, 1.76044, 136.01682, -336.08447])

		# distribute the radioactive particles
		weight = np.exp(-1.0*(thickness-deadlayerthickness)/(distributionlengthscale))
		weights.append(weight)
		
		losses = gennorm(ebins, *[mean, SD])
		biglosses2 = biglosses2 + losses * weight

	# normalize to 1
	biglosses2 = biglosses2/np.sum(biglosses2)
		
	ebins2 = decayenergy - ebins

	# now we gotta figure out how to add together the distributions
	# ebins aren't quite the same but they overlap somewhere
	# ebins1 should extend to the left of ebins2 because its respective peak
	# is lower
	shift = np.where(ebins1==ebins2[0])
	shift = shift[0][0]

	biglosses2 = np.insert(biglosses2, 0, np.zeros(shift))
	ebins = np.insert(ebins2, 0, ebins1[0:shift])

	# now biglosses2 has more elements than biglosses1, 
	# so we need to figure out how many zeros to add to
	# biglosses1
	toappend = biglosses2.shape[0] - biglosses1.shape[0]
	biglosses1 = np.append(biglosses1, np.zeros(toappend))

	ebinbegin = ebins[0]

	ebins = np.insert(ebins, 0, np.linspace(ebins[0]-100000, ebins[0]-1, 100000))
	ebins = np.append(ebins, np.linspace(ebins[-1]+1, ebins[-1]+100000, 100000))

	biglosses1 = np.insert(biglosses1, 0, np.zeros(100000))
	biglosses1 = np.append(biglosses1, np.zeros(100000))

	biglosses2 = np.insert(biglosses2, 0, np.zeros(100000))
	biglosses2 = np.append(biglosses2, np.zeros(100000))

	# now we need to weight the two peaks and add them together
	biglosses = biglosses1 * 0.2310 + biglosses2 * 0.7690

	# plt.figure()
	# plt.plot(thicknesses, weights)
	# plt.title('non-homog source weighting')
	# plt.xlabel('Thickness (angstrom)')
	# plt.ylabel('weight')
	# plt.show(block=False)

	if verbose:
		sourcemean1 = np.average(ebins, weights=biglosses1)
		sourcemean2 = np.average(ebins, weights=biglosses2)
		sourcemean = np.average(ebins, weights=biglosses)
		print('source mean 1: %0.2f'%(sourcemean1))
		print('source mean 2: %0.2f'%(sourcemean2))
		print('overall source mean: %0.2f'%(sourcemean))
	end = time.time()
	print('elapsed time: %i:%0.2f'%(int((end-start)/60),(end-start)%60.0))
	
	return biglosses, ebins, biglosses1, biglosses2
#########################################

# makes a piecewise function with a linear part and an exponential part
# defined to be exponential below a depth of 100 and linear above 100
def linpowsqrt(x, *params):
	a = params[0]
	c = params[1]
	m = params[2]
	b = params[3]

	lin = m * np.sqrt(x) + b
	power = a * np.sqrt(x)**c

	return lin * np.heaviside(x - 100.0, 0.5) + power * np.heaviside(100.0 - x, 0.5)

# generalized normal distribution that allows for different kurtosis
# it's normalized
# i've fixed the b parameter to 2.4, which has negative excess kurtosis
def gennorm(x, *params):
	#n = params[0]
	u = params[0]
	a = params[1]
	b = 2.4

	return (b)/(2*a*gamma(1.0/b)) * np.exp(-1.0 * ((np.abs(x - u))/(a))**b)*np.heaviside(x,1)

main()
