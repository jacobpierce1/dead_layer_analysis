# -*- coding: utf-8 -*-

# this script uses the functions from 

# date: 8.30.17
# author: jacob pierce



# my files 
from jacob_math import *
from jacob_pyplot import *
from jacob_utils import *
from deadlayer_functions import * 
import sql_db_manager



## includes 
import array
import numpy as np
import matplotlib.pyplot as plt
from math import log10, floor

import sys
import warnings
from scipy.optimize import OptimizeWarning
warnings.simplefilter("error", OptimizeWarning)




    



## constants
test_dir = '../data_from_mary/'
test_file = test_dir + "output/deadlayerdet3rt/deadlayerdet3rt_16_16.bin"
# test_file = test_dir + "output/deadlayerdet3cent/deadlayerdet3cent_1_1.txt"
# test_file = "output/test/test_29_25.bin"







###############################################################
## SCRIPTNG ###################################################
###############################################################

# x axis
x = np.array( range(4096) )

filename = test_file

try:
    f = open( filename, "rb" )
except IOError:
    print "ERROR: unable to open file: " + filename
    sys.exit(0)

# these are passed as reference to avoid unnecessary copying
efront_histo = [0] * 4096
eback_histo  = [0] * 4096

construct_histo_arrays( f, efront_histo, eback_histo )

f.close()


# plot the histogram without fit yet 
plt.clf()
ax = plt.axes()
plot_bounds = [2630, 3300 ]
plot_histo( ax, x, efront_histo, plot_bounds=plot_bounds, logscale=1,
                title = "Rotated Alpha Source, Pixel = (16,16)",
                xlabel = "Energy (keV)",
                ylabel = "Counts"       )



mu_all = []  # this shall store the mu values 
muerr_all = []
reduc_chisq_all = [] # store all reduced chisq values.




# first pair 
fit_bounds = [ 2660, 2760 ]
p0 = [ 6.0, 0.97, 20.2, 2.0 ]
p0 += [ 14000.0, 2730.0 ]
p0 += [ 30000.0, 2755.0 ]
fitfunc = n_fitfuncs_abstract( fitfunc_n_alpha_peaks, 2 )

# apply the laest squares fit.
ret = jacob_least_squares( x, efront_histo, np.sqrt(efront_histo), fit_bounds, p0, fitfunc )

# unpack if successful
if ret is not None:
    reduc_chisq, dof, pf, pferr = ret

add_fit_to_plot( ax, x, fit_bounds, pf, pferr, fitfunc )

mu_all.extend( [pf[5], pf[7]] )
muerr_all.extend( [pferr[5], pferr[7]] )
reduc_chisq_all.extend( [reduc_chisq] )




# second pair 
fit_bounds = [ 2810, 2935 ]
p0 = [ 6.0, 0.99, 42.0, 1.6 ]
p0 += [ 2000.0, 2895.0 ]
p0 += [ 4000.0, 2919.0 ]
fitfunc = n_fitfuncs_abstract( fitfunc_n_alpha_peaks, 2 )


# apply the laest squares fit.
ret = jacob_least_squares( x, efront_histo, np.sqrt(efront_histo), fit_bounds, p0, fitfunc )

# unpack if successful
if ret is not None:
    reduc_chisq, dof, pf, pferr = ret

add_fit_to_plot( ax, x, fit_bounds, pf, pferr, fitfunc )

mu_all.extend( [pf[5], pf[7]] )
muerr_all.extend( [pferr[5], pferr[7]] )
reduc_chisq_all.extend( [reduc_chisq] )




# third peak 
fit_bounds = [ 2990, 3095 ]
p0 = [ 4.0, 0.99, 30.0, 1.0 ]
p0 += [ 200000.0, 3088.0 ]
fitfunc = n_fitfuncs_abstract( fitfunc_n_alpha_peaks, 1 )

# apply the laest squares fit.
ret = jacob_least_squares( x, efront_histo, np.sqrt(efront_histo), fit_bounds, p0, fitfunc )

# unpack if successful
if ret is not None:
    reduc_chisq, dof, pf, pferr = ret

add_fit_to_plot( ax, x, fit_bounds, pf, pferr, fitfunc )

plt.show()

mu_all.extend( [pf[5]] )
muerr_all.extend( [pferr[5]] )
reduc_chisq_all.extend( [reduc_chisq] )





# add a description of the fit to the plot
# note: \mathrm{} is replacement for \text{}
fitstr = "$ f(E) = \\sum_{i,\pm} \\frac{A_i \eta_\pm}{2 \\tau_\pm}    \
                    \\cdot \\exp \\left[   \
                            \\frac{E-\\mu_i}{\\tau_\pm}   \
                            + \\frac{\\sigma^2}{2 \\tau_\pm^2 }  \
                    \\right] \
                    \\cdot \\mathrm{erfc} \\left[   \
                            \\frac{1}{\\sqrt{2}} \\left(  \
                                    \\frac{x-\\mu}{\\sigma}    \
                                    + \\frac{\\sigma}{\\tau_\pm}   \
                            \\right)   \
                    \\right] $"





### FORMATTING FOR THE CASE OF AN ETA=1 fit .
#A_measurements = format_measurement_vector( "A", pf[range(2,pf.size,2)].tolist(), pferr[range(2,pf.size,2)].tolist() )
#
#mu_measurements = format_measurement_vector( "\\mu", pf[range(3,pf.size,2)].tolist(), pferr[range(3,pf.size,2)].tolist() )
#
#tau_str = " $ \\tau = %s \\pm %s $" % tuple( sigfig( pf[0], pferr[0] ) )
#
#sigma_str = "$ \\sigma = %s \\pm %s $" % tuple( sigfig( pf[1], pferr[1] ) )
#
#chisq_str = "$ \\tilde{\\chi}^2 = %s \; (\mathrm{dof} = %d ) $" % ( sigfig(reduc_chisq, 0.01)[0], dof )
#
#fitstr += '\n' + '\n'.join( [ A_measurements, mu_measurements, tau_str, sigma_str, chisq_str ] )




### FORMATTING FOR A FREE ETA FIT
## p format: sigma, eta, tau1, tau2, A1, mu1, ..., A_n, mu_n
#A_measurements = format_measurement_vector( "A", pf[range(4,pf.size,2)].tolist(), pferr[range(4,pf.size,2)].tolist() )
#
#mu_measurements = format_measurement_vector( "\\mu", pf[range(5,pf.size,2)].tolist(), pferr[range(5,pf.size,2)].tolist() )
#
#tau1_str = " $ \\tau_1 = %s \\pm %s $" % tuple( sigfig( pf[2], pferr[2] ) )
#tau2_str = " $ \\tau_2 = %s \\pm %s $" % tuple( sigfig( pf[3], pferr[3] ) )
#
#sigma_str = "$ \\sigma = %s \\pm %s $" % tuple( sigfig( pf[0], pferr[0] ) )
#
#eta_str = "$ \\eta = %s \\pm %s $" % tuple( sigfig( pf[1], pferr[1] ) )
#
#chisq_str = "$ \\tilde{\\chi}^2 = %s \; (\mathrm{dof} = %d ) $" % ( sigfig(reduc_chisq, 0.01)[0], dof )
#
#fitstr += '\n' + '\n'.join( [ A_measurements, mu_measurements, sigma_str, eta_str,  tau1_str, tau2_str, chisq_str ] )



# formatting only with relevant parameters. 
# mu_str= format_measurement_vector( "\\mu", mu_all, muerr_all )
reduc_chisq_str = format_measurement_vector( "\\tilde{\\chi}^2", reduc_chisq_all, 0.01 )

# fitstr += '\n' + '\n'.join( [mu_str, reduc_chisq_str] )
fitstr += '\n' + reduc_chisq_str

ax.text( 0.02, .96, fitstr, transform=ax.transAxes, fontsize=18, verticalalignment='top')
    

# add_legend(ax)




# reduced_fitfunc = lambda x: fitfunc_n_alpha_peaks( 1, pf, x )
#reduced_fitfunc = lambda y: alpha_fit(pf[4+2*0],pf[4+2*0+1],pf[0],pf[1],pf[2],pf[3], y)
reduced_fitfunc = lambda x: alpha_fit(pf[4+2*0],pf[4+2*0+1],pf[0],pf[1],pf[2],pf[3], x)

## x = np.linspace( mu-100, mu+100, 1000 )
## y = reduced_fitfunc(x)
## fwhm, fwhm_x, fwhm_dx = fwhm_from_data( x,y )
#fwhm, fwhm_x, fwhm_dx = fwhm_from_func( reduced_fitfunc, np.linspace(mu_all[0]-30, mu_all[0]+30, 1000) )
#print "fwhm = " + str(fwhm)
#print "fwhm coords = " + str(fwhm_x) 
#print "fwhm dx = " + str(fwhm_dx)

# fwhm, fwhm_x, fwhm_dx = fwhm_from_func( fitfunc, np.linspace(mu-100, mu+100, 1000) )






###############################
# params: sigma, eta, tau1, tau2, A1, u1, ..., An, un )    
# p0 = np.array( [ 100., 10., 100., 10., 10., 10., 10., 10. ])
# p0 = np.array( [2.0]*8 )
# p0 = np.array( [10.0]*8 )
# p0 = np.array( [ 200.0, 2829.0, 200.0, 2853.0, 10.0, 0.5, 1.0, 1.0 ])
# p0 = [10.0, 0.5, 1.0, 1.0, 100.0, 2655.0, 100.0, 2678.0, 200.0, 2831.0, 200.0, 2853.0 ]
# p0 = [100.0, 0.5, 200.0, 200.0, 2000.0, 2829.0, 2000.0, 2853.0 ]
# p0 = [100.0, 0.7, 200.0, 200.0, 2000.0, 2853.0 ]


## parameterts for a fitfunc with n_alpha_peaks_eta1
## fit_bounds = [2710, 2758]
## fit_bounds = [ 2883, 2935 ]
#fit_bounds = [ 2819, 2935 ]
#
#p0 = [6.0, 5.5]
## p0 += [ 1000.0, 2719.0]
## p0 += [ 2000.0, 2744.0]
#p0 += [5000.0, 2896.0 ]
#p0 += [2000.0, 2919.0 ]
#p0 = np.array(p0)
#
## fitfunc = fitfunc_2_alpha_peaks
## fitfunc = fitfunc_n_alpha_peaks_abstract(1)  
## fitfunc = fitfunc_n_alpha_peaks_abstract(4)  
#fitfunc = n_fitfuncs_abstract( fitfunc_n_alpha_peaks_eta1, 2 )




## parameters for fitfunc with a free eta.
#fit_bounds = [ 2655.0, 2757.0 ]
#p0 = [ 4.0, 0.99, 50.0, 10.0]
#p0 += [ 1000.0, 2728.0 ]
#p0 += [ 2000.0, 2760.0 ] 
## p0 += [ 2000.0, 2895.0 ]
## p0 += [ 3000.0, 2919.0 ]
## p0 += [ 3000.0, 3080.0 ]
#fitfunc = n_fitfuncs_abstract( fitfunc_n_alpha_peaks, 2 )

#
#
## winner 
#fit_bounds = [ 2685, 2758 ]
#p0 = [ 4.0, 0.97, 5.2, 5.0 ]
#p0 += [ 2000.0, 2730.0 ]
#p0 += [ 3000.0, 2755.0 ]
#fitfunc = n_fitfuncs_abstract( fitfunc_n_alpha_peaks, 2 )





## winner 
##fit_bounds = [ 2660, 2934 ]
#fit_bounds = [ 2819, 2935 ]
#p0 = [ 6.0, 0.97, 20.2, 2.0 ]
##p0 += [ 14000.0, 2730.0 ]
##p0 += [ 30000.0, 2755.0 ]
#p0 += [ 2000.0, 2909.0 ]
#p0 += [ 3000.0, 2931.0 ]
#fitfunc = n_fitfuncs_abstract( fitfunc_n_alpha_peaks, 2 )



## winner 
#fit_bounds = [ 2660, 2934 ]
#p0 = [ 6.0, 0.97, 20.2, 2.0 ]
#p0 += [ 14000.0, 2730.0 ]
#p0 += [ 30000.0, 2755.0 ]
## p0 += [ 2000.0, 2909.0 ]
## p0 += [ 3000.0, 2931.0 ]
#fitfunc = n_fitfuncs_abstract( fitfunc_n_alpha_peaks, 2 )





# generate the ax and return so that we can keep modifying.
# ax = plot_histo_and_fit( x, efront_histo )



# pf = np.array([round_to_1 (x) for x in pf])
# pferr = np.array([round_to_1 (x) for x in pferr])




# textfit  = "$ f(E) = \\frac{1-\\eta}{\\tau_1} \\cdot \\exp \\left( \\frac{E-\\mu}{\\tau_1} + \\frac{\\sigma^2}{2 \\tau_1^2 }  \\right)    \
#                    \\cdot \\text{erf} \\left(  $" # \\frac{1}{\\sqrt{2}} \\left(  \\frac{x-\\mu}{\\sigma} + \\frac{\\sigma}{\\tau_1} \\right) \\right)  $"
                    #+ (eta/tau2) * np.exp( (x-mu)/tau2 + sigma**2.0/(2*tau2**2.0)) \
                    #* erf( ( (x-mu)/sigma + sigma/tau2) / np.sqrt(2) )  \
                    #  $"
#textfit = ""
#textfit = '$ N(t) = N_0 e^{ -x/x_0} + B $ \n' \
#           '$N_0 = %.0f \pm %.0f$ counts \n' \
#           '$x_0 = %.1f \pm %.1f$ channels \n' \
#           '$B = %.1f\pm %.1f$ counts \n' \
#           '$\\tilde{\chi}^2 = %.2f$ \n' \
#           % (pf[0], pferr[0], pf[1], pferr[1], pf[2], pferr[2], reduc_chisq )
#
# ax.text(0.7, .87, textfit, transform=ax.transAxes, fontsize=18,
#        verticalalignment='top')


