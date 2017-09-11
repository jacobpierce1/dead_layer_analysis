# -*- coding: utf-8 -*-

# the purpose of this script is to bulid off of fit_5_peaks to run through all the data. this poses a challenge because
# it is difficult to verify that all the fits actually worked. the strategy is to display two images of 16 plots
# (32 total) for each of the 32 pixel rows. alarming chisq values ( > 2 ) are reported, along with instances in which
# a fit fails to converge. all the pf's are written to a file for later processing. 

# the first step is to identify the peaks and use this to make guesses for p0, the parameter guess for the 
# least squares fit. 

# date: 9.06.17
# author: jacob pierce


# config 
DEBUG = 1   # this will make it do 2 plots instead of 32 * 2
SHOW_PLOT = 1
PRINT_PEAK_POSITIONS = 0
PRINT_PF = 0
PRINT_PFERR = 0
PRINT_FIT_STATUS = 0
PRINT_FILE_NAMES = 0
BREAK_AFTER_1_PLOT = 0
SAVE_DATA = 1



# my files 
from jacob_math import *
from jacob_pyplot import *
from jacob_utils import *
from deadlayer_functions import * 
from peakdetect import peakdetect 
# reload(jacob_math)



## includes 
import numpy as np
import matplotlib.pyplot as plt
import time 
import sys
import pickle
import sqlite3
# import json



def printvar( var ):
    print var, '=', repr(eval(var))
    


    

    #mode = 'a'
    #if( extract_and_write_data.first_time ):
    #    mode = 'w'
    #    extract_and_write_data.first_time = 0  
    #    
    #with open( outfile, mode ) as f:
    #    d = {
    #        'pixel_coords' : pixel_coords,
    #        'fit_id' : fit_id,
    #        'pf' : pf.tolist(),
    #        'pferr' : pferr.tolist(),
    #        'p0' : p0,
    #        'fit_bounds' : fit_bounds
    #    }
    #    json.dump(d, f)
        
#extract_and_write_data.first_time = 1
        



    
# this function estimates the remaining time for a for loop in which the 
# output is expected to take roughly the same amount of time per run.
# thiis is the equivaletn of a static variable in c 
# currently implemented to only be used once per script, could be changed by adding bool reset.
def estimate_time_left( current_iteration, total_iterations, start_time, num_updates=10 ):
    
    if( current_iteration * 1.0 / total_iterations >= estimate_time_left.counter * 1.0 / num_updates ): 
        current_time = time.time()
        estimate_time_left.counter += 1
        print "%d/%d complete, %f mins remaining" \
                    % ( estimate_time_left.counter, num_updates, \
                    (current_time - start_time) / 60.0 * (num_updates - estimate_time_left.counter ) / estimate_time_left.counter )
estimate_time_left.counter = 1


	#if( j*1.0 / num_entries >= num_hundredths_finished / 100.0 )
	#{
	#    time( &current_time );
	#    printf( "%d/100 complete, ETA: %f mins\n", num_hundredths_finished, ( difftime(current_time, start_time) / 60.0 ) * (100 - num_hundredths_finished ) / num_hundredths_finished );
	#    ++num_hundredths_finished;
	#}



# this function creates another guess for p0 in the case of a failed fit based on
# observations about what is gonig wrong. the particular attempts we use depend on the peak id.
def next_fit_attempt( p0, p0_attempt, fit_bounds, fit_bounds_attempt, npeaks, attempt ):
    
    p0_attempt[:] = p0[:]
    fit_bounds_attempt[:] = fit_bounds[:]    
    
    if( attempt==0 ):
        pass
    
    # try 2: reduce the left fit bound 
    elif( attempt == 1 ):
        fit_bounds_attempt[0] += 10
    
    # try 1: scale down both amplitudes  
    elif( attempt==2 ):
        for i in range(4,4+npeaks*2,2):
            p0_attempt[i] = 0.4 * p0[i] 
    
    # try 2: move back the mu values by 1.
    elif( attempt==3 ):
        for i in range( 5,5+npeaks*2,2):
            p0_attempt[i] -= 2

    else:
        return 0

    return 1







# apply peak fit, make sure it worked, try again if not, and then plot .
def apply_peak_fit( pixel_coords, ax, x, y, fit_bounds, p0, npeaks, mu_all, muerr_all, reduc_chisq_all, fit_id=0, outfile="" ):
    fitfunc = n_fitfuncs_abstract( fitfunc_n_alpha_peaks, npeaks )
    
    attempt=0
    p0_attempt = []
    fit_bounds_attempt = []
        
    status = 0
    
    while( not status ):
        if not next_fit_attempt( p0, p0_attempt, fit_bounds, fit_bounds_attempt, npeaks, attempt ):
            if PRINT_FIT_STATUS:
                print "WARNING: peak failed to converge."
            return 0
                
        status, reduc_chisq, dof, pf, pferr = jacob_least_squares( x, y, np.sqrt(y), fit_bounds_attempt, p0_attempt, fitfunc )

        if( status ):
            if PRINT_FIT_STATUS:
                print "INFO: success, breaking "
            break
        else:
            if PRINT_FIT_STATUS:
                print "Unsuccessful fit, trying new parameters..."
            attempt += 1

    # write all relevant data to the db
    if outfile:
        pass #extract_and_write_data( outfile, pixel_coords, fit_id, pf, pferr, p0, fit_bounds )


    # where we are after breaking    
    if PRINT_PF:
        print "pf = " + str(pf)
    
    if PRINT_PFERR:
        print "pferr = " + str(pferr)
    
    
    add_fit_to_plot( ax, x, fit_bounds, pf, pferr, fitfunc )
    
    mu_all.extend( [pf[ range( 5, 3+2*npeaks, 2) ] ]  )
    muerr_all.extend( [pferr[ range( 5, 3+2*npeaks, 2)] ] )
    reduc_chisq_all.extend( [reduc_chisq] )
    return 1
    
    



def process_file( ax, infile, pixel_coords, outfile="" ):
    
    coordsstr = "(%d, %d)" % pixel_coords

    # x axis
    x = np.array( range(4096) )
        
        
    # EXTRACT ALL DATA
    try:
        f = open( infile, "rb" )
    except IOError:
        print "ERROR: unable to open file: " + filename
        sys.exit(0)
    
    # these are passed as reference to avoid unnecessary copying
    efront_histo = [0] * 4096
    
    construct_histo_array( f, efront_histo )
    f.close()
    

    
    # detect the positions of the 5 highest peaks. peak_positions is an array of tuples (peakpos, value)
    # https://stackoverflow.com/questions/6910641/how-to-get-indices-of-n-maximum-values-in-a-numpy-array
    # def peakdetect(y_axis, x_axis = None, lookahead = 200, delta=0):
    peak_positions = peakdetect.peakdetect( efront_histo, lookahead=10 )[0]
    ind = np.argpartition( [z[1] for z in peak_positions ], -5 )[-5:]
    our_peaks = [ peak_positions[z][0] for z in sorted(ind) ]
    
    if PRINT_PEAK_POSITIONS:
        print "INFO: attempting to process peaks " + str(our_peaks)
    
    # do a check on peak values: energy differenc for the pairs should be constant. no check necessary on the 
    # largest peak since it always dominates, the only potential problem is really the smallest peak in the 
    # lowest energy pair.
    if( abs(23 - (our_peaks[1] - our_peaks[0] ) ) > 10 ):
        print "WARNING: invalid peak suspected for pair 1. " 
        # return -1
    
    if( abs(23 - (our_peaks[3] - our_peaks[2] ) ) > 10 ):
        print "WARNING: invalid peak suspected for pair 2. " 
        # sys.exit -1
    
    
    # set the x axis according to the peak position 
    plot_bounds = [our_peaks[0] - 60, our_peaks[-1]+100 ]
    # plot the histogram without fit yet 
    # plt.clf()
    # ax = plt.axes()
    # plot_bounds = [2400, 3300 ]
    plot_histo( ax, x, efront_histo, plot_bounds=plot_bounds, logscale=0,
                    title = "", 
                    xlabel = "", #"Energy (uncalibrated)",
                    ylabel = "" ) #"Counts"       )
    
    if SHOW_PLOT:
        plt.show()
    
    mu_all = []  # this shall store the mu values 
    muerr_all = []
    reduc_chisq_all = [] # store all reduced chisq values.
    
    
    # first pair 
    fit_bounds = [ our_peaks[0]-50, our_peaks[1]+13]
    p0 = [ 6.0, 0.97, 20.2, 2.0 ]
    #p0 += [ 14000.0, our_peaks[0]+8.0 ]
    #p0 += [ 30000.0, our_peaks[1]+8.0 ]
    p0 += [ 20000.0, our_peaks[0]+8.0 ]
    p0 += [ 60000.0, our_peaks[1]+8.0 ]
    apply_peak_fit( pixel_coords, ax, x, efront_histo, fit_bounds, p0, 2, mu_all, muerr_all, reduc_chisq_all, 0, outfile )
    
    
    # second pair 
    fit_bounds = [ our_peaks[2]-30, our_peaks[3] + 13 ]
    p0 = [ 6.0, 0.99, 42.0, 1.6 ]
    p0 += [ 50000.0, our_peaks[2]+8.0 ]
    p0 += [ 100000.0, our_peaks[3]+8.0 ]    
    #p0 += [ 2000.0, our_peaks[2]+8.0 ]
    #p0 += [ 4000.0, our_peaks[3]+8.0 ]
    apply_peak_fit( pixel_coords, ax, x, efront_histo, fit_bounds, p0, 2, mu_all, muerr_all, reduc_chisq_all, 1, outfile )
    
    
    
    # third peak 
    fit_bounds = [ our_peaks[4]-80, our_peaks[4]+14 ]
    p0 = [ 4.0, 0.99, 30.0, 1.0 ]
    p0 += [ 200000.0, our_peaks[4]+8 ]
    apply_peak_fit( pixel_coords, ax, x, efront_histo, fit_bounds, p0, 1, mu_all, muerr_all, reduc_chisq_all, 2, outfile )


    
    
    
    
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
#    fitstr += '\n' + reduc_chisq_str
    fitstr = '\n' + reduc_chisq_str
    #ax.text( 0.02, 0.92, fitstr, transform=ax.transAxes, fontsize=12, verticalalignment='top')
    #ax.text( 0.02, 0.93, coordsstr, transform=ax.transAxes, fontsize=12, verticalalignment='top' )
    
    ax.text( 0.02, 1.05, fitstr, transform=ax.transAxes, fontsize=12, verticalalignment='top')
    ax.text( 0.03, 0.80, coordsstr, transform=ax.transAxes, fontsize=12, verticalalignment='top' )
        
    

    
    # reduced_fitfunc = lambda x: fitfunc_n_alpha_peaks( 1, pf, x )
    #reduced_fitfunc = lambda y: alpha_fit(pf[4+2*0],pf[4+2*0+1],pf[0],pf[1],pf[2],pf[3], y)
    reduced_fitfunc = lambda x: alpha_fit(pf[4+2*0],pf[4+2*0+1],pf[0],pf[1],pf[2],pf[3], x)

    return 1




dimx = 4
dimy = 4
files = [ "deadlayerdet3rt", "deadlayerdet3cent" ]

# plt.gcf().clear()
def make_all_plots():
    plt.clf()
    f, axarr = plt.subplots(dimx, dimy, figsize=(20,10))
    start_time = time.time()
    
    for i in range(dimx):
        for j in range(dimy):

            estimate_time_left( i + j*dimy, dimx * dimy, start_time )
            
            coords = ( i+14, j+14 )
            current_file = files[0] + "/" + files[0] + "_%d_%d.bin" % coords
            if PRINT_FILE_NAMES:
                print "INFO: processing file: " + current_file
            current_file = "../data_from_mary/output/" + current_file
            title = "(%d, %d)" % coords

            process_file( axarr[i,j], current_file, coords )
            
            if( BREAK_AFTER_1_PLOT ):
                plt.show()
                return 1
                
    # plt.setp([a.get_xticklabels() for a in axarr[0, :]], visible=False)
    # plt.setp([a.get_yticklabels() for a in axarr[:, 1]], visible=False)

    # saveplot_med_quality( "../images/rotated_source_plots/", "test" )
    
    
    
    
    
def make_one_plot( coords1, coords2 ):
    pixel_coords = ( coords1, coords2 )
        
    plt.clf()
    ax = plt.axes()
    current_file = files[0] + "/" + files[0] + "_%d_%d.bin" % pixel_coords
    if PRINT_FILE_NAMES:
        print "INFO: processing file: " + current_file
    current_file = "../data_from_mary/output/" + current_file
    
    process_file( ax, current_file, pixel_coords, "../extracted_spectra_details/test.txt" )    
    plt.show()
    
    
    
make_one_plot(16,16)
# make_all_plots()
