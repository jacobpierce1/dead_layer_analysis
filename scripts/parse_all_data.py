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
SHOW_PLOT = 0
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
import sql_db_manager
# reload(jacob_math)



## includes 
import numpy as np
import matplotlib.pyplot as plt
import time 
import sys
import pickle
import sqlite3
# import json





    
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






# if db connection supplied, then we assume that the entry is supposed to be put in, 
# it makes more sense for the check on whether the data is in the db to be performed 
# before calling this function. when called: 
# apply peak fit, make sure it worked, try again if not, put it in db if successful,
# and then plot. return 1 if successful fit obtained, otherwise 0.

def apply_peak_fit( ax, pixel_coords, fit_id, x, y, fit_bounds, p0, npeaks, mu_all=None, \
        muerr_all=None, reduc_chisq_all=None, sql_conn=None ):
    
    fitfunc = n_fitfuncs_abstract( fitfunc_n_alpha_peaks, npeaks )
    
    attempt=0
    p0_attempt = []
    fit_bounds_attempt = []
        
    status = 0
    
    while( not status ):
        
        if not next_fit_attempt( p0, p0_attempt, fit_bounds, fit_bounds_attempt, \
                npeaks, attempt ):
            
            if PRINT_FIT_STATUS:
                print "WARNING: peak failed to converge."
            return 0
                
        ret = jacob_least_squares( x, y, np.sqrt(y), fit_bounds_attempt, p0_attempt, fitfunc )
        
        if ret is not None:
            if PRINT_FIT_STATUS:
                print "INFO: success, breaking "
            break
        else:
            if PRINT_FIT_STATUS:
                print "Unsuccessful fit, trying new parameters..."
            attempt += 1

    # now we have broken out of the loop after a successful fit. unpack the results.
    reduc_chisq, dof, pf, pferr = ret         
                              
    # write all relevant data to the db
    if sql_conn is not None:
        sql_db_manager.insert_fit_data_into_db( sql_conn, pixel_coords, fit_id, 1, attempt, reduc_chisq, pf, \
                pferr, p0, fit_bounds ) 

    # where we are after breaking    
    if PRINT_PF:
        print "pf = " + str(pf)
    
    if PRINT_PFERR:
        print "pferr = " + str(pferr)
    
    
    add_fit_to_plot( ax, x, fit_bounds, pf, pferr, fitfunc )
    
    mu_all.extend( [pf[ 5:3+2*npeaks:2 ] ]  )
    muerr_all.extend( [pferr[ 5:3+2*npeaks:2 ] ] )
    reduc_chisq_all.append( reduc_chisq )
    
    return 1
    
    


# do all behavior specific to our application: read file, fit peaks, make the plot.
# if sql_conn is supplied, they we assume that it is a valid open sqlite connection
# and check db for whether this file has been processed before.
def process_file( ax, infile, pixel_coords, sql_conn=None ):
    
    coordsstr = "(%d, %d)" % pixel_coords

    # x axis, needed no matter what.
    x = np.array( range(5000) )
        
    # EXTRACT ALL DATA
    try:
        f = open( infile, "rb" )
    except IOError:
        print "ERROR: unable to open file: " + infile
        sys.exit(0)
    
    # these are passed as reference to avoid unnecessary copying
    efront_histo = [0] * len(x)
    
    # creat the histogram
    construct_histo_array( f, efront_histo )
    f.close()
    
    
    
    # detect the positions of the 5 highest peaks. peak_positions is an array of tuples (peakpos, value)
    # https://stackoverflow.com/questions/6910641/how-to-get-indices-of-n-maximum-values-in-a-numpy-array
    # even if we found all the fits, for simplicity we still do the peak detection since it is the easiest
    # way to deteremine a good location for the plot boundaries. 
    
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
    
    
    # number of peaks for each fit, need this info before plotting anything so it is defined here.
    num_peaks = [ 2, 2, 1 ] 
    
    
    # determine which fits to perform, 1 = fit must be attempted. by default if no 
    # conn is supplied we process all the fits.
    fits_to_perform = [ 1, 1, 1 ]
    if sql_conn is not None:
        for i in range(len(fits_to_perform)):
            
            # extract
            result = sql_db_manager.read_data_from_db( sql_conn, pixel_coords, i )
            successful_fit, fit_attempt, reduc_chisq, pf, pferr, p0, fit_bounds, fwhm_data = result

            fits_to_perform[i] = not successful_fit 

            if successful_fit:
                fitfunc = n_fitfuncs_abstract( fitfunc_n_alpha_peaks, num_peaks[i] )
                add_fit_to_plot( ax, x, fit_bounds, pf, pferr, fitfunc )
                reduc_chisq_all.append( reduc_chisq )


    
    # initial fit bounds, these are reduced a bit if the fit fails.
    fit_bounds = [\
                    [ our_peaks[0]-50, our_peaks[1]+13 ],  \
                    [ our_peaks[2]-30, our_peaks[3]+13 ],  \
                    [ our_peaks[4]-80, our_peaks[4]+14 ] \
                ]
    
    # initial guesses which get modified later if the fit fails.
    p0 = [
            [ 6.0, 0.97, 20.2, 2.0 ] + [ 20000.0, our_peaks[0]+8.0 ] + [ 60000.0, our_peaks[1]+8.0 ],
            [ 6.0, 0.99, 42.0, 1.6 ] + [ 50000.0, our_peaks[2]+8.0 ] + [ 100000.0, our_peaks[3]+8.0 ],
            [ 4.0, 0.99, 30.0, 1.0 ] + [ 200000.0, our_peaks[4]+8 ]
        ]
   
    
    # loop through the fits that were not in the db and add them if successful.
    for i in range(len(num_peaks)):
        if fits_to_perform[i]:
            apply_peak_fit( ax, pixel_coords, i, x, efront_histo, fit_bounds[i], p0[i], \
                    num_peaks[i], mu_all, muerr_all, reduc_chisq_all, sql_conn=sql_conn )
        
    #
    ## first pair 
    #fit_bounds = [ our_peaks[0]-50, our_peaks[1]+13]
    #p0 = [ 6.0, 0.97, 20.2, 2.0 ]
    ##p0 += [ 14000.0, our_peaks[0]+8.0 ]
    ##p0 += [ 30000.0, our_peaks[1]+8.0 ]
    #p0 += [ 20000.0, our_peaks[0]+8.0 ]
    #p0 += [ 60000.0, our_peaks[1]+8.0 ]
    #apply_peak_fit( pixel_coords, ax, x, efront_histo, fit_bounds, p0, 2, mu_all, muerr_all, reduc_chisq_all, 0, outfile )
    #
    #
    ## second pair 
    #fit_bounds = [ our_peaks[2]-30, our_peaks[3] + 13 ]
    #p0 = [ 6.0, 0.99, 42.0, 1.6 ]
    #p0 += [ 50000.0, our_peaks[2]+8.0 ]
    #p0 += [ 100000.0, our_peaks[3]+8.0 ]    
    ##p0 += [ 2000.0, our_peaks[2]+8.0 ]
    ##p0 += [ 4000.0, our_peaks[3]+8.0 ]
    #apply_peak_fit( pixel_coords, ax, x, efront_histo, fit_bounds, p0, 2, mu_all, muerr_all, reduc_chisq_all, 1, outfile )
    #
    #
    #
    ## third peak 
    #fit_bounds = [ our_peaks[4]-80, our_peaks[4]+14 ]
    #p0 = [ 4.0, 0.99, 30.0, 1.0 ]
    #p0 += [ 200000.0, our_peaks[4]+8 ]
    #apply_peak_fit( pixel_coords, ax, x, efront_histo, fit_bounds, p0, 1, mu_all, muerr_all, reduc_chisq_all, 2, outfile )


    
    
    
    
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




# plt.gcf().clear()
def parse_all_data():

    dimx = 4
    dimy = 4
    
    totalx = 32
    totaly = 32
    
#    startx = 14
#    starty = 14
    
    files = [ "deadlayerdet3rt", "deadlayerdet3cent" ]
    
    start_time = time.time()
  
    with sqlite3.connect( sql_db_manager.db_filename ) as sql_conn:
            
        # loop through each file type.
        hack = 0
        for fileprefix in files:
                
            # these 2 loops loop over grids of dimx x dimy images. 
            for x in range( totalx / dimx):
                for y in range( 0,totaly / dimy ):
                        
                    print str([x,y])    
                    
                    estimate_time_left( x + y*totaly + hack, 2*totalx * totaly, start_time )

                    plt.clf()
                    plt.close()
                    f, axarr = plt.subplots(dimx, dimy, figsize=(20,10))    
            
                    # these loops create the 4x4 grid.
                    for i in range(dimx):
                        for j in range(dimy):
                            
                            coords = ( i + x * dimx, j + y * dimy )
                            current_file = fileprefix + "/" + files[0] + "_%d_%d.bin" % coords
            
                            if PRINT_FILE_NAMES:
                                print "INFO: processing file: " + current_file
            
                            current_file = "../extracted_ttree_data/" + current_file
                
                            process_file( axarr[i,j], current_file, coords, sql_conn=sql_conn)
                            
                            if( BREAK_AFTER_1_PLOT ):
                                plt.show()
                                return 1
                    
                    saveplot_med_quality( '../current_fit_images/', fileprefix + '_4x4_%d_%d' % (x,y) )
            hack += 32*32
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
    
    
    
    
    
# make_one_plot(16,16)
parse_all_data()
