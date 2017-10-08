# -*- coding: utf-8 -*-

# the purpose of this script is to bulid off of fit_5_peaks to run
# through all the data. this poses a challenge because it is difficult
# to verify that all the fits actually worked. the strategy is to
# display two images of 16 plots (32 total) for each of the 32 pixel
# rows. alarming chisq values ( > 2 ) are reported, along with
# instances in which a fit fails to converge. all the pf's are written
# to a file for later processing.

# the first step is to identify the peaks and use this to make guesses
# for p0, the parameter guess for the least squares fit.

# date: 9.06.17
# author: jacob pierce


# default config. do not adjust here, make a dedicated function that sets 
# config parameters in an abnormal situtation (e.g. debugging). any mistake in
# the db can be very costly in time, either through repairing db which will 
# probably required a specialized function or through rerunnnig the program.


SHOW_HISTO = 0
PRINT_PEAK_POSITIONS = 0
PRINT_PF = 0
PRINT_PFERR = 0
PRINT_FIT_STATUS = 0
PRINT_FILE_NAMES = 0
BREAK_AFTER_1_PLOT = 0
SAVE_DATA = 1
UPDATE_DB = 1
PRINT_FIT_ATTEMPT = 0




# my files 
import libjacob.jmath as jmath
import libjacob.jpyplot as jplt
import libjacob.jutils as jutils


import deadlayer_helpers.functions as dlf


from peakdetect import peakdetect 
import deadlayer_helpers.sql_db_manager as db



## includes 
import numpy as np
import matplotlib.pyplot as plt
import time 
import sys
import sqlite3
import os 




NUM_PEAKS_PER_FEATURE = [ 2, 2, 2 ] 
NUM_PEAKS = np.sum( NUM_PEAKS_PER_FEATURE )
NUM_FEATURES = len( NUM_PEAKS_PER_FEATURE )





# these are the functions that attempt to modify p0_attempt and fit_bounds_attempt.
def default( p0, p0_attempt, fit_bounds, fit_bounds_attempt, npeaks ):
    pass 
    
def increase_left_fit_bound( p0, p0_attempt, fit_bounds, fit_bounds_attempt, npeaks ):
    fit_bounds_attempt[0] += 10
    
def increase_left_fit_bound_more( p0, p0_attempt, fit_bounds, fit_bounds_attempt, npeaks ):
    fit_bounds_attempt[0] += 20
    
def smaller_A( p0, p0_attempt, fit_bounds, fit_bounds_attempt, npeaks ):
    for i in range( 4, 4+npeaks*2, 2 ):
        p0_attempt[i] = 0.4 * p0[i]

def larger_A( p0, p0_attempt, fit_bounds, fit_bounds_attempt, npeaks ):
    # fit_bounds_attempt[1] += 5
    for i in range( 4, 4+npeaks*2, 2 ):
        p0_attempt[i] = 3 * p0[i]

def next_fit( p0, p0_attempt, fit_bounds, fit_bounds_attempt, npeaks ):
    pass




# this function creates another guess for p0 in the case of a failed fit based on
# observations about what is gonig wrong. the particular attempts we use depend on the peak id.
def next_fit_attempt( p0, p0_attempt, fit_bounds, fit_id,
                      fit_bounds_attempt, npeaks,
                      attempt, sql_conn=None, pixel_coords=None ):

    
    # these will be modified and used later as the initial guess.
    p0_attempt[:] = p0[:]
    fit_bounds_attempt[:] = fit_bounds[:]    


        
    # array of functions that will be used to modify p0_attempt and fit_bounds_attempt
    if fit_id == 1:
        modifier_functions = [ default, increase_left_fit_bound, larger_A, next_fit ]

    else:
        modifier_functions = [ default, increase_left_fit_bound, 
                                    increase_left_fit_bound_more, smaller_A ]

        
    # this covers the case in which we have already used all the functions as recorded in the 
    # db, so we don't bother to try again. the only way to get out of this is reset the function
    # number in db or add another function (latter almost certainly what you need).
    # note that we don't update the db
    if( attempt >= len(modifier_functions ) ):
        return 0

    
    # this checks if we had not previously attempted all functions, but none of them succeeded
    # on this run. in that case put it in the DB. 
    if attempt == len(modifier_functions):
        
        if UPDATE_DB and sql_conn is not None:
            if pixel_coords is None or fit_id is None:

                msg = 'WARNING (%d, %d, %d): unable to insert bad fit attempt into db, \
                        pixel_coords or fit_id not specified.' % (pixel_coords[0], pixel_coords[1], fit_id ) 

                print( msg )

                if log_enabled:
                    log_message( logfile, msg )

                    return 0


                # end of log handling: add to the db
                db.insert_fit_data_into_db( sql_conn, pixel_coords, \
                    fit_id, 0, fit_attempt=attempt )

            # return 0 if last available attempt failed.
            return 0

    
    # otherwise we call the current attempt function. 
    modifier_functions[attempt]( p0, p0_attempt, fit_bounds, fit_bounds_attempt, npeaks )
    
    if PRINT_FIT_ATTEMPT:
        print( 'INFO: testing fit attempt ' + str(attempt) )
        print( 'p0_attempt = ' + str(p0_attempt) )
        print( 'fit_bounds_attempt = ' + str(fit_bounds_attempt) )
        
    return 1
    





# if db connection supplied, then we assume that the entry is supposed to be put in, 
# it makes more sense for the check on whether the data is in the db to be performed 
# before calling this function. when called: 
# apply peak fit, make sure it worked, try again if not, put it in db if successful,
# and then plot. return 1 if successful fit obtained, otherwise 0. last_attempt is the 
# last attempted + 1, i.e. it is the first fit that will be attempted. 

def apply_peak_fit( ax, pixel_coords, fit_id, x, y, fit_bounds,
                    p0, npeaks, last_attempt=-1, mu_all=None, 
                    muerr_all=None, reduc_chisq_all=None, sql_conn=None, peak_detect=None ):
    

    fitfunc = dlf.n_fitfuncs_abstract( dlf.fitfunc_n_alpha_peaks, npeaks )

    
    # these will be populated based on the current attempt number and the given p0 and fit_bounds,
    # which are used as the first guesses. they are passed as parameters so they can be used after
    # without returning them.
    p0_attempt = []
    fit_bounds_attempt = []

    
    # start attempt is the one after the most recently tested attempt.
    current_attempt = last_attempt + 1

    
    while( 1 ):

        all_attempts_failed = next_fit_attempt( p0, p0_attempt,
                                                fit_bounds, fit_id, fit_bounds_attempt, 
                                                npeaks, current_attempt, sql_conn=sql_conn,
                                                pixel_coords=pixel_coords )
        
        if not all_attempts_failed:
            
            if PRINT_FIT_STATUS:
                print( "WARNING: peak failed to converge." )

            return 0

        
        ret = jmath.jacob_least_squares( x, y, np.sqrt(y),
                                         p0_attempt, fitfunc, fit_bounds=fit_bounds_attempt,
                                         reduc_chisq_max = 6 )
        
        if ret is not None:

            if PRINT_FIT_STATUS:
                print( "INFO: success, breaking " )
            break

        
        else:
            if PRINT_FIT_STATUS:
                print( "WARNING: Unsuccessful fit, trying new parameters..." )
            current_attempt += 1


    # now we have broken out of the loop after a successful fit. unpack the results.
    reduc_chisq, dof, pf, pferr = ret         
                              
    # write all relevant data to the db
    if UPDATE_DB and sql_conn is not None:
        db.insert_fit_data_into_db( sql_conn, pixel_coords, fit_id,
                                    1, last_attempt, reduc_chisq,
                                    pf, pferr, p0_attempt, fit_bounds_attempt,
                                    peak_detect ) 

    # where we are after breaking    
    if PRINT_PF:
        print( "pf = " + str(pf) )
    
    if PRINT_PFERR:
        print( "pferr = " + str(pferr) )

        
    jplt.add_fit_to_plot( ax, x, fit_bounds_attempt, pf, pferr, fitfunc )


    # extend vectors as necessary. done here since this is the last time
    # we get access to such variables unless taking such an approach.
    if mu_all is not None:
        mu_all.extend( [pf[ 5:3+2*npeaks:2 ] ]  )

    if muerr_all is not None:
        muerr_all.extend( [pferr[ 5:3+2*npeaks:2 ] ] )

    if reduc_chisq_all is not None:
        reduc_chisq_all.append( reduc_chisq )


        
    return 1
    
    





# do all behavior specific to our application: read file, fit peaks, make the plot.
# if sql_conn is supplied, they we assume that it is a valid open sqlite connection
# and check db for whether this file has been processed before.
def process_file( ax, infile, pixel_coords, sql_conn=None, logfile=None ):

    log_enabled = logfile is not None
    
    coordsstr = "(%d, %d)" % pixel_coords

    # x axis, needed no matter what.
    x = np.arange( 5000 )


    # create and populate the histogram array 
    efront_histo = np.zeros( x.size )

    
    if not dlf.construct_histo_array( infile, efront_histo ):

        msg = "ERROR: unable to open file: " + infile 

        print( msg )
        
        if log_enabled:
            log_message( logfile, msg )

        sys.exit(0)
    
    
    # get the 6 suspected peaks
    our_peaks = [0] * NUM_PEAKS
    num_peaks_found = dlf.get_n_peak_positions( NUM_PEAKS, efront_histo, our_peaks )

    if num_peaks_found < 5:

        msg = 'WARNING (%d,%d,*): found less than 5 peaks, skipping this fit.'
        msg %= (pixel_coords[0], pixel_coords[1] )

        print( msg )

        if log_enabled:
            log_message( logfile, msg )
        
        if sql_conn is not None:
            for fit_id in range(3):
                db.insert_fit_data_into_db( sql_conn, pixel_coords, fit_id,
                                            successful_fit = 0 )


        return 0

    
    # debug 
    if PRINT_PEAK_POSITIONS:
        print( "INFO: found peaks " + str(our_peaks) )
    
                
    # set the x axis according to the peak position 
    plot_bounds = [our_peaks[0] - 60, our_peaks[-1]+100 ]

 
    # plot the histogram without fit yet 
    jplt.plot_histo( ax, x, efront_histo, plot_bounds=plot_bounds, logscale=0,
                     title = "", 
                     xlabel = "", 
                     ylabel = "" )


    if SHOW_HISTO:
        plt.show()



    # do a check on peak values: energy differenc for the pairs should
    # be constant. no check necessary on the largest peak since it
    # always dominates, the only potential problem is really the
    # smallest peak in the lowest energy pair. this occurs after
    # plotting so that we can return with the plot alread made.
    # basically, we are saying that if we cannot determine where the
    # peak positions are to 0th order, we are not going to bother
    # fitting them since it will definitely not work if we misidentify
    # a peak position.

    check_warnings = 1
    
    if check_warnings:
        if( abs(23 - (our_peaks[1] - our_peaks[0] ) ) > 45 ):

            msg = ( "WARNING (%d,%d,*): invalid peak suspected for pair 1. Skipping. "
                         % (pixel_coords[0], pixel_coords[1]) )

            print( msg ) 

            if log_enabled:
                log_message( logfile, msg )

            return 0
                         
    
        if( abs(23 - (our_peaks[3] - our_peaks[2] ) ) > 45 ):

            msg = ( "WARNING (%d,%d,*): invalid peak suspected for pair 2. Skipping"
                         % (pixel_coords[0], pixel_coords[1]) ) 

            print( msg ) 
                         
            if log_enabled:         
                log_message( logfile, msg )

            return 0



        
        # determine which fits to perform, 1 = fit must be attempted. by default if no 
    # conn is supplied we process all the fits.

    fit_attempts = [0, 0, 0 ] # first fit to try, which is the one we left off on.

    fits_to_perform = [ 1, 1, 1 ]

    #    mu_all = []  # this shall store the mu values 
    # muerr_all = []
    reduc_chisq_all = [] # store all reduced chisq values.
    
    
    if sql_conn is not None:
        for i in range(len(fits_to_perform)):
            
            # extract
            result = db.read_data_from_db( sql_conn, pixel_coords, i )
            successful_fit, fit_attempt, reduc_chisq, pf, pferr, p0, fit_bounds, peak_detect = result
            
            fits_to_perform[i] = not successful_fit
            
            fit_attempts[i] = fit_attempt
            
            if successful_fit:
                fitfunc = dlf.n_fitfuncs_abstract( dlf.fitfunc_n_alpha_peaks, NUM_PEAKS_PER_FEATURE[i] )
                jplt.add_fit_to_plot( ax, x, fit_bounds, pf, pferr, fitfunc )
                reduc_chisq_all.append( reduc_chisq )
                
                
        if not any( fits_to_perform ):
            return 1
                
 
    
        
    # number of peaks for each fit, need this info before plotting
    # anything so it is defined here. depends on if we found 5 or 6 peaks
    if num_peaks_found == 5:
        num_peaks = [ 2, 2, 1 ] 

    elif num_peaks_found == 6:
        num_peaks = [ 2, 2, 2 ]

        
    
   
    
    # initial fit bounds, these are reduced a bit if the fit fails.
    # the last one is a hack that only works in this particular situation.
    fit_bounds = [\
                    [ our_peaks[0]-50, our_peaks[1]+13 ],  \
                    [ our_peaks[2]-30, our_peaks[3]+13 ],  \
                    [ our_peaks[4]-80, our_peaks[-1]+16 ] \
                ]
    
    # initial guesses which get modified later if the fit fails.
    # the final term depends on whether we found 5 or 6 peaks.
    p0 = [
            [ 6.0, 0.97, 20.2, 2.0 ] + [ 20000.0, our_peaks[0] + 8.0 ] + [ 60000.0, our_peaks[1] + 8.0 ],
            [ 6.0, 0.99, 42.0, 1.6 ] + [ 50000.0, our_peaks[2] + 8.0 ] + [ 100000.0, our_peaks[3] + 8.0 ],
    ]

    if num_peaks_found == 5:
        p0 += [ [ 4.0, 0.99, 30.0, 1.0 ] + [ 200000.0, our_peaks[4]+8 ] ]

    elif num_peaks_found == 6:
        p0 += [ [ 6.0, 0.99, 42.0, 1.6 ] + [ 50000.0, our_peaks[4] + 8.0 ] + [ 200000.0, our_peaks[5] + 8.0 ] ] 


    peak_detect = [ our_peaks[0:2], our_peaks[2:4], our_peaks[4:] ]
        
    # loop through the fits that were not in the db and add them if successful.
    for i in range( NUM_FEATURES ):
        if fits_to_perform[i]:
                   
            apply_peak_fit( ax, pixel_coords, i, x, efront_histo, fit_bounds[i], p0[i],
                            num_peaks[i], fit_attempts[i], peak_detect[i],
                            reduc_chisq_all = reduc_chisq_all, sql_conn = sql_conn,
                            peak_detect = peak_detect[i] ) #, mu_all,
                            # muerr_all, reduc_chisq_all, sql_conn=sql_conn )
    
    
    # # add a description of the fit to the plot
    # # note: \mathrm{} is replacement for \text{}
    # fitstr = "$ f(E) = \\sum_{i,\pm} \\frac{A_i \eta_\pm}{2 \\tau_\pm}    \
    #                     \\cdot \\exp \\left[   \
    #                             \\frac{E-\\mu_i}{\\tau_\pm}   \
    #                             + \\frac{\\sigma^2}{2 \\tau_\pm^2 }  \
    #                     \\right] \
    #                     \\cdot \\mathrm{erfc} \\left[   \
    #                             \\frac{1}{\\sqrt{2}} \\left(  \
    #                                     \\frac{x-\\mu}{\\sigma}    \
    #                                     + \\frac{\\sigma}{\\tau_\pm}   \
    #                             \\right)   \
    #                     \\right] $"
    
    
    
    
    
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
    reduc_chisq_str = jutils.format_measurement_vector( "\\tilde{\\chi}^2", reduc_chisq_all, 0.01 )
    
    # fitstr += '\n' + '\n'.join( [mu_str, reduc_chisq_str] )
#    fitstr += '\n' + reduc_chisq_str
    fitstr = '\n' + reduc_chisq_str
    #ax.text( 0.02, 0.92, fitstr, transform=ax.transAxes, fontsize=12, verticalalignment='top')
    #ax.text( 0.02, 0.93, coordsstr, transform=ax.transAxes, fontsize=12, verticalalignment='top' )
    
    ax.text( 0.02, 1.05, fitstr, transform=ax.transAxes, fontsize=12, verticalalignment='top')
    ax.text( 0.03, 0.80, coordsstr, transform=ax.transAxes, fontsize=12, verticalalignment='top' )
        
    return 1





# main function, goes through all the data.
def fit_all_peaks():

    # dimesions of detector
    totalx = 32
    totaly = 32

    # dimensions of output images
    dimx = 4
    dimy = 4

    # used to name the image files
    files = db.all_db_ids

    databases = db.all_dbs_list

    # used for estimating time remaining.
    start_time = time.time()
  
    # loop through each file prefix / corresponding database.
    for a in range( len(files) ):

        # create the current db if it doesn't already exist.
        if not os.path.exists( databases[a] ):
            db.create_db( databases[a] )

        # this is the name associated with the current db, e.g. 'centered'.
        fileprefix = files[a]

        # make dir for output images if it doesn't exist 
        current_outdir = '../../images/current_fit_images/' + fileprefix + '/'

        if not os.path.exists( current_outdir ):
            os.mkdir( current_outdir )


        logfile = current_outdir + 'log.txt'

        with sqlite3.connect( databases[a] ) as sql_conn:
            
            # these 2 loops loop over grids of dimx x dimy images. 
            for x in range( totalx ):
                
                # these loops create the 4x4 grids. 2 per strip row.
                for k in range( totaly // (dimx * dimy ) ):

                    plt.clf()
                    plt.close()
                    f, axarr = plt.subplots( dimx, dimy, figsize=(20,10) )    
                    
                    for i in range(dimx):
                        for j in range(dimy):

                            y = ( i * dimx + j ) + ( k * dimx * dimy )
                                
                            print( str( [x,y] ) )    
                                
                            # this estimates the time remaining for the program to terminate
                            jutils.estimate_time_left( x * totaly + y
                                                       + ( a * totalx * totaly ), 
                                                       len(files) * totalx * totaly,
                                                       start_time, num_updates=100 )
                                
                            # current pixel coords
                            coords = ( x, y )
                            
                            # file containing the data
                            current_file = fileprefix + "_%d_%d.bin" % coords
                            
                            if PRINT_FILE_NAMES:
                                print( "INFO: processing file: " + current_file )
                                
                            current_file = ( "../../data/extracted_ttree_data/"
                                             + fileprefix + '/' + current_file )
                                
                            process_file( axarr[i,j], current_file, coords, sql_conn = sql_conn, logfile = logfile)
                                
                            if( BREAK_AFTER_1_PLOT ):
                                plt.show()
                                return 1
                                
                    jplt.saveplot_low_quality( current_outdir, fileprefix + '_x=%d.%d' % (x, k) )


                        
    # plt.setp([a.get_xticklabels() for a in axarr[0, :]], visible=False)
    # plt.setp([a.get_yticklabels() for a in axarr[:, 1]], visible=False)
    


    
    
    
# this function is to be used for debugging the fits. once you have found a fit that 
# does not converge, add a new fit and 
def make_one_plot( db_id, x, y ):
    
    if db_id not in db.all_db_ids:
        print( 'ERROR: db_id given is not in the list of db ids. ' )
        return 0
    
    set_globals_for_debugging()
    
    current_file = db_id +  "_%d_%d.bin" % (x,y)
    
    if PRINT_FILE_NAMES:
        print( "INFO: processing file: " + current_file )
    
    current_file = ( "../../data/extracted_ttree_data/"
                     + db_id + '/'  + current_file )
        
    plt.clf()
    ax = plt.axes()
    
    with sqlite3.connect( db.all_dbs_dict[ db_id ]  ) as sql_conn:
        
        process_file( ax, current_file, (x,y), sql_conn = sql_conn, logfile = None )

    plt.show()

    return 1
    
    
       


# print message to the console and write it to a log file.
def log_message( logfile, msg ):
    
    if log_message.first_msg or not os.path.exists( logfile ):
        log_message.first_msg = 0
        mode = 'w'
    else:
        mode = 'a'
        
    with open( logfile, mode ) as log:
        log.write( msg + '\n' )
    
log_message.first_msg = 1




# set debug globals, called by make_one_plot()
def set_globals_for_debugging():

    global UPDATE_DB
    UPDATE_DB = 0
    
    global PRINT_PEAK_POSITIONS
    PRINT_PEAK_POSITIONS = 1
    
    global PRINT_FIT_ATTEMPT
    PRINT_FIT_ATTEMPT = 1
    
    global PRINT_FIT_STATUS
    PRINT_FIT_STATUS = 1

    global SHOW_HISTO
    SHOW_HISTO = 1
     
       
# make_one_plot( db.moved_id, 2, 26 )
fit_all_peaks()


