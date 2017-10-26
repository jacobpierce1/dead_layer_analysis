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


# from peakdetect import peakdetect 
import deadlayer_helpers.sql_db_manager as dbmgr
import deadlayer_helpers.data_handler as data


import jspectroscopy as spec


# from lmfit import Model, Parameters 




## includes 
import numpy as np
import matplotlib.pyplot as plt
import time 
import sys
import sqlite3
import os
from enum import Enum 




NUM_PEAKS_PER_FEATURE = [ 2, 2, 2 ] 
NUM_PEAKS = np.sum( NUM_PEAKS_PER_FEATURE )
NUM_FEATURES = len( NUM_PEAKS_PER_FEATURE )




def make_model_and_apply_fit( npeaks, params_guess, x, y, dy ):
    pass


# these are the functions that attempt to modify p0_attempt and fit_bounds_attempt.
def default( p0, p0_attempt, fit_bounds, fit_bounds_attempt, npeaks ):
    pass 
    
def increase_left_fit_bound( p0, p0_attempt, fit_bounds, fit_bounds_attempt, npeaks ):
    fit_bounds_attempt[0] += 10
    
def increase_left_fit_bound_more( p0, p0_attempt, fit_bounds, fit_bounds_attempt, npeaks ):
    fit_bounds_attempt[0] += 20
    
def smaller_A( p0, p0_attempt, fit_bounds, fit_bounds_attempt, npeaks ):
    # print( npeaks ) 
    for i in range( npeaks ):
        p0_attempt['A' + str(i)] = 0.4 * p0[ 'A' + str(i) ]

def larger_A( p0, p0_attempt, fit_bounds, fit_bounds_attempt, npeaks ):
    # fit_bounds_attempt[1] += 5
    for i in range( npeaks ):
        p0_attempt[ 'A' + str(i) ] = 3 * p0[ 'A' + str(i) ]

def next_fit( p0, p0_attempt, fit_bounds, fit_bounds_attempt, npeaks ):
    pass



class fit_attempt_status( Enum ):
    SUCCESS = 1
    FAILURE = 2
    NO_UPDATE = 3



# this function creates another guess for p0 in the case of a failed fit based on
# observations about what is gonig wrong. the particular attempts we use depend on the peak id.
def next_fit_attempt( p0, p0_attempt, fit_bounds, fit_id,
                      fit_bounds_attempt, npeaks, attempt ):

    
    # these will be modified and used later as the initial guess.
    for key, val in p0.items():
        p0_attempt[ key ] = val     

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
        return fit_attempt_status.NO_UPDATE

    
    # this checks if we had not previously attempted all functions, but none of them succeeded
    # on this run. in that case put it in the DB. 
    if attempt == len(modifier_functions):
        return fit_attempt_status.FAIL

    
    # otherwise we call the current attempt function
    # to generate a new fit parameters attempt.
    modifier_functions[ attempt ] ( p0, p0_attempt, fit_bounds, fit_bounds_attempt, npeaks )
    
    if PRINT_FIT_ATTEMPT:
        print( 'INFO: testing fit attempt ' + str(attempt) )
        print( 'p0_attempt = ' + str(p0_attempt) )
        print( 'fit_bounds_attempt = ' + str(fit_bounds_attempt) )
        
    return 1
    





# if db connection supplied, then we assume that the entry is supposed
# to be put in, it makes more sense for the check on whether the data
# is in the db to be performed before calling this function. when
# called: apply peak fit, make sure it worked, try again if not, put
# it in db if successful, and then plot. return 1 if successful fit
# obtained, otherwise 0. last_attempt is the last attempted + 1,
# i.e. it is the first fit that will be attempted.

def apply_peak_fit( ax, x, y, fit_id,
                    xvals, yvals,
                    fit_bounds,
                    p0, npeaks, last_attempt=-1, 
                    params_bounds = None, reduc_chisq_all=None,
                    db=None, peak_detect=None ):
    
    
    fitfunc = spec.sum_n_fitfuncs( spec.fitfunc_n_alpha_peaks, npeaks )

    
    # these will be populated based on the current attempt number and
    # the given p0 and fit_bounds, which are used as the first
    # guesses. they are passed as parameters so they can be used after
    # without returning them.

    p0_attempt = {}
    fit_bounds_attempt = []

    
    # start attempt is the one after the most recently tested attempt.
    current_attempt = last_attempt + 1

    
    while( 1 ):

        status = next_fit_attempt( p0, p0_attempt,
                                   fit_bounds, fit_id, fit_bounds_attempt, 
                                   npeaks, current_attempt )

        # our fit attempts all failed in this case 
        if status == fit_attempt_status.FAILURE:

            
            if PRINT_FIT_STATUS:
                print( "INFO: fit failed, adding to DB " )
            
            if UPDATE_DB and db is not None:
                db.insert_fit_data( x, y, fit_id,
                                    last_attempt=attempt )

        # in this case, the fits all failed but it was already in DB 
        if status == fit_attempt_status.NO_UPDATE:
            
            if PRINT_FIT_STATUS:
                print( "INFO: peak failed to converge, already detected in db.." )

            return 0


        # if those 2 cases are not reached, then we have a new fit attempt.
        # try doing least squares. 
        model = jmath.jleast_squares( xvals, yvals, np.sqrt( yvals ),
                                      p0_attempt, fitfunc,
                                      fit_bounds = fit_bounds_attempt,
                                      reduc_chisq_max = 3.5,
                                      params_bounds = params_bounds,
                                      successful_fit_predicate = successful_alpha_fit_predicate,
                                      print_results = 0 )
        
        if model is not None:
            if PRINT_FIT_STATUS:
                print( "INFO: success, breaking " )
            break

        
        else:
            if PRINT_FIT_STATUS:
                print( "WARNING: Unsuccessful fit, trying new parameters...\n\n" )
            current_attempt += 1


    # now we have broken out of the loop after a successful fit. unpack the results.
    reduc_chisq = model.redchi

    # model.plot_fit( ax = ax, datafmt = '-r', numpoints = 100 * model.ndata )

    # dof = model.nfree
    # pf = model.params
    # pferr = model.         


    
    
    # write all relevant data to the db
    if UPDATE_DB and db is not None:
        db.insert_fit_data( x, y, fit_id,
                            successful_fit = 1,
                            npeaks = npeaks, 
                            last_attempt = current_attempt,
                            params_guess = p0_attempt,
                            fit_bounds = fit_bounds_attempt, 
                            peak_guesses = peak_detect,
                            model = model ) 
        

    # # where we are after breaking    
    # if PRINT_PF:
    #     print( "pf = " + str(pf) )
    
    # if PRINT_PFERR:
    #     print( "pferr = " + str(pferr) )


    # final_fitfunc = lambda x : model.eval( x=x ) 
    jplt.add_fit_to_plot( ax, xvals, fit_bounds_attempt, jmath.model_func( model ) )


    # # extend vectors as necessary. done here since this is the last time
    # # we get access to such variables unless taking such an approach.
    # if mu_all is not None:
    #     mu_all.extend( [pf[ 5:3+2*npeaks:2 ] ]  )

    # if muerr_all is not None:
    #     muerr_all.extend( [pferr[ 5:3+2*npeaks:2 ] ] )

    if reduc_chisq_all is not None:
        reduc_chisq_all.append( reduc_chisq )


        
    return 1
    
    




# only allow a successful fit if we have sub-2 channel precision
# on the mu values.

def successful_alpha_fit_predicate( params ):

    i = 0
    keys = params.keys()

    while( 1 ) :

        mu = 'mu' + str(i)
        
        if mu in keys: 
            if params[ mu ].stderr > 4:
                return 0
            i += 1 
            
        else:
            break
            
    return 1 






# do all behavior specific to our application: read file, fit peaks,
# make the plot.  if db is supplied, they we assume that it is a valid
# open sqlite connection and check db for whether this file has been
# processed before.

def process_file( ax, infile, x, y, db = None, logfile=None ):

    log_enabled = logfile is not None
    
    # x axis, needed no matter what.
    xaxis = np.arange( 5000 )


    # create and populate the histogram array 
    efront_histo = np.zeros( xaxis.size )

    
    if not data.construct_histo_array( infile, efront_histo ):

        msg = "ERROR: unable to open file: " + infile 

        print( msg )
        
        if log_enabled:
            log_message( logfile, msg )

        sys.exit(0)
    
    
    # get the 6 suspected peaks
    our_peaks = [0] * NUM_PEAKS
    our_peaks, num_peaks_found = jmath.get_n_peak_positions( NUM_PEAKS, efront_histo )

    if num_peaks_found < 5:

        msg = 'WARNING (%d,%d,*): found less than 5 peaks, skipping this fit.'
        msg %= ( x, y )

        print( msg )

        if log_enabled:
            log_message( logfile, msg )
        
        if db is not None:
            for fit_id in range(3):
                db.insert_fit_data( x, y, fit_id, successful_fit = 0 )


        return 0

    
    # debug 
    if PRINT_PEAK_POSITIONS:
        print( "INFO: found peaks " + str(our_peaks) )
    
                
    # set the x axis according to the peak position 
    plot_bounds = [our_peaks[0] - 60, our_peaks[-1]+100 ]

 
    # plot the histogram without fit yet 
    jplt.plot_histo( ax, xaxis, efront_histo, plot_bounds=plot_bounds, logscale=0,
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
                         % ( x, y ) )

            print( msg ) 

            if log_enabled:
                log_message( logfile, msg )

            return 0
                         
    
        if( abs(23 - (our_peaks[3] - our_peaks[2] ) ) > 45 ):

            msg = ( "WARNING (%d,%d,*): invalid peak suspected for pair 2. Skipping"
                         % ( x, y ) ) 

            print( msg ) 
                         
            if log_enabled:         
                log_message( logfile, msg )

            return 0

        # make sure that we didn't detect the small alpha peak after the main
        # highest energ peak. to do this, make sure that peak 5 has a higher
        # number of counts that peak 4

        if num_peaks_found == 6:
            if efront_histo[ our_peaks[5] ]  < efront_histo[ our_peaks[4] ]:
                our_peaks = our_peaks[0:5]
                num_peaks_found = 5
                
            
        
        
    # determine which fits to perform, 1 = fit must be attempted. by default if no 
    # conn is supplied we process all the fits.

    fit_attempts = [ -1, -1, -1 ] # first fit to try, which is the one we left off on.

    fits_to_perform = [ 1, 1, 1 ]

    #    mu_all = []  # this shall store the mu values 
    # muerr_all = []
    reduc_chisq_all = [] # store all reduced chisq values.
    

    # option to pass None as the DB 
    if db is not None:

        for i in range(len(fits_to_perform)):
            
            # extract
            db_data = db.read_fit_data( x, y, i )

            successful_fit = db_data[ 'successful_fit' ]

            # signal that this fit will have to be re-attempted.
            fits_to_perform[i] = not successful_fit

            # last attempt that we left off on 
            fit_attempts[i] = db_data[ 'last_attempt' ]
            
            if successful_fit:
                # fitfunc = spec.sum_n_fitfuncs( spec.fitfunc_n_alpha_peaks, NUM_PEAKS_PER_FEATURE[i] )
                
                model = db_data[ 'model' ] 
                jplt.add_fit_to_plot( ax, xaxis, db_data[ 'fit_bounds' ], jmath.model_func( model ) )
                reduc_chisq_all.append( model.redchi )
                

        # no further processing required if all the fits converged.
        if not any( fits_to_perform ):
            _add_text( ax, x, y, reduc_chisq_all )
            return 1
                
 
    
        
    # number of peaks for each fit, need this info before plotting
    # anything so it is defined here. depends on if we found 5 or 6 peaks
    if num_peaks_found == 5:
        num_peaks = [ 2, 2, 1 ] 

    elif num_peaks_found == 6:
        num_peaks = [ 2, 2, 2 ]

        
    
    # construct guesses for the fitbounds and parameters.
        
    fit_bounds = [\
                    [ our_peaks[0]-50, our_peaks[1]+13 ], 
                    [ our_peaks[2]-30, our_peaks[3]+13 ], 
                    [ our_peaks[4]-80, our_peaks[-1]+16 ]  # hack over here 
                ]
    
    sigma_guess = 5.0
    eta_guess = 0.97
    tau1_guess =  35.0
    tau2_guess = 1.5

    # construct guesses for the A parameter.
    A_array_guess = [ [ 20000.0, 60000.0 ], [ 50000.0, 100000. ] ]
    mu_array_guess = [ [ 8 + our_peaks[ 2*i ], 8 + our_peaks[ 2*i + 1 ] ] for i in [0, 1]  ]
    
    if num_peaks_found == 6:
        A_array_guess.append( [ 10000.0, 200000.0 ] )
        mu_array_guess.append( [ our_peaks[4] + 8, our_peaks[5] + 8 ] )
    else: 
        A_array_guess.append( [ 200000.0 ] ) 
        mu_array_guess.append( [ our_peaks[4] + 8 ] )
        

    # these store the initial parameter guesses and 
    # bounds on each parameter ( slightly inefficient ).

    p0 = [0] * NUM_FEATURES
    params_bounds = [0] * NUM_FEATURES 

    for i in range( NUM_FEATURES ):

        p0[i], params_bounds[i] = spec.construct_n_alpha_peaks_params(
            sigma_guess, eta_guess, tau1_guess,
            tau2_guess, A_array_guess[i],
            mu_array_guess[i] )


    # list of detected peaks, to be added to DB (not yet implemented) 
    peak_detect = [ our_peaks[0:2], our_peaks[2:4], our_peaks[4:] ]

    
    # loop through the fits that were not in the db and add them if successful.
    for i in range( NUM_FEATURES ):
        if fits_to_perform[i]:

            #            print( '\nINFO: fitting feature ' + str(i) + '...' ) 
            
            apply_peak_fit( ax, x, y, i,
                            xaxis, efront_histo,
                            fit_bounds[i], p0[i],
                            num_peaks[i], fit_attempts[i],
                            reduc_chisq_all = reduc_chisq_all, db = db,
                            peak_detect = peak_detect[i], 
                            params_bounds = params_bounds[i] ) 
    
    
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
    

    _add_text( ax, x, y, reduc_chisq_all ) 
    
        
    return 1






def _add_text( ax, x, y, reduc_chisq_all ):
    
    coordsstr = "(%d, %d)" % ( x, y )

    reduc_chisq_str = jutils.format_measurement_vector( "\\tilde{\\chi}^2", reduc_chisq_all, 0.01 )

    fitstr = '\n' + reduc_chisq_str
    
    ax.text( 0.02, 1.05, fitstr, transform=ax.transAxes, fontsize=12, verticalalignment='top')
    ax.text( 0.03, 0.80, coordsstr, transform=ax.transAxes, fontsize=12, verticalalignment='top' )


    



# main function, goes through all the data.
def fit_all_peaks():

    # dimesions of detector
    totalx = 32
    totaly = 32

    # dimensions of output images
    dimx = 4
    dimy = 4

    a = -1

    start_time = time.time() 

    for db in dbmgr.all_dbs:

        # count number of db's we have processed. 
        a += 1 
        print( db.name )
        

        # construct db if not already existing 
        if not db.exists():
            db.create()

        db.connect()

       
        # make dir for output images if it doesn't exist 
        current_outdir = '../../images/current_fit_images/' + db.name + '/'

        if not os.path.exists( current_outdir ):
            os.mkdir( current_outdir )


        logfile = current_outdir + 'log.txt'
        
        
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
                                                   len( dbmgr.all_dbs ) * totalx * totaly,
                                                   start_time, num_updates=100 )
                        
                        # current pixel coords
                        # coords = ( x, y )
                        
                        # file containing the data
                        current_file = db.name + "_%d_%d.bin" % ( x, y )
                        
                        if PRINT_FILE_NAMES:
                            print( "INFO: processing file: " + current_file )
                            
                        current_file = ( "../../data/extracted_ttree_data/"
                                         + db.name + '/' + current_file )
                        
                        process_file( axarr[i,j], current_file, x, y,
                                      db = db, logfile = logfile ) # sql_conn = sql_conn, logfile = logfile)
                        
                        if( BREAK_AFTER_1_PLOT ):
                            plt.show()
                            return 1
                        
                jplt.saveplot_low_quality( current_outdir, db.name + '_x=%d.%d' % (x, k) )
                            
            
        db.disconnect()
                        # plt.setp([a.get_xticklabels() for a in axarr[0, :]], visible=False)
                            # plt.setp([a.get_yticklabels() for a in axarr[:, 1]], visible=False)
                            
                            

    
    
    
# this function is to be used for debugging the fits. once you have found a fit that 
# does not converge, add a new fit and 
def make_one_plot( db, x, y, test_db = 0 ):
    
    if db not in dbmgr.all_dbs :
        print( 'ERROR: db_id given is not in the list of db ids. ' )
        return 0
    
    set_globals_for_debugging( test_db )
    
    current_file = db.name +  "_%d_%d.bin" % (x,y)
    
    if PRINT_FILE_NAMES:
        print( "INFO: processing file: " + current_file )
    
    current_file = ( "../../data/extracted_ttree_data/"
                     + db.name + '/'  + current_file )
        
    ax = plt.axes()
    

    if test_db:

        # construct db if not already existing 
        if not db.exists():
            db.create()

        db.connect()
        process_file( ax, current_file, x, y, db = db )
        db.disconnect()

    else:
        process_file( ax, current_file, x, y ) 

        
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
def set_globals_for_debugging( test_db = 0 ):

    if not test_db:
        global UPDATE_DB
        UPDATE_DB = 0
    
    global PRINT_PEAK_POSITIONS
    PRINT_PEAK_POSITIONS = 1
    
    global PRINT_FIT_ATTEMPT
    PRINT_FIT_ATTEMPT = 1
    
    global PRINT_FIT_STATUS
    PRINT_FIT_STATUS = 1

    global SHOW_HISTO
    SHOW_HISTO = 0
     
       
# make_one_plot( dbmgr.centered, 16, 16, test_db = 1 )
fit_all_peaks()


