# fit all peaks for data with Cm and Gd sources


import matplotlib
matplotlib.use('Agg')
# -*- coding: utf-8 -*-

import jspectroscopy as spec
import numpy as np
#import deadlayer_helpers.data_handler as data
import matplotlib.pyplot as plt 
import jutils

import bpt
import sys


debug = 0



# inputs = spectrum_fitter_inputs( (32,32), data_fetcher )

# num_groups = 3

rel_plot_bounds = [ -100, 120 ] 


group_ranges = np.array( [ [-95,20], [-70,20] ] )
peak_locations = [ [-63,0], [-21, 0] ]
peak_sizes_guesses = [ [300, 15000], [3000, 15000] ]

peak_types = [ ['a'] * len(x) for x in peak_locations ]

peak_width_guesses = [ [5.0] * len(x) for x in peak_locations ]

det_params_guesses = [ { 'a' : [ 0.95, 40.0, 2.4 ] } ] * 2

peak_mu_offset = 7.0

num_peaks_to_detect = 2


coords = ( 0, 28, 25 ) 




def peak_params_guesser( i, primary_peak, chans, counts ) :

    # gd
    peak_params = np.zeros((2,3))
    peak_params[:,2] = 5.0

    # print( primary_peak ) 
        
    if i == 0 :
        counts = np.sum( counts[ ( chans > ( primary_peak - 140 ) )
                               & ( chans < ( primary_peak + 32 ) ) ] )
        peak_params[:,0] = 2 * np.array( [0.02, 0.98] ) * counts
        peak_params[:,1] = primary_peak + peak_locations[i] + peak_mu_offset

    # cm 
    elif i == 1 :
        counts = np.sum( counts[ ( chans > primary_peak - 110 )
                                 & ( chans < primary_peak + 50 ) ] )
        peak_params[:,0] = 2 * np.array( [0.23, 0.77] ) * counts
        peak_params[:,1] = primary_peak + peak_locations[i] + peak_mu_offset

    # print( peak_params ) 
    
    return peak_params




def fit_bounds_guesser( i, primary_peak, chans, counts ) :

    if i == 0 :
        return primary_peak + group_ranges[0]

    else :
        return primary_peak + group_ranges[1] 




    
# a function that is guarenteed to detect the reference peak
# of each group given a list of the peaks detected in the spectrum

def primary_peak_detector( peaks_detected, counts ) :
    
    if debug : 
        print( peaks_detected ) 

    if len( peaks_detected ) < 2 :
        return None

    if( debug ) :
        print( 'peaks_detected', peaks_detected )

    if ( peaks_detected[0] > 2000
         or peaks_detected[1] < 2000
         or peaks_detected[0] < 1000 ) :

        if debug :
            print( "peak detect failed" )
        
        return None
    
    return peaks_detected 





def params_shuffler() :
    return 1 



    

def fit_acceptor( x, y, dy, spec_fitter_result ) :

    if spec_fitter_result.pvalue < 0.05 :
        return 0 

    mu_values = [ spec_fitter_result.peak_params[i][1]
                  for i in range( spec_fitter_result.num_peaks ) ]

    # print( mu_values )
    # print( sorted( mu_values ) )
    
    if sorted( mu_values ) != mu_values :
        return 0

    mu_guess = [ spec_fitter_result.peak_params_guess[i][1]
                 for i in range( spec_fitter_result.num_peaks ) ]

    if spec_fitter_result.group_num == 2 : 
        if np.any( np.abs( np.array( mu_values ) - np.array( mu_guess ) ) > 7 ) :
            return 0 

    # make sure mu values do not differ much from guess
    
    return 1





def one_spectrum( db_name, detnum, x, y ) :
    
    db = spec.spectrum_db( db_name, '../../../storage/' )

    global debug
    debug = 1
        
    chans, counts, dcounts = db.get_histo( detnum, x, y )
    
    dcounts = np.sqrt( counts )
    dcounts[ ( dcounts == 0 ) ] = 1 
    
    plt.figure(figsize=(10,12))
    ax = plt.axes()

    primary_peaks = db.load_dill( 'primary_peaks' )
    primary_peaks = primary_peaks[ :, detnum, x, y ] 

    db.disconnect()

    # peak_params_guesses = [ peak_params_guesser( i, detnum, x, y )
    #                         for i in range(2) ] 

    spec.auto_fit_spectrum( chans, counts, dcounts,
                            fit_bounds_guesser,
                            peak_params_guesser, det_params_guesses,
                            fit_acceptor = fit_acceptor,
                            params_shuffler = params_shuffler,
                            primary_peaks = primary_peaks,
                            rel_plot_bounds = rel_plot_bounds,
                            ax = ax, print_output = 1  )
    
    # spec.auto_fit_spectrum( x, y, dy,
    #                         group_ranges, peak_locations,
    #                         num_peaks_to_detect, primary_peak_detector,
    #                         peak_sizes_guesses, peak_width_guesses,
    #                         det_params_guesses, peak_mu_offset,
    #                         fit_acceptor = fit_acceptor,
    #                         params_shuffler = params_shuffler,
    #                         ax = ax,
    #                         rel_plot_bounds = rel_plot_bounds,
    #                         print_output = 1)
    

    plt.show()



    

def all_spectra( db_names ) :
    
    constrain_det_params = { 'a' : 1 }

    time_estimator = jutils.time_estimator( len(db_names) * 4 * 32 * 32, 20 )

    for name in db_names :

        # dets_used = bpt.dets_used( bpt_data_path, name )

        db = spec.spectrum_db( name, '../../../storage/', (32,32),
                               peak_types, constrain_det_params )

        mask = np.zeros((3,32,32))
        mask[0,0:15,:] = 1
        mask[0,23:,:] = 1 

        spec.auto_fit_many_spectra( db, (4,4),
                                    # group_ranges, peak_locations,
                                    # num_peaks_to_detect, primary_peak_detector,
                                    # peak_sizes_guesses, peak_width_guesses,
                                    # det_params_guesses,
                                    fit_bounds_guesser,
                                    peak_params_guesser, det_params_guesses,
                                    fit_acceptor = fit_acceptor,
                                    params_shuffler = params_shuffler,
                                    rel_plot_bounds = rel_plot_bounds,
                                    logscale = 1,
                                    time_estimator = time_estimator,
                                    print_output = 0,
                                    debug_peaks = 0,
                                    overwrite = 0,
                                    mask = mask )
        
        # mu_path = '../../storage/mu_values/%s_%d_mu_values.bin' % ( name, detnum ) 

            
        db.disconnect()





# all_spectra( [ 'alpharun20-30', 'alpharun11-19' ] ) 
# one_spectrum( 'full_bkgd_tot', 1, 3, 28 );



# one_spectrum( 'full_bkgd_tot', *coords );
all_spectra( [ 'full_bkgd_tot' ] )

# fit_all_spectra( dbmgr.all_dbs )


