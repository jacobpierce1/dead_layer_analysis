# -*- coding: utf-8 -*-

import jspectroscopy as spec
import numpy as np
import deadlayer_helpers.data_handler as data
import matplotlib.pyplot as plt 
import libjacob.jutils as jutils

import bpt


# inputs = spectrum_fitter_inputs( (32,32), data_fetcher )

# num_groups = 3

rel_plot_bounds = [ -100, 120 ] 


# group_ranges = [ [ -65, 15 ], [-75, 16], [-70,85] ] 


# peak_locations = [ [ -50, -20, 0 ], [ -57, -20, 0 ], [ -50, -30, 0, 22, 46, 64 ] ]

# peak_sizes_guesses = [ [ 2000., 30000.0, 90000.0 ],
#                        [ 3000., 50000.0, 100000. ],
#                        [ 5000.0, 10000.0, 100000.0, 100.0, 500.0, 500.0 ] ]

peak_locations = [ [ -42, -20, 0 ], [ -57, -20, 0 ], [ -30, 0, 22, 46, 64 ] ]

peak_sizes_guesses = [ [ 2000., 30000.0, 90000.0 ],
                       [ 3000., 50000.0, 100000. ],
                       [ 10000.0, 100000.0, 100.0, 500.0, 500.0 ] ]

group_ranges = [ [ -65, 15 ], [-75, 16], [-40, 85] ]

peak_types = [ ['a'] * len(x) for x in peak_locations ]

peak_width_guesses = [ [5.] * len(x) for x in peak_locations ]

det_params_guesses = [ { 'a' : [ 0.97, 35.0, 1.5 ] } ] * 3 

peak_mu_offset = 8 # 9 

num_peaks_to_detect = 6

# primary_peak_ids = None

# self.peak_structures = None 



                   






# a function that is guarenteed to detect the reference peak
# of each group given a list of the peaks detected in the spectrum

def primary_peak_detector( peaks_detected, histo ) :

    # print( peaks_detected ) 

    if len( peaks_detected ) < 6 :
        return None

    indices = np.empty(3, dtype = int )

    ret = np.empty( 3, dtype = float ) 
    
    # primary peaks 0 and 2 are the largest peaks from the left and right, respectively.

    # we are counting on detection of the first peak. without that
    # all is lost
    
    if abs( peaks_detected[1] - peaks_detected[0] ) > 45 :
        return None 
        
    indices[0] = 1
    
    peak0 = peaks_detected[1]
    ret[0] = peak0
    
    next_peak_offsets = [ 174, 320 ] 
    
    for i in range(2) :
        ret[ i+1 ] = ( peak0 + next_peak_offsets[i] - 15
                       + np.argmax( histo[ peak0 + next_peak_offsets[i] - 15 :
                                           peak0 + next_peak_offsets[i] + 15 ] ) )
    
    # print( ret )
    
    # if peaks_detected[3] - peaks_detected[2] < 34 :

    #     if histo[ peaks_detected[3] ] - histo[ peaks_detected[2] ] < 5 :
    #         indices[1] = 3
            
    #     else :

    #         if histo[ peaks_detected[4] ] - histo[ peaks_detected[3] ] < 5 : 
    #             indices[1] = 4
    #             indices[2] = 5

    #         else :
    #             return None
            
    # else :
    #     indices[1] = 2
    #     indices[2] = 3 + np.argmax( histo[ peaks_detected ][3:] ) 

    # ret = peaks_detected[ indices ] 

    # if ret[1] - ret[0] > 250 :
    #     return None

    # print( ret )
    
    return ret 





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




def one_spectrum( coords ) : 
    x, y, dy = data_fetcher( 'angled', *coords ) 


    dy = np.sqrt( y )
    dy[ ( dy == 0 ) ] = 1 
    
    plt.figure(figsize=(10,12))
    ax = plt.axes()     

    spec.auto_fit_spectrum( x, y, dy,
                            group_ranges, peak_locations,
                            num_peaks_to_detect, primary_peak_detector,
                            peak_sizes_guesses, peak_width_guesses,
                            det_params_guesses, peak_mu_offset,
                            fit_acceptor = fit_acceptor,
                            params_shuffler = params_shuffler,
                            ax = ax,
                            rel_plot_bounds = rel_plot_bounds,
                            print_output = 1)
    

    plt.show()



    

def all_spectra( db_names ) :

    constrain_det_params = { 'a' : 1 }

    time_estimator = jutils.time_estimator( len(db_names) * 32 * 32, 20 )

    for name in db_names : 
        
        db = spec.spectrum_db( '../../storage/databases/' + name, (32,32),
                               peak_types, constrain_det_params )


        data_retriever = lambda x, y : data_fetcher( name, x, y ) 

        spec.auto_fit_many_spectra( db, data_retriever,
                                    '../../images/current_fit_images/' + name + '/', (4,4),
                                    group_ranges, peak_locations,
                                    num_peaks_to_detect, primary_peak_detector,
                                    peak_sizes_guesses, peak_width_guesses, det_params_guesses,
                                    peak_mu_offset,
                                    fit_acceptor = fit_acceptor,
                                    params_shuffler = params_shuffler,
                                    rel_plot_bounds = rel_plot_bounds,
                                    logscale = 1, time_estimator = time_estimator,
                                    print_output = 0 )

        mu_path = '../../storage/mu_values/' + name + '_mu_values.bin'

        db.write_mu_values( mu_path )

        db.disconnect()


 
# one_spectrum( (20,19) ) 

# all_spectra( [ 'moved', 'angled', 'centered', 'flat' ] ) #  [ 'moved', 'angled', 'centered', 'flat', 'det3_cent', 'det3_moved' ] )  

# all_spectra( [ 'alpharun20-30', 'alpharun11-19'



all_spectra( [ 'det3_cent', 'det3_moved' ] )


# # plot the histogram without fit yet 
    # if nice_format : 
    #     jplt.plot_histo( ax, xaxis, efront_histo,
    #                      plot_bounds = None, logscale = 1,
    #                      title = "Example Spectrum", xlabel = "Channel",
    #                      ylabel = "Counts" )

    #     labels = [ '$^{240}\mathrm{Pu}$', '$^{238}\mathrm{Pu}$',
    #                    '$^{249}\mathrm{Cf}$' ]

    #     for l in range( 3 ) :

    #         peakpos = our_peaks[ 2*l ]
    #         ax.text( peakpos - 30, 30 + efront_histo[ peakpos ],
    #                  labels[l] ) # , xycoords = 'data' )

    #         # l.draggable()
            
    # else :s
    #     jplt.plot_histo( ax, xaxis, efront_histo,
    #                      plot_bounds = None, logscale = 1,
    #                      title = "", xlabel = "",
    #                      ylabel = "" )
        
# fit_all_spectra( dbmgr.all_dbs )
