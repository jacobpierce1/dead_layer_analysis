# fit all peaks for data with gadolinium source



# -*- coding: utf-8 -*-

import jspectroscopy as spec
import numpy as np
import deadlayer_helpers.data_handler as data
import matplotlib.pyplot as plt 
import libjacob.jutils as jutils

debug = 0



# inputs = spectrum_fitter_inputs( (32,32), data_fetcher )

# num_groups = 3

rel_plot_bounds = [ -100, 120 ] 


# group_ranges = [ [ -65, 15 ], [-75, 16], [-70,85] ] 


# peak_locations = [ [ -50, -20, 0 ], [ -57, -20, 0 ], [ -50, -30, 0, 22, 46, 64 ] ]

# peak_sizes_guesses = [ [ 2000., 30000.0, 90000.0 ],
#                        [ 3000., 50000.0, 100000. ],
#                        [ 5000.0, 10000.0, 100000.0, 100.0, 500.0, 500.0 ] ]

# peak_locations = [ [ -42, -20, 0 ], [ -57, -20, 0 ], [ -30, 0, 22, 46, 64 ] ]

# peak_sizes_guesses = [ [ 2000., 30000.0, 90000.0 ],
#                        [ 3000., 50000.0, 100000. ],
#                        [ 10000.0, 100000.0, 100.0, 500.0, 500.0 ] ]

# group_ranges = [ [ -65, 15 ], [-75, 16], [-40, 85] ]




group_ranges = [ [-35,12], [-40,12] ]
peak_locations = [ [0], [-21, 0] ]
peak_sizes_guesses = [ [5000], [1500, 5000] ]

peak_types = [ ['a'] * len(x) for x in peak_locations ]

peak_width_guesses = [ [5.0] * len(x) for x in peak_locations ]

det_params_guesses = [ { 'a' : [ 0.9, 40.0, 2.4 ] } ] * 2

peak_mu_offset = 7.0 # 9 

num_peaks_to_detect = 2


coords = ( 16, 17 ) 




# create and populate the histogram array 
def data_fetcher( name, detnum, x, y ) :

    xaxis = np.arange( 5000 )

    if detnum is not None : 
        infile = ( '../../data/extracted_ttree_data/'
                   + name + '/%s_%d_%d_%d.bin' % ( name, detnum, x, y ) )

    else :
        infile = ( '../../data/extracted_ttree_data/'
                   + name + '/%s_%d_%d.bin' % ( name, x, y ) )
    
    # print( infile ) 
    
    efront_histo = np.zeros( xaxis.size )

    if not data.construct_histo_array( infile, efront_histo, 0 ) :
        print( 'error: couldnt open file' ) 

    dy = np.sqrt( efront_histo )
    dy[ dy==0 ] = 1 
        
    return ( xaxis, efront_histo, dy ) 






# a function that is guarenteed to detect the reference peak
# of each group given a list of the peaks detected in the spectrum

def primary_peak_detector( peaks_detected, histo ) :
    
    if debug : 
        print( peaks_detected ) 

    if len( peaks_detected ) < 2 :
        return None

    # indices = np.empty(2, dtype = int )

    # ret = np.empty( 3, dtype = float ) 
    
    # # primary peaks 0 and 2 are the largest peaks from the left and right, respectively.

    # # we are counting on detection of the first peak. without that
    # # all is lost
    
    # if abs( peaks_detected[1] - peaks_detected[0] ) > 45 :
    #     return None 
        
    # indices[0] = 1
    
    # peak0 = peaks_detected[1]
    # ret[0] = peak0
    
    # next_peak_offsets = [ 174, 320 ] 
    
    # for i in range(2) :
    #     ret[ i+1 ] = ( peak0 + next_peak_offsets[i] - 15
    #                    + np.argmax( histo[ peak0 + next_peak_offsets[i] - 15 :
    #                                        peak0 + next_peak_offsets[i] + 15 ] ) )
    
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




def one_spectrum( db_name, detnum, coords ) :
    
    x, y, dy = data_fetcher( db_name, detnum, *coords ) 

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

    time_estimator = jutils.time_estimator( len(db_names) * 4 * 32 * 32, 20 )

    for name in db_names :

        print( name ) 

        for detnum in [1, 2, 3, 4] : 

            print( detnum ) 
            
            db = spec.spectrum_db( '../../storage/databases/%s_%d' % ( name, detnum ), (32,32),
                                   peak_types, constrain_det_params )


            data_retriever = lambda x, y : data_fetcher( name, detnum, x, y ) 

            spec.auto_fit_many_spectra( db, data_retriever,
                                        '../../images/current_fit_images/' + name + '/', (4,4),
                                        group_ranges, peak_locations,
                                        num_peaks_to_detect, primary_peak_detector,
                                        peak_sizes_guesses, peak_width_guesses,
                                        det_params_guesses,
                                        peak_mu_offset,
                                        fit_acceptor = fit_acceptor,
                                        params_shuffler = params_shuffler,
                                        rel_plot_bounds = rel_plot_bounds,
                                        logscale = 1, time_estimator = time_estimator,
                                        print_output = 0 )

            mu_path = '../../storage/mu_values/%s_%d_mu_values.bin' % ( name, detnum ) 

            db.write_mu_values( mu_path )
            db.disconnect()




def compute_success_rates( db_names ) :

    constrain_det_params = { 'a' : 1 }

    
    for name in db_names :

        print( name ) 

        for detnum in [1, 2, 3, 4] : 

            print( detnum )
            
            mu_path = '../../storage/mu_values/%s_%d_mu_values.bin' % ( name, detnum ) 
            
            db = spec.spectrum_db( '../../storage/databases/%s_%d' % ( name, detnum ), (32,32),
                                   peak_types, constrain_det_params )
            
            # compute yield
            values = db.read_values( mu_path )

            total = 0
            success = 0 
            
            for i in range( len( peak_locations ) ) :
                total += 32 * 32
                success += np.count_nonzero( np.isnan( values[i][0].x ) )
            
            success_rate = success / total 
            print( 'success rate: %.2f' % success_rate ) 



# one_spectrum( 'alpharun20-30', 2, coords  ) 


all_spectra( [ 'alpharun20-30', 'alpharun11-19' ] ) 

compute_success_rates( [ 'alpharun20-30', 'alpharun11-19' ] )





# one_spectrum( (20,19) ) 

# all_spectra( [ 'moved', 'angled', 'centered', 'flat' ] ) #  [ 'moved', 'angled', 'centered', 'flat', 'det3_cent', 'det3_moved' ] )  

# all_spectra( [ 'alpharun20-30', 'alpharun11-19'



# all_spectra( [ 'det3_cent', 'det3_moved' ] )


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
