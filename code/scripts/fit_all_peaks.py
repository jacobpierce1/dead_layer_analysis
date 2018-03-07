# -*- coding: utf-8 -*-

import jspectroscopy as spec
import numpy as np
import deadlayer_helpers.data_handler as data
import matplotlib.pyplot as plt 



# inputs = spectrum_fitter_inputs( (32,32), data_fetcher )

# num_groups = 3

rel_plot_bounds = [ -100, 120 ] 

group_ranges = [ [ -65, 15 ], [-75, 16], [-70,85] ] 

# peak_structures = [ [2], [3], [5] ] 

peak_locations = [ [ -50, -20, 0 ], [ -57, -20, 0 ], [ -50, -30, 0, 22, 46, 64 ] ]

peak_types = [ ['a'] * len(x) for x in peak_locations ]

peak_sizes_guesses = [ [ 2000, 30000.0, 90000.0 ],
                     [ 3000, 50000.0, 100000. ],
                     [ 5000.0, 10000.0, 100000.0, 100.0, 500.0, 500.0 ] ]

det_params_guesses = [ { 'a' : [ 5.0, 0.97, 35.0, 1.5 ] } ] * 3 

peak_mu_offset = 9

num_peaks_to_detect = 6

# primary_peak_ids = None

# self.peak_structures = None 



                   

# create and populate the histogram array 
def data_fetcher( name, indices ) :

    xaxis = np.arange( 5000 )
        
    x, y = indices 

    infile = ( '../../data/extracted_ttree_data/'
               + name + '/%s_%d_%d.bin' % ( name, x,y ) )

    print( infile ) 
    
    efront_histo = np.zeros( xaxis.size )

    if not data.construct_histo_array( infile, efront_histo ) :
        print( 'error: couldnt open file' ) 

        # msg = "ERROR: unable to open file: " + infile 

        # print( msg )
        
        # if log_enabled:
        #     log_message( logfile, msg )

        # sys.exit(0)

    # print( max( efront_histo )  ) 
        
    return ( xaxis, efront_histo ) 




# a function that is guarenteed to detect the reference peak
# of each group given a list of the peaks detected in the spectrum

def primary_peak_detector( peaks_detected, histo ) :

    # shifts = np.zeros( 5 ) 
    
    ret = np.empty( 3 )

    ret[0] = peaks_detected[1]

    if histo[ peaks_detected[3] ] / histo[ peaks_detected[2] ] > 5 :
        ret[1] = peaks_detected[4]
        ret[2] = peaks_detected[5]

    else :
        ret[1] = peaks_detected[3]
        
        if histo[ peaks_detected[4] ] < histo[ peaks_detected[5] ] :
            ret[2] = peaks_detected[5]

        else :
            ret[2] = peaks_detected[4]
            
    return ret 





def params_shuffler() :
    return 1 



    

def fit_acceptor( x, y, dy, spec_fitter_result ) :

    return 1



x, y = data_fetcher( 'angled', (20,16) ) 


dy = np.sqrt( y )
dy[ ( dy == 0 ) ] = 1 


plt.figure(figsize=(10,12))
ax = plt.axes() 



test_db = spec.spectrum_db( '../../storage/databases/test_db', peak_types, (32,32) ) 

    

# spec.auto_fit_spectrum( x, y, dy,
#                         group_ranges, peak_locations,
#                         num_peaks_to_detect, primary_peak_detector,
#                         peak_sizes_guesses, det_params_guesses, peak_mu_offset,
#                         fit_acceptor = fit_acceptor,
#                         params_shuffler = params_shuffler,
#                         ax = ax,
#                         rel_plot_bounds = rel_plot_bounds )


spec.auto_fit_many_spectra( test_db, '../../images/current_fit_images/det1_angled_test',
                            group_ranges, peak_locations,
                            num_peaks_to_detect, primary_peak_detector,
                            peak_sizes_guesses, det_params_guesses, peak_mu_offset,
                            fit_acceptor = fit_acceptor,
                            params_shuffler = params_shuffler,
                            ax = ax,
                            rel_plot_bounds = rel_plot_bounds )



plt.show()




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
            
    # else :
    #     jplt.plot_histo( ax, xaxis, efront_histo,
    #                      plot_bounds = None, logscale = 1,
    #                      title = "", xlabel = "",
    #                      ylabel = "" )
        
# fit_all_spectra( dbmgr.all_dbs )
