# the purpose of this script is to read in data for the energies of known alphas
# in the spectra of the 3 sources used in the calibration and plot it, using calibration 


## includes 
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

# my lib 
import libjacob.jacob_math as jmath
import libjacob.jacob_pyplot as jplt
import libjacob.jacob_utils as jutils 
import libjacob.jacob_stats as jstats
import deadlayer_analysis as dla
import deadlayer_functions as dlf

import sql_db_manager
import deadlayer_analysis



# configurable test file for example
test_dir = '../data/extracted_ttree_data/'
test_file = test_dir + '/deadlayerdet3cent/deadlayerdet3cent_30_31.bin'



# # configurable test file for example()
# test_dir = '../data/extracted_ttree_data/'
# test_file = test_dir + '/deadlayerdet3rt/deadlayerdet3cent_31_31.bin'


def plot_suspected_peaks():

    # need a db / file where every fit converged so that we can perform the calibration,
    # but the suspected features are visible.
    db_name = sql_db_manager.centered_db
    df = dla.read_db_into_df( db_name )

    coords = (30,31) 
    mu_values = dla.get_mu_values( df, coords )
    mu_delta_values = dla.get_mu_delta_values( df, coords )

    # check if nan 
    if any( np.isnan( mu_values ) ):
        print 'ERROR: one of the fits did not converge, use different pixel.'
        return 0


    plt.clf()
    f, axarr = plt.subplots( 2, sharex='all' )
    
    p0 = [ 0.5, -118.0 ]


    # this is a constant array of the peak energies, we have to flattne though. 
    energies = np.asarray( jutils.flatten_list( deadlayer_analysis.peak_energies ) )
    
    # perform linear calibration, this adds a plot
    cal = jstats.linear_calibration( energies, mu_values, mu_delta_values, p0, print_fit_data=0, ax=axarr[1], invert=1 )
    if cal is None:
        print 'ERROR: linear calibration failed. '
        return 0
        
    # read in the histogram and add to the left plot
    x = np.arange(5000)
    histo = np.zeros( x.size )
    dlf.construct_histo_array( test_file, histo )

    jplt.plot_histo( axarr[0], cal.f(x), histo, plot_bounds=None, #=map( lambda x: cal.f(x), [2700,3300] ),
                     ylabel='Counts', title='Suspected Alpha Peaks' ) # , xlabel='Calibrated Energy (keV)' )


    # find the 5 peak positions to nearest integer 
    peakpos_guesses_arr = [0]*5
    if not dlf.get_5_peak_positions( peakpos_guesses_arr, histo ):
        print 'ERROR: unable to find 5 peaks. '
        return 0

    
    
    
    # now read in all peak values. all_spectra is a series containing a dataframe for each source.
    # each dataframe has a col of energy and associated intensity. all values are then added to plot
    all_spectra = dla.get_all_alpha_spectra()
    for source, df in all_spectra.iteritems():
        for index, data in df.iterrows():
            pass
    
    # customize 
    plt.show()      
    
    return 1

    
    
    
    



    
    
plot_suspected_peaks()
