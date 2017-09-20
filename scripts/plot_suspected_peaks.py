# the purpose of this script is to read in data for the energies of known alphas
# in the spectra of the 3 sources used in the calibration and plot it, using calibration 


## includes
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import scipy.optimize

# my lib 
import libjacob.jacob_math  as jmath
import libjacob.jacob_pyplot as jplt
import libjacob.jacob_utils as jutils 
import libjacob.jacob_stats as jstats
import deadlayer_analysis as dla
import deadlayer_functions as dlf

import sql_db_manager
import deadlayer_analysis
import sqlite3


# configurable test file for example
PRINT_FIT_PARAMS = 1
PRINT_PEAKPOS_GUESS = 1
PRINT_PEAKPOS_ARR = 1
PRINT_MAX_PEAKPOS_SIM_RESULTS = 1
PRINT_BASE_PEAKPOS = 1
DISPLAY_FIT_0 = 0
DISPLAY_PEAKPOS_HIST = -1
ADJUST_PEAKPOS_GUESSES_ARR = 0



test_dir = '../data/extracted_ttree_data/'
test_file = test_dir + '/deadlayerdet3cent/deadlayerdet3cent_30_31.bin'
# test_file = test_dir + '/deadlayerdet3rt/deadlayerdet3rt_31_31.bin'



# # configurable test file for example()
# test_dir = '../data/extracted_ttree_data/'
# test_file = test_dir + '/deadlayerdet3rt/deadlayerdet3cent_31_31.bin'


def plot_suspected_peaks():

    # clear plot and create array of subplots
    plt.clf()
    f, axarr = plt.subplots( 2 ) # , sharex='all' )


    
    # need a db / file where every fit converged so that we can perform the calibration,
    # but the suspected features are visible.

    db_name = sql_db_manager.centered_db
    df = dla.read_db_into_df( db_name )
    coords = (30,31)

    # db_name = sql_db_manager.rotated_db
    # df = dla.read_db_into_df( db_name )
    # coords = (31,31)

    # # read in mu and mu uncertainty values
    # mu_values = dla.get_mu_values( df, coords )
    # mu_delta_values = dla.get_mu_delta_values( df, coords )

    
    # read in all values and uncertainties for the 5 peaks
    all_fit_params = dla.get_all_single_peak_fit_parameters( df, coords )
    if PRINT_FIT_PARAMS:
        print 'all_fit_params: '
        print( np.array( all_fit_params['pf_arr']) )
        print ''
        print( np.array( all_fit_params['pf_delta_arr']) )

        #         print 'all_fit_params: ' + str( all_fit_params )
        # print_df( all_fit_params )
    
    # check if nan 
    if np.any( np.isnan( all_fit_params['pf_arr'] ) ):
        print 'ERROR: one of the fits did not converge, use different pixel.'
        return 0


    # read in the histogram. this is needed for the peak position estimate.
    x = np.arange(5000)
    histo = np.zeros( x.size )
    dlf.construct_histo_array( test_file, histo )

    
    # find the 5 peak positions to nearest integer 
    peakpos_guesses_arr = [0]*5
    if not dlf.get_5_peak_positions( peakpos_guesses_arr, histo ):
        print 'ERROR: unable to find 5 peaks. '
        return 0

    
    if ADJUST_PEAKPOS_GUESSES_ARR:
        adjust_peakpos_guesses_arr( peakpos_guesses_arr )

    if PRINT_PEAKPOS_GUESS:
        print 'peakpos_guesses_arr: ' + str( peakpos_guesses_arr )

        
    # now use those as the guesses to estimate each of the peak positions.
    # a monte carlo simulation is used.
    alpha_fitfunc = dlf.n_fitfuncs_abstract( dlf.fitfunc_n_alpha_peaks, 1 )
    peakpos_arr = []
    peakpos_delta_arr = []
    num_iterations = 300

    for peaknum in range(5):
        print ''
        peakpos_sim_results = jmath.estimate_peakpos( alpha_fitfunc, all_fit_params['pf_arr'][peaknum],
                                all_fit_params['pf_delta_arr'][peaknum], peakpos_guesses_arr[peaknum],
                                              num_iterations = num_iterations )

        peakpos_arr.append( np.mean( peakpos_sim_results ) )
        peakpos_delta_arr.append( np.std( peakpos_sim_results ) / np.sqrt( num_iterations ) )

        if PRINT_MAX_PEAKPOS_SIM_RESULTS:
            print 'max( peakpos_sim_results ): ' + str( max( peakpos_sim_results ) )
            print 'min( peakpos_sim_results ): ' + str( min( peakpos_sim_results ) )

        # print the peakpos detected with no random sampling.
        if PRINT_BASE_PEAKPOS:
            inverted_f = lambda x_: 0 - alpha_fitfunc( all_fit_params['pf_delta_arr'][peaknum], x_ )  
            result = scipy.optimize.fmin( inverted_f, peakpos_guesses_arr[peaknum], disp=0 )
            print 'base peakpos: ' + str( result )

        if PRINT_PEAKPOS_ARR:
            print 'peakpos mean: ' + str( peakpos_arr[peaknum] )
            print 'peakpos std of mean: ' + str( peakpos_delta_arr[peaknum] )

        # histogram the result of the simulation
        if DISPLAY_PEAKPOS_HIST != -1 and peaknum == DISPLAY_PEAKPOS_HIST:
            axarr[1].hist( peakpos_sim_results, bins=num_iterations // 15 )
            plt.show()
            return 1


        
    # peakpos_arr is an array containing the estimates positions of the 5 peaks. peakpos_delta_arr
    # is the uncertainty on each estimate. note that the uncertainty is due entirely to uncertainty
    # in fit parameters. peakvals_arr is the value of the function evaluated at each point.
    # peakvals_arr is used to estimate the number of counts expected for each peak.
    peakpos_arr = np.asarray( peakpos_arr )
    peakpos_delta_arr = np.asarray( peakpos_delta_arr )
    peakvals_arr = np.asarray( [ alpha_fitfunc( all_fit_params['pf_arr'][peaknum],
                                            [ peakpos_arr[peaknum] ] )[0] for peaknum in range(5) ] )
    if PRINT_PEAKPOS_ARR:
        print 'peakvals_arr: ' + str(peakvals_arr)
    

    
    
                                              
    # this is a constant array of the main peak energies, we have to flatten though.
    # used for the calibration.
    print ''
    energies = np.asarray( jutils.flatten_list( deadlayer_analysis.peak_energies ) )
    
    # perform linear calibrati# on, this adds a plot
    p0 = [ 0.5, -118.0 ]
    # # cal = jstats.linear_calibration( energies, mu_values, mu_delta_values,
    #                                  # p0, print_fit_data=0, ax=axarr[1], invert=1 )
    cal = jstats.linear_calibration( energies, peakpos_arr, peakpos_delta_arr,
                                    p0, print_fit_data=1, ax=axarr[1], invert=0, scatter=1 )

    # try a quadratic calibration
    # p0 = [0.01, 0.5, -120.0 ]
    # cal = jstats.quadratic_calibration( energies, peakpos_arr, peakpos_delta_arr,
                                     # p0, print_fit_data=1, ax=axarr[1], invert=0, scatter=1 )
    
    if cal is None:
        print 'ERROR: linear calibration failed. '
        return 0
    axarr[1].set_title('Calibration', fontsize=22)
    axarr[1].set_xlabel( 'Energy of Alpha', fontsize=18 )
    axarr[1].set_ylabel( 'Peak Value of Alpha', fontsize=18 )
        
    # add histogram of detector data to plot.
    jplt.plot_histo( axarr[0], x, histo, plot_bounds=[peakpos_arr[0]-100, peakpos_arr[4]+300], #=map( lambda x: cal.f(x), [2700,3300] ),
                     ylabel='Counts', title='Suspected Alpha Peaks', logscale=1 ) # , xlabel='Calibrated Energy (keV)' )

    
    # now read in all peak values. all_spectra is a series containing a dataframe for each source.
    # each dataframe has a col of energy and associated intensity. all values are then added to plot
    all_spectra = dla.get_all_alpha_spectra()
    
    # intensities plotted relative to these amplitudes, the 3 highest peaks which are numbered 1, 3, 4
    # starting from 0
    max_peak_counts = pd.Series( peakvals_arr[ [1,3,4] ], index=['pu240', 'pu238', 'cf249'] )
    colors = pd.Series( ['r', 'g', 'y'],  index=['pu240', 'pu238', 'cf249'] )
    
    for source, df in all_spectra.iteritems():
        for index, data in df.iterrows(): # index is a dummy, just the col number.
            predicted_counts = data['intensity'] * max_peak_counts[source] / dla.max_peak_intensities[source]
            predicted_chan = cal.f( data['energy'] )
            axarr[0].plot( [ predicted_chan, predicted_chan ], [ 0, predicted_counts ],
                           color=colors[ source ], label=source, linewidth=3 )
            jplt.add_legend( axarr[0], 1 )
    

            
    # this was used for debugging. for some reason the extrapolated peak positions were quite far off
    # for the first peak. i hypothesized that the reason is that the sum of the 2 peaks was a good fit
    # but not the individual fits. this test refuted my hypothesis: the individual fits were clearly
    # very good.
    if DISPLAY_FIT_0:

        # get fit data, could be reconstrcuted from the data that was already read but this is
        # faster to type.
        with sqlite3.connect( db_name ) as sql_conn:
            result = sql_db_manager.read_data_from_db( sql_conn, coords, 0 )
            successful_fit = result[0]
            if successful_fit:
                pf_fit0 = result[3]
                pf_delta_fit0 = result[4]
                fit_bounds = result[6]

        # add plot of the 2 peaks 
        two_peaks_fitfunc = dlf.n_fitfuncs_abstract( dlf.fitfunc_n_alpha_peaks, 2 )
        jplt.add_fit_to_plot( axarr[0], x, fit_bounds, pf_fit0, pf_delta_fit0, two_peaks_fitfunc )

        # add single peak plots
        jplt.add_fit_to_plot( axarr[0], x, fit_bounds, all_fit_params['pf_arr'][0],
                              all_fit_params['pf_delta_arr'][0], alpha_fitfunc, color='--g' )
        jplt.add_fit_to_plot( axarr[0], x, fit_bounds, all_fit_params['pf_arr'][1],
                              all_fit_params['pf_delta_arr'][1], alpha_fitfunc, color='--y' )


        
    plt.show()      
    
    return 1

    

    
def print_df( df ):
    with pd.option_context( 'display.max_rows', None, 'display.max_columns', 3 ):
        print( df )




def adjust_peakpos_guesses_arr( peakpos_guesses_arr ):
    for i in range(len(peakpos_guesses_arr)):
        peakpos_guesses_arr[i] -= 1



        
plot_suspected_peaks()
