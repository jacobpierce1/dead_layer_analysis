# contains functions that were useful at one point but not anymore.




def nancount( arr ):
    return np.count_nonzero( ~np.isnan( arr ) )




# this function takes an array of observed channels x and converts to
# energies of the alphas as they were originally emitted, accounting
# for loss of energy in both the detector and source dead layers and
# assuming a linear calibration of the detector.
#
# format: p[ si_deadlayer_efficiency, a, b, detector_dl, pu238_dl,
# pu240_dl, cf249_dl ]
#
# m and b are the parameters of a linear fit y =
# mx + b identified_energy is something extracted from the peaks that
# is expected to be proportional to the energy of the alpha as it
# entered the detector, such as the extrapolated peak position or the
# mu value obtained from the fit. i am going to try using both
# separately.
#
# the x array must have one entry for each source.
def calibrated_energy( p, x, stop_power_interp_functions ):

    E_current = np.empty( x.shape )
    
    # first use linear calibration to get the energy of the alpha as
    # recorded by the detector:
    for i in np.arange(len(E_current)):

        # this is the energy of the alpha entering the detector, assumed to be
        # independent of the distance traveled. such a loss could be incorporated
        # later.

        E_current[i] = (
            _energy_entering_detector( p[0], p[1], p[2], p[4],
                                      stop_power_interp_funcs[-1],
                                      x[i], cosine_matrices[i] ) )
        
    # todo: account for distance traveled through air by the alpha
    #
    #
    
    # # next: account for energy lost in the source dead layer.
    # for i in np.arange(len(E_current)):
    #     E_current[i] += 


    return E_current



# map from E_0 to a spectrum feature, e.g. mu value or peak value.
# p[0] = linear coefficient, p[1] = const offset, p[2] = distance of the
# deadlayre in mm
def energy_to_feature( p, E0, energy_loss_func ):
    return p[0]*( E0 - energy_loss_func( E0, p[2] ) ) + p[1] 




# the procedure is complicated a bit by the fact that we measure both charge from
# the detector (assuming 100% collection efficiency) along with charge from the
# dead layer which by hypothesis does not have perfect collection efficiency.
# this function uses Newton's method to obtain what the energy entering the
# dead layer of the detector must have been given an interpolation function
# of the stopping power, the channel number, and guesses for m, b, and the
# collection efficiency of the dead layer.
#
# in short, solve ax + b = det_dl_efficiency * ( E0 - interp(E0) ) + interp(E0)
# for E0.
def _energy_entering_detector( det_dl_efficiency, a, b, stop_power_interp_func, channel, cosine_matrix ):

    # guess for solution: calibration ignoring the dead layer
    E0_guess = a*x + b

    
    func = lambda E0: ( det_dl_efficiency * (E0 - stop_power_interp_func( E0 ) ) +
                        stop_power_interp_func(E0) - ( a * channel + b ) )

    return scipy.optimize.newton( func, E0_guess )




# correct for the energy lost in the source deadlayer. use Newton's method. 
def _energy_leaving_source( E_entering_det, ):

    E0_guess = E_entering_det

    



def _populate_stop_power_interp_funcs( stop_power_interp_funcs, all_objects_deadlayer_ids ):
    for i in range( len( stop_power_interp_funcs ) ):

        current_file = (
            _dir_stopping_power_data + 'alpha_stopping_power_' +
            _deadlayer_ids[i] + '.txt'
        )
        
        stop_power_interp_funcs[i] = (
            stop.stopping_power_interpolation.from_nist_file( current_file,
                                                              _deadlayer_densities[i]
            ).interp
        )



# def _estimate_deadlayer( cosine_matrices, alpha_energies, p0 ):
    


def _main():

    # this takes a while to populate. use debug=1 option when developing.
    cosine_matrices = geom.get_cosine_matrices( debug=0 )

    # read in grid of the mu values 
    mu_grid_center = anal.get_mu_grid_where_valid( 2, get_db_contents(
        sql_db_manager.centered_db ) )

    # anal.get_mu_grid_where_valid( 
    
    # guess 100 nm depth for all dead layers, both detector and source.
    # calculation is 100 nm * ( 1 m / 10^9 nm ) * 1000 mm / nm 

    # deadlayer_depth_guess = [ 100.0 / 1e9 * 1000 ] * len( _all_objects )
    det_deadlayer_depth_guess = [ 100.0 / 1e9 * 1000 ]
    
    # # guess for the fraction of energy that is collected in the detector deadlayer
    # deadlayer_efficiency_guess = 0.5

    # actual energies of the 5 peaks we are looking at:
    alpha_energies = jutils.flatten_list( anal.peak_energies )

    # construct the stopping power interpolation functions
    _populate_stop_power_interp_funcs( _stop_power_interp_funcs,
                                       _deadlayer_ids )

    # guess for all params: A, B, det deadlayer guess,
    p0 = [ 2.0, 50.0, det_deadlayer_depth_guess ]
    #    estimate_deadlayer( cosine_matrices, alpha_energies, p0 )

    current_row_cosines = cosine_matrices[2,0,0,:]


    calibration = lambda _p, _x: calibrated_energy( _p, _x,
                                                    stop_power_interp_functions )

    
    ret = jmath.jacob_least_squares( current_row_cosines.x, mu_values, mu_values_delta,
                                     p0, calibrated_energy )


    


    
# in principle we should be able to see the expected effect for any source
# because of small angular variation throughout pixels.
# this function makes an estimate using this fact. for each strip we
# look at small angular variations and fit to a line to extract the
# parameter identified as the dead layer distance. then we look at the
# distribution of results obtained using this method.
#
# results: mostly unsuccessful. the plot shows the opposite trend,
# namely increasing mu with cos theta. the angles are just too
# sensitive to the position of the sources and detector, which are
# simply not known that well.  proceeding to another test.

# _DEBUG = 0

# def example_estimate_for_one_source():

#     strip = 16

#     # this takes a while to populate. use debug=1 option when developing.
#     # [i,j,k,l] = [ sourcenum, det or source (0 or 1 ), xcoord, ycoord ] 
#     cosine_matrices = geom.get_cosine_matrices( debug=0 )

#     # read in grid of the mu values as meas class
#     mu_grid_center = dbmgr.rotated.get_mu_grid_where_valid( 3 )

#     if _DEBUG:
#         print_strip_stds( mu_grid_center )
    
#     # get 1 / mu and 1 / cosine for the current strip 
#     cosine_recip_row = meas.abs( 1 / meas.meas.from_list( cosine_matrices[ 3,0,strip,: ] ) )
#     mu_row = mu_grid_center[strip,:]

#     # make a plot and display it.
#     ax = plt.axes()

#     jplt.plot( ax, cosine_recip_row.x, mu_row.x,
#                # xerr = cosine_recip_row.dx,
#                yerr = mu_row.dx,
#                title = r'Sample Plot: Angular Variation of $\mu$ for One Strip',
#                xlabel = r'$ 1 / \cos \theta $' ,
#                ylabel = r'$\mu$ (Alpha Energy Cutoff)' )

#     plt.show()
    
       
#     ret = jmath.jacob_least_squares( current_row_cosines.x, mu_values, mu_values_delta,
#                                      p0, calibrated_energy )


    



# plot_sim_results: optional arg to plot the histogram of the
# simulation results. if not -1, then will plot that number peak
# out of the 6 peaks. so plot_sim_results = 1 will plot peak
# number 1 (0-indexed)

def get_peakvals( db, x, y, plot_sim_results = -1 ):    

    alpha_fitfunc = spec.sum_n_fitfuncs( spec.fitfunc_n_alpha_peaks, 1 )
    peakpos_arr = []
    peakpos_delta_arr = []
    num_iterations = 300

    histo = np.zeros( 5000 )
    data.get_pixel_histo( db, x, y, histo )

    
    # find the 5 peak positions to nearest integer

    peakpos_guesses_arr = [0] * 6
    
    peakpos_guesses_arr, npeaks_found = jmath.get_n_peak_positions( 6, histo )

    if npeaks_found != 6:
        print( 'ERROR: unable to find 6 peaks. ' ) 
        return 0
    
    single_fit_params = db.get_single_peak_fit_parameters( x, y )

    for peaknum in range(6):
                
        # this occurs if there was not a convergent fit
        # recorded in the DB
        
        if np.isnan( single_fit_params.x[ peaknum ][ 0 ] ):
            peakpos_arr.append( np.nan )
            peakpos_delta_arr.append( np.nan )
            continue

        # print( single_fit_params.x[ peaknum ] )

        print( single_fit_params[ peaknum ] )
        
        peakpos_sim_results = jmath.estimate_peakpos( alpha_fitfunc,
                                                      single_fit_params.x[ peaknum ],
                                                      single_fit_params.dx[ peaknum ],
                                                      peakpos_guesses_arr[ peaknum ],
                                                      num_iterations = num_iterations )

        peakpos_arr.append( np.mean( peakpos_sim_results ) )
        peakpos_delta_arr.append( np.std( peakpos_sim_results ) / np.sqrt( num_iterations ) )

        # if PRINT_MAX_PEAKPOS_SIM_RESULTS:
        #     print( 'max( peakpos_sim_results ): ' + str( max( peakpos_sim_results ) ) )
        #     print( 'min( peakpos_sim_results ): ' + str( min( peakpos_sim_results ) ) )

        # print the peakpos detected with no random sampling.
        # if PRINT_BASE_PEAKPOS:
        #     inverted_f = lambda x_: 0 - alpha_fitfunc( all_fit_params['pf_delta_arr'][peaknum], x_ )  
        #     result = scipy.optimize.fmin( inverted_f, peakpos_guesses_arr[peaknum], disp=0 )
        #     print( 'base peakpos: ' + str( result ) )

        # if PRINT_PEAKPOS_ARR:
        #     print( 'peakpos mean: ' + str( peakpos_arr[peaknum] ) )
        #     print( 'peakpos std of mean: ' + str( peakpos_delta_arr[peaknum] ) )

        # histogram the result of the simulation
        if plot_sim_results != -1 and peaknum == plot_sim_results:
            ax = plt.axes()
            ax.hist( peakpos_sim_results, bins=num_iterations // 15 )
            plt.show()
            # return 1

            
    # peakpos_arr is an array containing the estimates positions of the 5 peaks. peakpos_delta_arr
    # is the uncertainty on each estimate. note that the uncertainty is due entirely to uncertainty
    # in fit parameters. peakvals_arr is the value of the function evaluated at each point.
    # peakvals_arr is used to estimate the number of counts expected for each peak.

    peakpos_arr = np.asarray( peakpos_arr )
    peakpos_delta_arr = np.asarray( peakpos_delta_arr )
    peakvals_arr = np.asarray( [ alpha_fitfunc( single_fit_params.x[ peaknum ],
                                                [ peakpos_arr[peaknum] ] )[0]
                                 for peaknum in range(5) ] )

    return ( meas.meas( peakpos_arr, peakpos_delta_arr ), peakvals_arr )










# do the analysis for a particular strip to see
# the relationship between difference of secants and
# difference in mu values.

def preview_secant_differences():

    # the point of this function is that we analyze a particular strip,
    # in this case x == 16 
    x = 16

    moved_peakvals = [0] * 32
    centered_peakvals = [0] * 32 


    # populate the previous two arrays. doing this to reduce duplicate code.
    peakvals_arr = [ moved_peakvals, centered_peakvals ]
    dbs = [ dbman.moved, dbman.centered ]
    for db in dbs:
        
        db.connect()
        for y in range( len( moved_peakvals ) ):
            moved_peakvals[y] = get_peakvals( db, x, y, plot_sim_results = 2 )  
        db.disconnect()


    print( moved_peakvals )
    print( centered_peakvals )
    

        
    f, axarr = plt.subplots( 2 )
    
    cosine_matrices = geom.get_cosine_matrices( debug=0 )


    # construct the sec differences from geometry
    secant_differences = meas.meas.from_list( 1 / cosine_matrices[ 'pu_238_moved' ][ 0, 16, : ]
                                            - 1 / cosine_matrices[ 'pu_238_centered' ][ 0, 16, : ] )

    # now construct mu differences for peaks 2 and 3 at
    # a particular strip, in this case 16.
    
    peaknums = [ 2, 3 ]
    mu_differences = [0] * 2
    
    for i in range( 2 ):

        moved_mu = dbman.get_mu_grid_where_valid( db.moved_db, peaknums[i] )[ 16, : ] 
        centered_mu = dbman.get_mu_grid_where_valid( db.centered_db, peaknums[i] )[ 16, : ] 

        print( 'moved_mu: ' + str( moved_mu ) )
        print( 'centered_mu : ' + str(centered_mu ) ) 
        
        mu_differences[i] = moved_mu - centered_mu 

        print( 'mu_differences[i]: ' + str( mu_differences[i] ) )


        jplt.plot( axarr[i], secant_differences.x,
                   mu_differences[i].x,
                   xerr = secant_differences.dx,
                   yerr = mu_differences[i].dx ) 

        
        
    
    # jplt.plot( ax, secant_differences.x, mu_differences.x,
    #            xerr = secant_differences.dx,
    #            yerr = mu_differences.dx,
    #            xlabel = r'$ \Delta( 1 / \cos( \theta_i ) ) $',
    #            ylabel = r'$ \Delta( \mu_i )$' )

        
    plt.show() 
