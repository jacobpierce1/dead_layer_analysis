# this script reads in dataframes containing: 1. the penetration angle
# of each pixel ( constructed in deadlayer_geometry.py ), 2. the mu
# values of each peak fit ( originally detected in parse_all_data.py
# and read into a matrix in deadlayer_analysis.py ) it then reads in a
# file of stopping power data and interpolates it. using all this data
# we make estimates of both the source dead layer depths and
# the pixel dead layer depth.


import libjacob.jpyplot as jplt
import libjacob.jmeas as meas
import libjacob.jutils as jutils
import libjacob.jmath as jmath
import libjacob.jstats as jstats

# import deadlayer_helpers.stopping_power_interpolation as stop
import deadlayer_helpers.stopping_power_interpolation as stop
import deadlayer_helpers.geometry as geom
import deadlayer_helpers.sql_db_manager as dbmgr
import deadlayer_helpers.analysis as anal
import deadlayer_helpers.data_handler as data

import jspectroscopy as spec


import matplotlib.pyplot as plt
import pandas as pd
import numpy as np
import scipy.optimize


peak_energies = [ [ 5123.68, 5168.17 ], [ 5456.3, 5499.03 ], [ 5759.5, 5813.10,  ] ]
flattened_peak_calibration_energies = jutils.flatten_list(
    [ peak_energies[j] for j in [0,2] ] )


# stopping power for the 5456.3 and 5499.03 keV peaks alphas
# https://physics.nist.gov/cgi-bin/Star/ap_table-t.pl
stopping_power_energies = np.array( [ 5.802E+02 , 5.802E+02 ] )  # MeV cm^2 / g 
density_si = 2.328 # g / cm^2
stopping_power_energies *= density_si * 1000 * 100 / 1e9   # convert to MeV / cm

# modifiable config
_dir_stopping_power_data = '../../data/stopping_power_data/'



# constants, need density of each deadlayer id



# global vars
_sources = geom.sources
_all_objects = geom.all_objects
_deadlayer_ids = [ 'si', 'si', 'si', 'si', 'si', 'si', 'si' ]
_deadlayer_densities = [ 2.328 ] * len(_deadlayer_ids)


# the identity of each dead layer for each unique source. to
# be populated later
_stop_power_interp_funcs = [0] * len(_deadlayer_ids)


#############################################################################




# print std of x and y strips. the interpretation of this is that the
# smaller std values correspond to the same detector properties
# and admit valid comparison. turns out that stripy is the one with a much lower
# standard deviation (about 0.7)

def print_strip_stds( mu_grid ):

    stripx = mu_grid[ :, 15 ]
    stripy = mu_grid[ 15, : ]  # this one is constant roughly.
    
    # stripx = meas.meas.from_list( mu_grid[ :, 15 ] )
    # stripy = meas.meas.from_list( mu_grid[ 15, : ] )

    print( 'stripx: ' + str( stripx ) )
    print( 'stripy: ' + str( stripy ) )
    
    print( 'stripx mean: ' + str( stripx.nanmean() ) )
    print( 'stripx std: ' + str( stripx.nanstd() ) )
    print( 'stripy mean: ' + str( stripy.nanmean() ) )
    print( 'stripy std: ' + str( stripy.nanstd() ) )
    

    




    
    
# DESCRIPTION: with this function we make a few significant
# improvements over example_estimate_for_one_source(). for one, we are
# no longer sensitive to small variations in the detector / source
# position.  as another (which could have been incorporated in the
# previous function but was not) is that we look at the extrapolated
# peak positoins in the alpha spectra instead of the mu values. the
# peak positions are more likely to be correlated with the energy of
# the alpha.
#
# RESULTS:

def estimate_deadlayer_from_energy_differences():
    
    f, axarr = plt.subplots( 2 )


    # lists to store all the energy difference measurements
    # and geometry measurements for 2 middle peaks

    moved_calibrated_energies = meas.meas.empty( (2,32,32 ) )
    centered_calibrated_energies = meas.meas.empty( (2,32,32) )

    calibrated_energy_array = [ centered_calibrated_energies,
                                moved_calibrated_energies ] 
    
    # construct the sec differences from geometry

    cosine_matrices = geom.get_cosine_matrices( debug=0 )
    
    secant_differences = meas.meas.from_array(
        1 / cosine_matrices[ 'pu_238_moved' ][ 0 ]
        - 1 / cosine_matrices[ 'pu_238_centered' ][ 0 ] )

    
    for x in range( 32 ):

        print( x )

        centered_mu = dbmgr.centered.get_mu_for_x_strip( x )
        moved_mu = dbmgr.moved.get_mu_for_x_strip( x )
        
        # secant_differences = meas.meas.from_list(
        #     1 / cosine_matrices[ 'pu_238_moved' ][ 0, x, : ]
        #     - 1 / cosine_matrices[ 'pu_238_centered' ][ 0, x, : ] )

        
        # do least squares fit on each pixel
        
        mu_vals_array = [ centered_mu, moved_mu ]
        
        
        for i in range( 2 ):
            for j in range( 32 ):
                
                if mu_vals_array[i][1][0][j].x == np.nan: 
                    for l in range(2):
                        calibrated_energy_array[i][l][x][j] = meas.nan
                    continue

                                
                # construct arrays of mu values for features 0 and 2, ie the calibration
                # sources. 
                mu_vals = [ mu_vals_array[i][k][l][j].x
                            for k in [0,2]
                            for l in [0,1] ] 
                
                mu_deltas = [ mu_vals_array[i][k][l][j].dx
                              for k in [0,2]
                              for l in [0,1] ] 
                
                if np.count_nonzero( ~np.isnan( mu_vals ) ) < 3:
                    for l in range(2):
                        calibrated_energy_array[i][l][x][j] = meas.nan
                    
                    continue

                # print( mu_vals )
                # print( mu_deltas ) 
                
                linear_fit = jstats.linear_calibration( flattened_peak_calibration_energies,
                                                        mu_vals, mu_deltas, print_results = 0,
                                                        invert = 1 )

                if linear_fit is not None:
                    

                    m, b, f = linear_fit
                        
                    for l in range(2):

                        calibrated_energy_array[i][l][x][j] = meas.meas(
                            m.x * mu_vals_array[i][1][l][j].x + b.x,
                            m.x * mu_vals_array[i][1][l][j].dx )
                        
                        #energy_differences_flat[l].append( calibrated_energy )
                                                
                            
                else:
                    for l in range(2):
                        calibrated_energy_array[i][l][x][j] = meas.nan
                    continue
                
                    # f( mu_vals_array[i][1][l][j] ) 

                    

    # load interpolation
    # si_interpolation = stop.stopping_power_interpolation( 'si', [ 5.40, 5.55 ] ) 
                    
    energy_differences = centered_calibrated_energies - moved_calibrated_energies
    
    x = secant_differences.x.flatten()
    xerr = secant_differences.dx.flatten()

    
    for j in range( 2 ):
                        
                
        y = energy_differences[j].x.flatten()          
        yerr = energy_differences[j].dx.flatten()

        indices = ~ np.isnan( y ) 
        
        energy_differences_mean = meas.meas( y, yerr ).nanmean( option = 'weighted' )
    
        secant_differences_mean = meas.meas( x[ indices ],
                                             xerr[ indices ] ).nanmean( option = 'weighted' )

        
        
        deadlayer = ( energy_differences_mean
                      / ( secant_differences_mean
                          * stopping_power_energies[i] ) )

        print( 'deadlayer: ' + str( deadlayer ) ) 
        
        print( 'energy_differences_mean: ' + str( energy_differences_mean ) )

        print( 'secant_differences_mean: ' + str( secant_differences_mean ) )

                                             

        
        # print( x[ ~np.isnan( y ) ] )
        # print( y[ ~np.isnan( y ) ] )

        jplt.plot( axarr[j], x, y, xerr = xerr, yerr = yerr,
                   ylabel = r'$\Delta E$',
                   xlabel = r'$\Delta \cos \theta$' ) 
                   
        axarr[j].text( 0.1, 0.9, 'Peak %d' % (j,),
                         transform = axarr[j].transAxes,
                         fontsize=12,
                         verticalalignment='top' )

        # jstats.linear_calibration( ax = axarr[j] ) 
        

    plt.show()
    
    return 0
                
                








    









    
# print_strip_stds( dbmgr.moved.get_mu_grids_where_valid()[ 1][0] ) 
    
# example_estimate_for_one_source()
estimate_deadlayer_from_energy_differences()

# preview_secant_differences()

