# This script reads in dataframes containing: 1. the penetration angle
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

from sklearn import linear_model 


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
    



def plot_3d_errorbar( ax, xyz, dxyz ) :

    ax.plot( xyz[0], xyz[1], xyz[2], linestyle="None", marker="o")

    
    # for i in np.arange(0, len(xyz)):
    #     ax.plot( [ xyz[i][0] + dxyz[i][0], xyz[i][0] - dxyz[i][0] ],
    #              [ xyz[i][1] ] * 2,
    #              [ xyz[i][2] ] * 2,
    #              marker="_")

    #     ax.plot( [ xyz[i][0] ] * 2,
    #              [ xyz[i][1] + dxyz[i][1], xyz[i][1] - dxyz[i][1] ],
    #              [ xyz[i][2] ] * 2,
    #              marker="_")
        
    #     ax.plot( [ xyz[i][0] ] * 2,
    #              [ xyz[i][1] ] * 2,
    #              [ xyz[i][2] + dxyz[i][2], xyz[i][2] - dxyz[i][2] ],
    #              marker="_") 
                
        # #configure axes
        # ax.set_xlim3d(0.55, 0.8)
        # ax.set_ylim3d(0.2, 0.5)
        # ax.set_zlim3d(8, 19)
        

    
    
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

def estimate_source_deadlayer():

    fig = plt.figure()
    axarr = [ fig.add_subplot(2, 1, i, projection='3d') for i in [1,2] ]

    # lists to store all the energy difference measurements
    # and geometry measurements for 2 middle peaks

    flat_calibrated_energies = meas.meas.empty( (2,32,32) )
    angled_calibrated_energies = meas.meas.empty( (2,32,32 ) )

    calibrated_energy_array = [ flat_calibrated_energies,
                                angled_calibrated_energies ] 
    
    # construct the sec differences from geometry

    all_cosine_matrices = geom.get_cosine_matrices( compute_source_costheta = 1 )

   
    secant_matrices = [ 1 / meas.meas.from_array( all_cosine_matrices[i][0] ) 
                        for i in [ 'pu_238_flat', 'pu_238_angled' ] ]

    source_secant_matrices =  [ 1 / meas.meas.from_array( all_cosine_matrices[i][1] ) 
                                for i in [ 'pu_238_flat', 'pu_238_angled' ] ]
        
    for x in range( 32 ):

        print( x )

        flat_mu = dbmgr.flat.get_mu_for_x_strip( x )
        angled_mu = dbmgr.angled.get_mu_for_x_strip( x )

        
        # do least squares fit on each pixel
        
        mu_vals_array = [ flat_mu, angled_mu ]
        
        
        for i in range( 2 ):
            for y in range( 32 ):
                
                if meas.isnan( mu_vals_array[i][1][0][y] ) : 
                    for l in range(2):
                        calibrated_energy_array[i][l][x][y] = meas.nan
                    continue

                                
                # construct arrays of mu values for features 0 and 2, ie the calibration
                # sources. 
                mu_vals = [ mu_vals_array[i][k][l][y].x
                            for k in [0,2]
                            for l in [0,1] ] 
                
                mu_deltas = [ mu_vals_array[i][k][l][y].dx
                              for k in [0,2]
                              for l in [0,1] ] 
                
                if np.count_nonzero( ~np.isnan( mu_vals ) ) < 3:
                    for l in range(2):
                        calibrated_energy_array[i][l][x][y] = meas.nan
                    continue

                                              
                linear_fit = jstats.linear_calibration( flattened_peak_calibration_energies,
                                                        mu_vals, mu_deltas, print_results = 0,
                                                        invert = 1 )

                if linear_fit is not None:                    

                    # print('success') 
                    m, b, f = linear_fit
                        
                    for l in range(2):

                        # calibrated_energy_array[i][l][x][y] = m * mu_vals_array[i][1][l][y] + b
                        energy = meas.meas(
                            m.x * mu_vals_array[i][1][l][y].x + b.x,
                            m.x * mu_vals_array[i][1][l][y].dx )
                            
                        calibrated_energy_array[i][l][x][y] = energy
                        
                else:
                    for l in range(2):
                        calibrated_energy_array[i][l][x][y] = meas.nan
                    continue
                
                    

    axarr[0].set_title( 'Source and Detector: Combined Dead Layer Estimate' ) 

    
    # these contain 1D jmeas arrays: index 0 for flat, index 1
    # for moved.
    x = meas.meas.from_list( [ secant_matrices[i].flatten() for i in [0,1] ] ) 
    y = meas.meas.from_list( [ source_secant_matrices[i].flatten() for i in [0,1] ] )

    xflat = x.flatten()
    yflat = y.flatten()

    colors = [ 'r', 'b' ]


    for j in range(2):

        z = meas.meas.from_list( [ calibrated_energy_array[i][j].flatten() for i in [0,1] ] )
        zflat = z.flatten() 
        
        # plot it 
        for l in range(2):        
            plot_3d_errorbar( axarr[j], [ x.x[l], y.x[l], z.x[l] ],
                              [ x.dx[l], y.dx[l], z.dx[l] ] )

            
        clf = linear_model.LinearRegression()

        mask = ~ ( meas.isnan( xflat ) | meas.isnan( yflat ) | meas.isnan( zflat ) )

        clf.fit( xflat.x[mask].reshape(-1,1), zflat.x[mask].reshape(-1,1) )          
        print( clf.coef_ )
        print( 'dl estimate: ' + str( np.abs( clf.coef_[0] ) / stopping_power_energies[j] )  )


        train2 = np.array( [ xflat.x[mask].reshape(-1,1), yflat.x[mask] ] )
        
        clf2 = linear_model.LinearRegression()
        clf2.fit( train2, zflat.x[mask].reshape(-1,1) )          
        print( clf2.coef_ )
        print( 'dl estimate: ' + str( np.abs( clf2.coef_[0] ) / stopping_power_energies[j] )  )

        
    plt.show()
        
    return 0
    
        
    # for j in range(2) :

   
                               
    #     jplt.plot( axarr[j], x1, y1, xerr = x1err, yerr = y2err,
    #                ylabel = r'$E$',
    #                xlabel = r'$\sec \theta$',
    #                color = colors[0],
    #                leglabel = 'Flat' )

    #     jplt.plot( axarr[j], x2, y2, xerr = x2err, yerr = y2err,
    #                color = colors[1], leglabel = 'Angled' )
                   
        
    #     axarr[j].text( 0.1, 0.9, 'Peak %d' % (j,),
    #                    transform = axarr[j].transAxes,
    #                    fontsize=12,
    #                    verticalalignment='top' )

        
        
        # x = [ x1, x2 ]
        # xerr = [x1err, x2err ] 
        # y = [ y1, y2 ]
        # yerr = [ y1err, y2err ]

        # labels = [ 'Flat Fit', 'Angled Fit' ]
        # fitcolors = [ 'y', 'k' ]
        
        
        # for k in range(2):
            
        #     ret = jstats.linear_calibration( x[k], y[k], dy = yerr[k],
        #                                      ax = axarr[j], color = fitcolors[k],
        #                                      leglabel = labels[k],
        #                                      linestyle = '--' )
            
        #     if ret is not None:
                
        #         # deadlayer = ret[0] / stopping_power_energies[i]

        #         slopes[k] = ret[0]
                
        #         print( ret[2].fit_report() )  
                
        #         # print( 'deadlayer: ' + str( deadlayer ) )
                
        #         # axarr[j].text( 0.1, 0.2, r'Effective DL = $ %d \pm %d $ nm ' % ( abs( deadlayer.x ) , deadlayer.dx ),
        #         #                transform = axarr[j].transAxes,
        #         #                fontsize=12,
        #         #                verticalalignment='top' )
        #     print( '\n\n' )
            
        # print( slopes )

        # jplt.add_legend( axarr[j], 1 ) 
                
            
        
    # plt.show()
        
    # return 0

                








    









    
# print_strip_stds( dbmgr.moved.get_mu_grids_where_valid()[ 1][0] ) 
    
# example_estimate_for_one_source()
estimate_source_deadlayer()

# preview_secant_differences()

