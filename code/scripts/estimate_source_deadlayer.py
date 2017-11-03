import libjacob as libj


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




def estimate_source_deadlayer():

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

                        print( str(x) + ' ' + str(j) ) 
                        
                        calibrated_energy_array[i][l][x][j] = meas.meas(
                            m.x * mu_vals_array[i][1][l][j].x + b.x,
                            m.x * mu_vals_array[i][1][l][j].dx )
                        
                        #energy_differences_flat[l].append( calibrated_energy )
                                                
                            
                else:
                    for l in range(2):
                        calibrated_energy_array[i][l][x][j] = meas.nan
                    continue
                
                    # f( mu_vals_array[i][1][l][j] ) 
                    
            
                    
    for j in range( 2 ):
                        
                        # print( i, j ) 

            # moved_mu_tmp = moved_mu[ i, j ]# [ x, : ]
            # centered_mu_tmp = centered_mu[ i, j ] #[ x, : ]
            
            # print( 'moved_mu: ' + str( moved_mu_tmp ) )
            # print( 'centered_mu : ' + str(centered_mu_tmp ) ) 
            
        # mu_differences = moved_mu_tmp - centered_mu_tmp 
            
#            print( 'mu_differences: ' + str( mu_differences ) )
        # print( centered_calibrated_energies )

        # print( '\n' )
        # print (  moved_calibrated_energies )

        energy_differences = centered_calibrated_energies - moved_calibrated_energies

        x = secant_differences.x.flatten()
        xerr = secant_differences.dx.flatten()
        
        y = energy_differences[j].x.flatten()          
        yerr = energy_differences[j].dx.flatten()

        
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
