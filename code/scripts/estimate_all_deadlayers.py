# This script reads in dataframes containing: 1. the penetration angle
# of each pixel ( constructed in deadlayer_geometry.py ), 2. the mu
# values of each peak fit ( originally detected in parse_all_data.py
# and read into a matrix in deadlayer_analysis.py ) it then reads in a
# file of stopping power data and interpolates it. using all this data
# we make estimates of both the source dead layer depths and
# the pixel dead layer depth.

import sys 
import time

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
import lmfit


peak_energies = np.array( [ [ 5123.68, 5168.17 ],
                            [ 5456.3, 5499.03 ],
                            [ 5759.5, 5813.10,  ] ] )

flattened_peak_calibration_energies = peak_energies.flatten()



si_stopping_powers = np.array( [ [ 6.080E+02, 6.046E+02 ],
                                 [ 5.832E+02, 5.802E+02 ],
                                 [ 5.626E+02, 5.591E+02 ] ] )

density_si = 2.328 # g / cm^2

si_stopping_powers *= density_si * 1000 * 100 / 1e9   # convert to keV / nm from mev / cm



# stopping power for the 5456.3 and 5499.03 keV peaks alphas
# https://physics.nist.gov/cgi-bin/Star/ap_table-t.pl
# stopping_power_energies = np.array( [ 5.802E+02 , 5.802E+02 ] )  # MeV cm^2 / g 




# constants, need density of each deadlayer id



# # global vars
# _sources = geom.sources
# _all_objects = geom.all_objects
# _deadlayer_ids = [ 'si', 'si', 'si', 'si', 'si', 'si', 'si' ]
# _deadlayer_densities = [ 2.328 ] * len(_deadlayer_ids)

# # the identity of each dead layer for each unique source. to
# # be populated later
# _stop_power_interp_funcs = [0] * len(_deadlayer_ids)


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
        

    




def get_flattened_secant_angles() :
    
    # construct the sec differences from geometry

    all_cosine_matrices = geom.get_cosine_matrices( compute_source_costheta = 1 )

    keys = all_cosine_matrices.keys()

    print( keys ) 
    
    flattened_det_secant_theta = [0] * len( keys )
    flattened_source_secant_theta = [0] * len( keys ) 
   
    for i in range( len( keys ) ):
        flattened_det_secant_theta[ i ] = ( 1 / all_cosine_matrices[ keys[i] ][0] ).flatten()
        flattened_source_secant_theta[ i ] = ( 1 / all_cosine_matrices[ keys[i] ][i] ).flatten()
        
    return flattened_det_secant_theta, flattened_source_secant_theta, keys


    




def get_flattened_e_from_mu_calibration( db ):

    db_calibrated_energies = meas.meas.empty( (2,32,32) )

    mu_vals_array = db.get_all_mu_grids( 1 )
        
    for x in range( 32 ):
        
        for y in range( 32 ):
            
            if meas.isnan( mu_vals_array[1][0][x,y] ) : 
                for l in range(2):
                    db_calibrated_energies[l][x][y] = meas.nan
                continue
                    
                    
            # construct arrays of mu values for features 0 and 2, ie the calibration
            # sources. 
            mu_calibration_vals = [ mu_vals_array[k][l][x,y].x
                                    for k in [0,2]
                                    for l in [0,1] ] 
            
            mu_calibration_deltas = [ mu_vals_array[k][l][x,y].dx
                                          for k in [0,2]
                                          for l in [0,1] ] 
                
            if np.count_nonzero( ~np.isnan( mu_calibration_vals ) ) < 3:
                for l in range(2):
                    db_calibrated_energies[l][x][y] = meas.nan
                continue
                
                
            linear_fit = jstats.linear_calibration( flattened_peak_calibration_energies,
                                                    mu_calibration_vals,
                                                    mu_calibration_deltas,
                                                    print_results = 0,
                                                    invert = 1 )
            
            if linear_fit is not None:                    
                
                # print('success') 
                m, b, f = linear_fit
                
                for l in range(2):
                    
                    # calibrated_energy_array[i][l][x][y] = m * mu_vals_array[i][1][l][y] + b
                    energy = meas.meas(
                        m.x * mu_vals_array[1][l][x,y].x + b.x,
                        m.x * mu_vals_array[1][l][x,y].dx )
                    
                    calibrated_energy_array[l][x][y] = energy
                    
                else:
                    for l in range(2):
                        calibrated_energy_array[l][x][y] = meas.nan
                    continue
                    

    # flatten everything and return.
    energies = [ calibrated_energy_array[i][j].flatten() for j in [0,1] ] 
    return energies

    
                











def estimate_source_deadlayer_using_mu_calibration_sklearn() :

    ret = get_e_and_angles( dbmgr.flat, dbmgr.angled ) 

    fig = plt.figure()
    axarr = [ fig.add_subplot(2, 1, i, projection='3d') for i in [1,2] ]

    axarr[0].set_title( 'Source and Detector: Combined Dead Layer Estimate' )

    x = ret[ 'det_sectheta' ]
    y = ret[ 'source_sectheta' ]
    energies = ret[ 'energies' ]
    
    xflat = x.flatten()
    yflat = y.flatten()
    
    colors = [ 'r', 'b' ]

    for j in range(2):

        z = meas.meas.from_list( energies[j] )
        zflat = z.flatten()
        
        # plot it 
        for l in range(2):        
            plot_3d_errorbar( axarr[j],
                              [ x.x[l], y.x[l], z.x[l] ],
                              [ x.dx[l], y.dx[l], z.dx[l] ] )

            
        clf = linear_model.LinearRegression()

        mask = ~ ( meas.isnan( xflat ) | meas.isnan( yflat ) | meas.isnan( zflat ) )

        clf.fit( xflat.x[mask].reshape(-1,1), zflat.x[mask].reshape(-1,1) )          
        print( clf.coef_ )
        print( 'dl estimate: ' + str( np.abs( clf.coef_[0] ) / stopping_power_energies[j] )  )
        print( 'intercept: ' + str( clf.intercept_ ) )

        train2 = np.array( [ xflat.x[mask], yflat.x[mask] ] ).T
        
        clf2 = linear_model.LinearRegression()
        clf2.fit( train2, zflat.x[mask].reshape(-1,1) )          
        print( clf2.coef_ )
        print( 'dl estimate: ' + str( np.abs( clf2.coef_[0] ) /
                                      stopping_power_energies[j] )  )
        print( 'intercept: ' + str( clf.intercept_ ) )

        
    plt.show() 








def estimate_source_deadlayer_using_mu_calibration_lmfit() :
    
    ret = get_e_and_angles( dbmgr.flat, dbmgr.angled ) 

    f, axarr = plt.subplots( 2 )

    axarr[0].set_title( 'Source and Detector: Combined Dead Layer Estimate' )

    x = ret[ 'det_sectheta' ]
    y = ret[ 'source_sectheta' ]
    energies = ret[ 'energies' ]
    
    xflat = x.flatten()
    yflat = y.flatten()
    
    colors = [ 'r', 'b' ]
    fitcolors = [ 'g', 'y' ]

    for j in range(2):

        z = meas.meas.from_list( energies[j] )
        zflat = z.flatten()
                 
        jplt.plot( axarr[j], x.x[0], z.x[0],
                   xerr = x.dx[0], yerr = z.dx[0],
                   ylabel = r'$E$',
                   xlabel = r'$\sec \theta$',
                   color = colors[0],
                   leglabel = 'Flat' )
        
        jplt.plot( axarr[j], x.x[1], z.x[1],
                   xerr = x.dx[1], yerr = z.dx[1],
                   color = colors[1], leglabel = 'Angled' )
                   
        
        axarr[j].text( 0.1, 0.9, 'Peak %d' % (j,),
                       transform = axarr[j].transAxes,
                       fontsize=12,
                       verticalalignment='top' )


        labels = [ 'Flat Fit', 'Angled Fit' ]

        for k in range(2):
        
            cal = jstats.linear_calibration( x[k].x, z[k].x,
                                             dy = z[k].dx,
                                             ax = axarr[j],
                                             color = fitcolors[k],
                                             leglabel = labels[k],
            linestyle = '--' )
            
            if cal is not None:
                
                deadlayer = cal[0] / stopping_power_energies[j]
                                
                print( cal[2].fit_report() )  
                
                print( 'deadlayer: ' + str( deadlayer ) )
                
                axarr[j].text( 0.1, 0.2 + k * 0.1,
                               labels[k] + ' ' + r'DL = $ %d \pm %d $ nm ' % ( abs( deadlayer.x ),
                                                                       deadlayer.dx ),
                               transform = axarr[j].transAxes,
                               fontsize=12,
                               verticalalignment='top' )
                print( '\n\n' )
                
                
        jplt.add_legend( axarr[j], 1 ) 

        
    plt.show() 









    

def estimate_deadlayers_using_all_4_positions_mu_calibration_sklearn() :


    fig = plt.figure()
    axarr = [ fig.add_subplot(2, 1, i, projection='3d') for i in [1,2] ]

    axarr[0].set_title( 'Source and Detector: Combined Dead Layer Estimate' )

    
    sectheta_matrices = get_flattened_secant_matrices()

    keys = [ db.name for db in dbmgr.all_dbs ]

    energies = [ get_flattened_e_from_mu_calibration( db )
                 for db in dbmgr.all_dbs ]

    x = meas.meas.from_array( [ sectheta_matrices[ 'pu_238_' + key ][0] for key in keys ] ) 
    y = meas.meas.from_array( [ sectheta_matrices[ 'pu_238_' + key ][1] for key in keys ] )
    
    xflat = x.flatten()
    yflat = y.flatten()
    
    colors = [ 'r', 'g', 'b', 'y' ]


    for j in range( 2 ):

        z = meas.meas.from_array( energies[j] )
        zflat = z.flatten()
        
        # plot it 
        for l in range( len( x ) ):        
            plot_3d_errorbar( axarr[j],
                              [ x.x[l], y.x[l], z.x[l] ],
                              [ x.dx[l], y.dx[l], z.dx[l] ] )

            
        clf = linear_model.LinearRegression()

        mask = ~ ( meas.isnan( xflat ) | meas.isnan( yflat ) | meas.isnan( zflat ) )

        clf.fit( xflat.x[mask].reshape(-1,1), zflat.x[mask].reshape(-1,1) )          
        print( clf.coef_ )
        print( 'dl estimate: ' + str( np.abs( clf.coef_[0] ) / stopping_power_energies[j] )  )
        print( 'intercept: ' + str( clf.intercept_ ) )

        train2 = np.array( [ xflat.x[mask], yflat.x[mask] ] ).T
        
        clf2 = linear_model.LinearRegression()
        clf2.fit( train2, zflat.x[mask].reshape(-1,1) )          
        print( clf2.coef_ )
        print( 'dl estimate: ' + str( np.abs( clf2.coef_[0] ) /
                                      stopping_power_energies[j] )  )
        print( 'intercept: ' + str( clf.intercept_ ) )

        
    plt.show() 






    



def energy_from_mu_lmfit( params, mu, det_sectheta, source_sectheta,
                          db_name, x, i, j ) :

    a = params[ 'a_' + db_name + '_%d' % ( x, ) ].value
    b = params[ 'b_' + db_name + '_%d' % ( x, ) ].value

    det_deadlayer = params[ 'det_deadlayer' ].value

    source_constant = params[ 'source_constant_%d_%d' % (i,j) ].value

    return energy_from_mu( i, j, mu, det_sectheta, source_sectheta,
                           a, b, det_deadlayer, source_constant )
    





def energy_from_mu( i, j, mu, det_sectheta, source_sectheta,
                    a, b, det_deadlayer, source_constant ) :

    return ( a * mu + b
             + det_deadlayer * si_stopping_powers[ i, j ] * det_sectheta
             + source_constant * source_sectheta )

    




def objective( params, mu_matrices, secant_matrices, actual_energies ):

    # start_time = time.time() 
    
    resid = np.empty( ( 4, 3, 2, 32, 32 ) )

    db_names = [ db.name for db in dbmgr.all_dbs ]

    for db_num in range( len( db_names ) ) :

        db_name = db_names[ db_num ]

        source_names = [ 'pu_240',
                         'pu_238_' + db_name,
                         'cf_249' ]

        det_sectheta, source_sectheta = [ np.asarray( [ secant_matrices[ source ][k] 
                                                        for source in source_names ] )
                                          for k in range(2) ]
                
        for i in range(3):
            for j in range(2):
                for x in range(32):

                    computed_energy = energy_from_mu_lmfit(
                        params,
                        mu_matrices[ db_name ][ i ][ j ][ x ],
                        det_sectheta[i,x],
                        source_sectheta[i,x],
                        db_name, x, i, j )
                    
                
                    resid[ db_num, i, j, x ] = ( actual_energies[i,j] - computed_energy )

                    
    ret = resid.flatten()

    # print( 'objective: %f' % ( time.time() - start_time, ) )
           
    return ret 




                        



# do a fit of the form E = A * mu + b + s * sec(phi) + deadlayer_distance * si_stopping_power * sec(theta)  

def linear_calibration_on_each_x_strip():

    # prepare a giant Parameters() for the fit
    
    fit_params = lmfit.Parameters()


    fit_params.add( 'det_deadlayer', value = 100.0 ) #, min = 0.0 )

    
    for i in range( 3 ) :
        for j in range( 2 ) :
            fit_params.add( 'source_constant_%d_%d' % (i,j), value = 2.0 ) #, min = 0.0 )


    for db in dbmgr.all_dbs :
        for x in range( 32 ) :
            fit_params.add( 'a_' + db.name + '_%d' % ( x, ), value = 1.71 )
            fit_params.add( 'b_' + db.name + '_%d' % ( x, ), value = 700.0 )    

            
    # at this point we have many shared parameters. now read the data
    # and perform the minimization.
    
    mu_matrices = { db.name : [ [ db.get_all_mu_grids( 1 )[i][j].x
                                  for j in range(2) ]
                                for i in range(3) ]
                    for db in dbmgr.all_dbs }
    
    secant_matrices = geom.get_secant_matrices( 1 )

    secant_matrices = { k : v.x for k, v in secant_matrices.items() }
    
    result = lmfit.minimize( objective,
                             fit_params,
                             args = ( mu_matrices, secant_matrices, peak_energies ),
                             nan_policy = 'omit' )

    lmfit.report_fit( result.params, show_correl = 0 ) 
        
                    
            
            




    
  

    


    
    
# example_estimate_for_one_source()



# estimate_source_deadlayer_using_mu_calibration_lmfit()

# estimate_deadlayers_using_all_4_positions_mu_calibration_sklearn()

# preview_secant_differences()


linear_calibration_on_each_x_strip()
