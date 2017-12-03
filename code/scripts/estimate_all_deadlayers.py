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


import heapq

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


# https://physics.nist.gov/PhysRefData/Star/Text/ASTAR-t.html
# 5.12368
# 5.16817
# 5.4563
# 5.49903
# 5.7595
# 5.81310

si_stopping_powers = np.array( [ [ 6.080E+02, 6.046E+02 ],
                                 [ 5.832E+02, 5.802E+02 ],
                                 [ 5.626E+02, 5.591E+02 ] ] )

si_dioxide_stopping_powers = np.array( [ [ 6.489E+02, 6.452E+02 ],
                                         [ 6.222E+02, 6.190E+02 ],
                                         [ 6.001E+02, 5.963E+02   ] ] )

density_si = 2.328 # g / cm^2
density_si_dioxide = 2.65

si_stopping_powers *= density_si * 1000 * 100 / 1e9   # convert to keV / nm from mev / cm
si_dioxide_stopping_powers *= density_si_dioxide * 1000 * 100 / 1e9


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



def nth_largest(n, iter):
    return heapq.nlargest(n, iter)[-1]






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
                          db_name, x, i, j, vary_det_deadlayer = 0,
                          const_source_deadlayer = 0,
                          quadratic_source = 0,
                          calibrate_each_pixel = 0,
                          y = 0,
                          compute_weight = 0 ) :

    if not calibrate_each_pixel : 
        a = params[ 'a_' + db_name + '_%d' % ( x, ) ].value.item()
        b = params[ 'b_' + db_name + '_%d' % ( x, ) ].value.item()

    else: 
        a = params[ 'a_' + db_name + '_%d_%d' % ( x,y ) ].value.item()
        b = params[ 'b_' + db_name + '_%d_%d' % ( x,y ) ].value.item()

        
    if not vary_det_deadlayer:
        det_constant = params[ 'det_deadlayer' ].value.item() * si_stopping_powers[ i, j ]

    else:
        det_constant = params[ 'det_constant_%d_%d' % (i,j) ].value.item()

        
    source_constant = params[ 'source_constant_%d_%d' % (i,j) ].value

    # shift by a constant term if this option is enabled.
    if const_source_deadlayer :
        b += source_constant
        source_constant = 0

    if quadratic_source :
        source_constant2 = params[ 'source_constant2_%d_%d' % (i,j) ].value.item()

    else:
        source_constant2 = 0.0

        
    energy = energy_from_mu( mu.x, det_sectheta.x, source_sectheta.x,
                             a, b, det_constant, source_constant, source_constant2 )

   
    if compute_weight : 
        weight = 1 / ( a * mu.dx )
        return ( energy, weight ) 

    return energy 

                         





def energy_from_mu( mu, det_sectheta, source_sectheta,
                    a, b, det_constant, source_constant, source_constant2 ) :
    
    return ( a * mu + b
             + det_constant * det_sectheta
             + source_constant * source_sectheta 
             + source_constant2 * (source_sectheta ** 2) )
    




def objective( params, mu_matrices, secant_matrices, actual_energies,
               dbs, source_indices, vary_det_deadlayer = 0,
               const_source_deadlayer = 0,
               quadratic_source = 0,
               calibrate_each_pixel = 0 ):

    # start_time = time.time() 
    
    resid = np.zeros( ( 4, 3, 2, 32, 32 ) )

    db_names = [ db.name for db in dbs ]

    for db_num in range( len( db_names ) ) :

        db_name = db_names[ db_num ]

        source_names = [ 'pu_240',
                         'pu_238_' + db_name,
                         'cf_249' ]

        det_sectheta, source_sectheta = [ [ secant_matrices[ source ][k] 
                                            for source in source_names ]
                                          for k in range(2) ]

        
        for i in range( len( source_indices ) ) :
            for j in source_indices[i] :
                for x in range(32) :

                    # if not calibrate_each_pixel :

                    computed_energy, weight = energy_from_mu_lmfit(
                        params,
                        mu_matrices[ db_name ][ i ][ j ][ x ],
                        det_sectheta[i][x],
                        source_sectheta[i][x],
                        db_name, x, i, j, vary_det_deadlayer,
                        const_source_deadlayer,
                        quadratic_source,
                        compute_weight = 1 )
                    
                    residual = actual_energies[i,j] - computed_energy
                                        
                    resid[ db_num, i, j, x ] = residual * weight
                    
                    # else:

                    #     for y in range( 32 ) :
                        
                    #         computed_energy = energy_from_mu_lmfit(
                    #             params,
                    #             mu_matrices[ db_name ][ i ][ j ][ x,y ],
                    #             det_sectheta[i,x,y],
                    #             source_sectheta[i,x,y],
                    #             db_name, x, i, j, vary_det_deadlayer,
                    #             const_source_deadlayer,
                    #             quadratic_source,
                    #             calibrate_each_pixel,
                    #             y = y )
                            
                    #         resid[ db_num, i, j, x, y ] = ( actual_energies[i,j]
                    #                                         - computed_energy )

                                        
    ret = resid.flatten()

    # print( 'objective: %f' % ( time.time() - start_time, ) )
           
    return ret 




                        



# do a fit of the form E = A * mu + b + s * sec(phi) + deadlayer_distance * si_stopping_power * sec(theta)  

def linear_calibration_on_each_x_strip( dbs,
                                        source_indices,
                                        vary_det_deadlayer = 0,
                                        const_source_deadlayer = 0,
                                        quadratic_source = 0,
                                        calibrate_each_pixel = 0 ):

    
    # these are all the dbs we will look at, for quick scripting to look
    # at subsets of the data. these get passed to objective.

    # dbs = dbmgr.all_dbs
    # dbs  = dbmgr.normal_dbs
    # dbs = [ dbmgr.centered, dbmgr.angled ]
    # dbs = [ dbmgr.flat ] 
    # dbs = [ dbmgr.centered ] 
    
    # these are the sources we will look at. neglect the others.
    # normal operation is [ [0,1], [0,1], [0,1] ]

    # source_indices = [ [0,1], [0,1], [1] ]
    # source_indices = [ [], [0,1], [] ]

    
    # prepare a giant Parameters() for the fit
    
    fit_params = lmfit.Parameters()

    if not vary_det_deadlayer : 
        fit_params.add( 'det_deadlayer', value = 0, vary = 1 ) #, min=0 ) #, min = 0.0 )

    for i in range( len( source_indices ) ) :
        for j in source_indices[i] :

            fit_params.add( 'source_constant_%d_%d' % (i,j), value = 0.0,
                            vary = ( not vary_det_deadlayer
                                     and dbmgr.angled in dbs ) )
            if quadratic_source :
                fit_params.add( 'source_constant2_%d_%d' % (i,j), value = 0.0, vary = 1  )
                    
            if vary_det_deadlayer :
                fit_params.add( 'det_constant_%d_%d' % (i,j), value = 0.0, vary = 1  )


    
    for db in dbs :

        if not calibrate_each_pixel : 
        
            for x in range( 32 ) :
                fit_params.add( 'a_' + db.name + '_%d' % ( x, ), value = 1.99 )
                fit_params.add( 'b_' + db.name + '_%d' % ( x, ), value = -200 )    

        else :
            for x in range( 32 ) :
                for y in range( 32 ) :
                    fit_params.add( 'a_' + db.name + '_%d_%d' % ( x,y ), value = 1.99 )
                    fit_params.add( 'b_' + db.name + '_%d_%d' % ( x,y ), value = -200 )    

                
    # at this point we have many shared parameters. now read the data
    # and perform the minimization.

    peak_matrices = { db.name : db.get_all_peak_grids( 1 )
                      for db in dbs }
    
    # mu_matrices = { db.name : db.get_all_mu_grids( 1 )
    #                 for db in dbs }
    
    secant_matrices = geom.get_secant_matrices( 1 )

    result = lmfit.minimize( objective,
                             fit_params,
                             args = ( mu_matrices, secant_matrices, peak_energies,
                                      dbs, source_indices,
                                      vary_det_deadlayer,
                                      const_source_deadlayer,
                                      quadratic_source,
                                      calibrate_each_pixel ),
                             nan_policy = 'omit' )

    lmfit.report_fit( result.params, show_correl = 0 )

    print( 'reduced chisq: ' + str( result.redchi ) )

    plot_energy_vs_sectheta( result, secant_matrices, mu_matrices,
                             dbs, source_indices,
                             vary_det_deadlayer,
                             quadratic_source ) 
    
    # plot_results( result,
    #               secant_matrices, mu_matrices,
    #               dbs[0], 17, source_indices,
    #               vary_det_deadlayer,
    #               const_source_deadlayer,
    #               quadratic_source,
    #               calibrate_each_pixel )


    




    


def plot_results( lmfit_result,
                  secant_matrices, mu_matrices,
                  test_db, row,
                  source_indices,
                  vary_det_deadlayer = 0,
                  const_source_deadlayer = 0,
                  quadratic_source = 0,
                  calibrate_each_pixel = 0 ) :

    if calibrate_each_pixel :
        return 
    
    f, axarr = plt.subplots( 3, 2 )
    
    axarr[0][0].set_title( r'Absolute calibrated $\mu$ vs. $\sec \theta $ For Each Peak' ) 

    for i in range(len( source_indices ) ) :
        for j in source_indices[i] :
                
            if i == 0 :
                source = 'pu_240'
                
            elif i == 1 :
                source = 'pu_238_' + test_db.name
            
            elif i == 2 :
                source = 'cf_249'
                            
            x = secant_matrices[source][0][row]
            y = secant_matrices[source][1][row]
            z = mu_matrices[ test_db.name ][i][j][row]

            # print( source ) 
            # print( 'x: ' + str( x ) )
            # print( 'y: ' + str( y ) )
            # print( 'z: ' + str( z ) )
            # print( '\n\n' )
            
            
            axarr[i,j].errorbar( x.x, z.x, xerr = x.dx, yerr = z.dx,
                                 fmt = 'o', color='b'  )
            
            test_id = '_' + test_db.name + '_%d' % ( row, )
            
            a = lmfit_result.params[ 'a' + test_id ].value
            b = lmfit_result.params[ 'b' + test_id ].value
            
            if vary_det_deadlayer:
                det_constant = lmfit_result.params[ 'det_constant_%d_%d' % (i,j) ].value
                print( 'effective dl: ' + str ( det_constant / si_stopping_powers[i][j] ) )
                
            else:
                dl = lmfit_result.params[ 'det_deadlayer' ].value
                det_constant = dl * si_stopping_powers[i][j]

            source_constant = np.asscalar( lmfit_result.params[ 'source_constant_%d_%d' % (i,j) ] )
                
            yfit = ( peak_energies[i][j] - b 
                     - det_constant * x.x ) / a 
            
            if const_source_deadlayer :
                yfit -= source_constant / a 
            else:
                yfit -= source_constant * y.x / a

            if quadratic_source :
                yfit -= ( lmfit_result.params[ 'source_constant2_%d_%d' % (i,j) ]
                          * (y ** 2 ) / a ) 

            mask = ~ meas.isnan( z )

            # print(x)
            # print( x[mask] )

            # print(yfit)
            # print(yfit[mask])
            
            axarr[i,j].plot( x[ mask ].x, yfit[ mask ], c = 'r' ) 
                    

    plt.show()

    # fig = plt.figure()
    # axarr = [ fig.add_subplot(2, 1, i, projection='3d') for i in [1,2] ]

    


    


    

    

def plot_energy_vs_sectheta( lmfit_result, secant_matrices, mu_matrices,
                             dbs, source_indices,
                             vary_det_deadlayer = 0,
                             quadratic_source = 0 ) :

    
    f, axarr = plt.subplots( 3, 2 )
    
    axarr[0][0].set_title( r'Absolute $E$ vs. $\sec \theta $ For Each Peak' )


    for i in range(len( source_indices ) ) :
        for j in source_indices[i] :
            
                           
            energies = np.empty( ( len(dbs), 32, 32 ) )
            
            for d in range( len( dbs ) ) :

                db = dbs[ d ] 

                if i == 0 :
                    source = 'pu_240'
                
                elif i == 1 :
                    source = 'pu_238_' + db.name
            
                elif i == 2 :
                    source = 'cf_249'

                    
                if vary_det_deadlayer:
                    det_constant = lmfit_result.params[ 'det_constant_%d_%d' % (i,j) ].value.item()
                    print( 'effective dl: ' + str ( det_constant / si_stopping_powers[i][j] ) )
                    
                else:
                    dl = lmfit_result.params[ 'det_deadlayer' ].value.item()
                    det_constant = dl * si_stopping_powers[i][j]
                    
                source_constant = lmfit_result.params[ 'source_constant_%d_%d' % (i,j) ].value
                    
                for row in range(32):
                    
                    x = secant_matrices[source][0][row,:]
                    y = secant_matrices[source][1][row,:]
                    z = mu_matrices[ db.name ][i][j][row,:]                   
                                        
                    test_id = '_' + db.name + '_%d' % ( row, )

                    a = lmfit_result.params[ 'a' + test_id ].value.item()
                    b = lmfit_result.params[ 'b' + test_id ].value.item()

                    E = a * z + b
                                       
                    energies[ d, row, : ] = E.x

                    calibrated_E = energy_from_mu_lmfit( lmfit_result.params, z,
                                                         x,
                                                         y,
                                                         db.name, row, i, j,
                                                         vary_det_deadlayer = vary_det_deadlayer,
                                                         quadratic_source = quadratic_source )


                    Efit = ( peak_energies[i][j]
                             - ( calibrated_E - E.x ) )

                    axarr[i,j].errorbar( x.x, E.x, xerr = x.dx, yerr = E.dx,
                                         fmt = 'o', color='b', zorder = 1  )


                    mask = ~ meas.isnan( z )
                                
                    
                    axarr[i,j].plot( x[ mask ].x, Efit[ mask ], c = 'r', zorder = 2 ) 


            # remove outliers (not from fit, just the plot )
            flattened_energies = energies.flatten()
            mask = ~ np.isnan( flattened_energies )            
            sorted_energies = sorted( flattened_energies[ mask ] )
            
            newmin = sorted_energies[10]
            newmax = sorted_energies[ len(sorted_energies) - 5 ]
            
            axarr[i][j].set_ylim( newmin - 10, newmax + 10 )

                    
    plt.show()

                             




    
  

    


    
    
# example_estimate_for_one_source()



# estimate_source_deadlayer_using_mu_calibration_lmfit()

# estimate_deadlayers_using_all_4_positions_mu_calibration_sklearn()

# preview_secant_differences()


linear_calibration_on_each_x_strip( [ dbmgr.angled ],
                                    [ [0,1], [0,1], [1] ],
                                    vary_det_deadlayer = 1,
                                    quadratic_source = 0,
                                    calibrate_each_pixel = 0 )
