



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

