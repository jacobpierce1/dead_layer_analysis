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


from scipy.interpolate import interp1d

import heapq

from mpldatacursor import datacursor 

import matplotlib.pyplot as plt
import pandas as pd
import numpy as np
import scipy.optimize

from sklearn import linear_model
import lmfit







epsilon_0 = 3.67 / 1000  # keV / electron-hole pair 



peak_energies = np.array( [ [ 5123.68, 5168.17 ],
                            [ 5456.3, 5499.03 ],
                            [ 5759.5, 5813.10,  ] ] )   # all in keV

flattened_peak_calibration_energies = peak_energies.flatten()



# these are the stopping powers at the above peak energies.

si_stopping_powers = np.array( [ [ 6.080E+02, 6.046E+02 ],
                                 [ 5.832E+02, 5.802E+02 ],
                                 [ 5.626E+02, 5.591E+02 ] ] )

si_dioxide_stopping_powers = np.array( [ [ 6.489E+02, 6.452E+02 ],
                                         [ 6.222E+02, 6.190E+02 ],
                                         [ 6.001E+02, 5.963E+02   ] ] )



# more data for an interpolation.

si_interp_energies = np.array( [ 5.050E+00, 5.090E+00, 5.130E+00, 5.170E+00,
                                 5.210E+00, 5.250E+00, 5.290E+00, 5.330E+00,
                                 5.370E+00, 5.410E+00, 5.450E+00, 5.490E+00,
                                 5.530E+00, 5.570E+00, 5.610E+00, 5.650E+00,
                                 5.690E+00, 5.730E+00, 5.770E+00, 5.810E+00,
                                 5.81310,] )

si_interp_stopping_powers = np.array( [ 6.138E+02, 6.107E+02, 6.075E+02, 6.044E+02,
                                        6.014E+02, 5.983E+02, 5.953E+02, 5.924E+02,
                                        5.894E+02, 5.866E+02, 5.837E+02, 5.809E+02,
                                        5.781E+02, 5.753E+02, 5.726E+02, 5.699E+02,
                                        5.672E+02, 5.645E+02, 5.619E+02, 5.593E+02,
                                        5.591E+02 ] ) 


# big_si_interp_energies = np.linspace( 0.1, 6.0, 60 )

# big_si_interp_stopping_powers = np.array( [ 9.107E+02, 1.232E+03, 1.364E+03, 1.412E+03,
#                                             1.420E+03, 1.408E+03, 1.386E+03, 1.358E+03,
#                                             1.328E+03, 1.296E+03, 1.265E+03, 1.234E+03,
#                                             1.204E+03, 1.175E+03, 1.148E+03, 1.121E+03,
#                                             1.095E+03, 1.070E+03, 1.046E+03, 1.024E+03,
#                                             1.002E+03, 9.815E+02, 9.618E+02, 9.429E+02,
#                                             9.248E+02, 9.073E+02, 8.904E+02, 8.742E+02,
#                                             8.584E+02, 8.432E+02, 8.284E+02, 8.141E+02,
#                                             8.003E+02, 7.868E+02, 7.738E+02, 7.611E+02,
#                                             7.487E+02, 7.367E+02, 7.251E+02, 7.137E+02,
#                                             7.027E+02, 6.919E+02, 6.815E+02, 6.715E+02,
#                                             6.618E+02, 6.524E+02, 6.433E+02, 6.346E+02,
#                                             6.261E+02, 6.179E+02, 6.099E+02, 6.021E+02,
#                                             5.946E+02, 5.873E+02, 5.802E+02, 5.732E+02,
#                                             5.665E+02, 5.600E+02, 5.536E+02, 5.474E+02 ] )




density_si = 2.328 # g / cm^2
density_si_dioxide = 2.65

si_stopping_powers *= density_si * 1000 * 100 / 1e9   # convert to keV / nm from mev / cm
si_dioxide_stopping_powers *= density_si_dioxide * 1000 * 100 / 1e9
si_interp_stopping_powers *= density_si * 1000 * 100 / 1e9 


si_interp_energies *= 1000   # keV / MeV




# stopping power for the 5456.3 and 5499.03 keV peaks alphas
# https://physics.nist.gov/cgi-bin/Star/ap_table-t.pl
# stopping_power_energies = np.array( [ 5.802E+02 , 5.802E+02 ] )  # MeV cm^2 / g 




# constants, need density of each deadlayer id




def nth_largest(n, iter):
    return heapq.nlargest(n, iter)[-1]






    


# this object stores all the information required to state the model
# being used to describe the relationship between energy entering
# detector and the secant of penetration angle

class dead_layer_model_params( object ) :

    def __init__( self,
                  vary_det_deadlayer = 0,
                  quadratic_source = 0,
                  quadratic_det = 0,
                  calibrate_each_pixel = 0,
                  interp_stopping_power = 1,
                  mu = 1,
                  average_over_source = 1,
                  pulse_height_defect = 0,
                  # ignore_outer_pixels = 0,
                  fstrips = None,
                  bstrips = None ) :

        # various options that change the model the data is fit
        # to. see energy_from_mu_lmfit() for what they do.
        self.vary_det_deadlayer = vary_det_deadlayer
        self.quadratic_source = quadratic_source
        self.quadratic_det = quadratic_det
        self.calibrate_each_pixel = calibrate_each_pixel
        self.pulse_height_defect = pulse_height_defect
        
        # construct a stopping power interpolation.
        # it is not that much extra computation to use this
        # instead of the less complex approximation, so
        # i recommend using it. 
        self.interp_stopping_power = interp_stopping_power

        # bool: 1 to use mu values of each peak, 0 to use
        # peak values for each peak. i recommend using mu values.
        self.mu = mu

        # use a modified version of secant of penetration angle
        # for the sources
        self.average_over_source = average_over_source
        
        self.si_stopping_power_interpolation = None

        # self.ignore_outer_pixels = ignore_outer_pixels

        # these are the rows that will be considered when
        # performing the optimization. others will be ignored.
        
        if fstrips is None:
            self.fstrips = np.arange(32)
        else:
            self.fstrips = fstrips

        if bstrips is None:
            self.bstrips = np.arange(32)
        else:
            self.bstrips = bstrips
        
        return None

        

    

    




    

# interpolate stopping power of alphas in Si using data
# in the range emin, emax. data is from astar program of NIST.
# uses scipy.interpolate.interp1d. if from_0_keV is 1, then
# interpolate over a lot of data to obtain large interpolation
# and set the 

def construct_si_stopping_power_interpolation( from_0_keV = 0, plot = 0 ) :

    data = np.loadtxt( '../../data/stopping_power_data/alpha_stopping_power_si.txt',
                       skiprows = 10, unpack = 1 )

    emax = peak_energies.max()

    energy = data[0] * 1000
    
    energy = energy[ energy <= emax ] 
    
    stopping_power = data[3][ 0 : len( energy ) ]
    stopping_power *= density_si * 1000 * 100 / 1e9


    # add particular data points of interest to the interpolation

    energy = np.append( energy, si_interp_energies )
    stopping_power = np.append( stopping_power, si_interp_stopping_powers )
    
    interp = scipy.interpolate.interp1d( energy, stopping_power, kind = 'cubic' )

    
    if plot :

        ax = plt.axes()

        interp_axis = np.linspace( min( energy ),
                                   max( energy ),
                                   100 )
        
        ax.scatter( energy, stopping_power, color='r' )

        ax.plot( interp_axis,
                 interp( interp_axis ),
                 color = 'b' )
        
        plt.show()

        return 1
        

    return interp
    






    


# the energy as measured by the detector, energy_det = m * mu + b,
# does not equal the initial energy of the alpha because of several
# possible losses. the rest of the function adds possible losses and
# returns a prediction for the actual initial energy.  params are the
# parameters of the fit which change as the regression is run,
# model_params are constant parameters saying what type of fit we are
# using (e.g. with a quadratic term in sec theta )

def energy_from_mu_lmfit( params,
                          mu, det_sectheta, source_sectheta,
                          db_name, x, i, j,
                          model_params,
                          compute_weight = 0 ) :

    # print( 'test' ) 

    a = params[ 'a_' + db_name + '_%d' % ( x, ) ].value.item()
    b = params[ 'b_' + db_name + '_%d' % ( x, ) ].value.item()
    source_constant = params[ 'source_constant_%d_%d' % (i,j) ].value
    
    energy_det = a * mu + b 

    # in this case, we are using the actual depth of the
    # detector dead layer as a parameter.
    
    if not model_params.vary_det_deadlayer:

        det_deadlayer = params[ 'det_deadlayer' ].value.item()
        
        if model_params.si_stopping_power_interpolation is not None :

            S = model_params.si_stopping_power_interpolation( peak_energies[i][j]
                                                 - source_constant * source_sectheta.x )
            
            det_constant = det_deadlayer * S 
            
        else:
            det_constant = det_deadlayer * si_stopping_powers[ i, j ]

    else:
        det_constant = 0


        
    # compute how much energy was lost in the source and det deadlayers.
       
    combined_deadlayer_losses =  ( det_constant * det_sectheta.x
                                   + source_constant * source_sectheta.x )

    if model_params.quadratic_source :
        source_constant2 = params[ 'source_constant2_%d_%d' % (i,j) ].value.item()

        combined_deadlayer_losses += source_constant2 * (source_sectheta.x ** 2)

    if model_params.quadratic_det :
        det_constant2 = params[ 'det_constant2_%d_%d' % (i,j) ].value.item()

        combined_deadlayer_losses += det_constant2 * ( det_sectheta.x - 1 ) ** 2
        

    # ignoring the pulse height defect, just add back the energy
    # that was lost in the source and detector dead layers.
    
    if not model_params.pulse_height_defect :
        
        energy = energy_det + combined_deadlayer_losses

        
    # otherwise we compute it using the model of lennard et al (1986) 
    
    else:
        k = params['k'].value.item()

        energy = np.empty( len( mu ) )

        integrand = lambda E : 1 / ( epsilon_0 - k * si_stopping_power_interpolation( E ) )
                
        for k in range( len( mu ) ) :

            # print( combined_deadlayer_losses ) 
            
            energy[k] = scipy.integrate.quad( integrand,
                                              model_params.si_stopping_power_interpolation_emin,
                                              peak_energies[i,j] - combined_deadlayer_loss )[0]

        energy *= epsilon_0 


    # the weight is the uncertainty of the total computed energy
    # under the model. under all circumstances, this is dominated
    # by the uncertainty in mu.
        
    if compute_weight : 
        weight = 1 / ( a * mu.dx )
        return ( energy, weight ) 

    return energy 










def objective( params, mu_matrices, secant_matrices, actual_energies,
               dbs, source_indices, model_params ):

    # start_time = time.time()
    
    resid = np.empty( (len(dbs)
                       * len( jutils.flatten_list( source_indices) )
                       * len( model_params.fstrips ) 
                       * len( model_params.bstrips ) , ) ) # np.zeros( ( len(dbs), 3, 2, 32, 32 ) )

    db_names = [ db.name for db in dbs ]

    # keep track of where the 1D residual array has been
    # filled up to.
    resid_idx = 0
    num_bstrips = len( model_params.bstrips )

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
                for x in model_params.fstrips :

                    computed_energy, weight = energy_from_mu_lmfit(
                        params,
                        mu_matrices[ db_name ][i][j][x][ model_params.bstrips ],
                        det_sectheta[i][x][ model_params.bstrips ],
                        source_sectheta[i][x][ model_params.bstrips ],
                        db_name, x, i, j,
                        model_params,
                        compute_weight = 1 ) 

                    # strange bug: when adding meas and scalar, scalar
                    # must be on the right. must fix asap.

                    residual = - computed_energy + actual_energies[i,j]

                    # resid[ db_num, i, j, x ] = residual.x * weight
                    resid[ resid_idx : resid_idx + num_bstrips ] = residual.x * weight
                    resid_idx += num_bstrips 
                    
    # ret = resid.flatten()
    
    # print( 'objective: %f' % ( time.time() - start_time, ) )
           
    # return ret 
    return resid









# do a fit of the form E = A * mu + b + s * sec(phi) +
# deadlayer_distance * si_stopping_power * sec(theta)

def linear_calibration_on_each_x_strip( dbs,
                                        source_indices,
                                        model_params,
                                        annotate = 0,
                                        subtitle = '',
                                        view_pixel = None ) : 


    
    # next, read the data that will go into the regression. 
    if model_params.mu : 
        mu_matrices = { db.name : db.get_all_mu_grids( 1 )
                        for db in dbs }

        peak_matrices = mu_matrices

    else:
        peak_matrices = { db.name : db.get_all_peak_grids( 1 )
                          for db in dbs }

    # filter out bad data 
    for db_name in [ db.name for db in dbs ] :
        for i in range(3) :
            for j in range( 2 ) :
                mask = peak_matrices[db_name][i][j].dx > 2
                peak_matrices[db_name][i][j][ mask ] = meas.nan
                
                
    secant_matrices = geom.get_secant_matrices( compute_source_sectheta = 1,
                                                average_over_source = model_params.average_over_source,
                                                reset = not model_params.average_over_source )

    # cut out large sec theta 
    for key, val in secant_matrices.items() :
        mask = ( val[0].x >= 1.45 )
        secant_matrices[key][0][mask] = meas.nan

        
    if view_pixel is not None :
        plot_energy_vs_sectheta( None, secant_matrices, peak_matrices,
                                 dbs, source_indices, model_params,
                                 subtitle = subtitle,
                                 view_pixel = view_pixel )
        return 1

    
        
    if model_params.pulse_height_defect :
        interp_stopping_power = 1

    
    # prepare a giant lmfit.Parameters() for the fit
    
    fit_params = lmfit.Parameters()

    # this says whether we will put a lower bound on the following fit
    # parameters.

    # set_lower_bound = model_params.quadratic_source or model_params.quadratic_det
    set_lower_bound = 0 

    if not model_params.vary_det_deadlayer :
        fit_params.add( 'det_deadlayer', value = 100, vary = 1 ) #, min=0 ) #, min = 0.0 )
        if set_lower_bound :
            fit_params[ 'det_deadlayer' ].min = 0

        
    for i in range( len( source_indices ) ) :
        for j in source_indices[i] :

            fit_params.add( 'source_constant_%d_%d' % (i,j), value = 6.0 )
            
            if set_lower_bound :
                fit_params[ 'source_constant_%d_%d' % (i,j) ].min = 0
            
            if model_params.quadratic_source :
                fit_params.add( 'source_constant2_%d_%d' % (i,j), value = 0.0, vary = 1  )

            if model_params.quadratic_det :
                fit_params.add( 'det_constant2_%d_%d' % (i,j), value = 0.0, vary = 1  )
                
    
    for db in dbs :
        for x in model_params.fstrips :
            fit_params.add( 'a_' + db.name + '_%d' % ( x, ), value = 1.99 )
            fit_params.add( 'b_' + db.name + '_%d' % ( x, ), value = -50 )    

            
    # this is the k parameter of the lennard model. rather
    # than assuming their value for alphas, we let it be a free
    # parameter and will compare to theirs.
    
    if model_params.pulse_height_defect :
        fit_params.add( 'k', value = 0 )



    
    # construct an interpolation if required.
    
    if model_params.interp_stopping_power :
        interp = construct_si_stopping_power_interpolation( 1, 0 )
        emin = min( interp.x ) 

        model_params.si_stopping_power_interpolation = interp 
        model_params.si_stopping_power_interpolation_emin = emin
            
    else:
        model_params.si_stopping_power_interpolation = None

        
    result = lmfit.minimize( objective,
                             fit_params,
                             args = ( peak_matrices, secant_matrices, peak_energies,
                                      dbs, source_indices, model_params ),
                             nan_policy = 'omit' )

    
    lmfit.report_fit( result.params, show_correl = 0 )

    print( 'reduced chisq: ' + str( result.redchi ) )
    print( 'chisq: ' + str( result.chisqr ) )
    print( 'ndata: ' + str( result.ndata ) )
    print( 'nfree: ' + str( result.nfree ) ) 

    plot_energy_vs_sectheta( result, secant_matrices, peak_matrices,
                             dbs, source_indices,
                             model_params,
                             annotate = annotate,
                             subtitle = subtitle,
                             view_pixel = view_pixel ) 
    

    




    

    

    

def plot_energy_vs_sectheta( lmfit_result, secant_matrices, mu_matrices,
                             dbs, source_indices,
                             model_params,
                             annotate = 0,
                             subtitle = '',
                             view_pixel = None ) :
    
    f, axarr = plt.subplots( 3, 2 )

    # \mathrm instead of \text
    if view_pixel is None :
        
        f.suptitle( r'Absolute $ E_\mathrm{det}$ vs. $\sec ( \theta_\mathrm{det} ) $ For Each Peak'
                    + '\n' + subtitle + ', ' + r'$ \tilde{\chi}^2 = '
                    + '%.2f' % ( lmfit_result.redchi , ) + '$' )


    # f.text( 0.95, 0.95, r'$ \tilde{\chi}^2 = ' + str( lmfit_result.redchi ) + '$' )
    
    # create common axis labels
    
    f.add_subplot(111, frameon=False)

    # hide tick and tick label of the big axes

    plt.tick_params(labelcolor='none', top='off', bottom='off', left='off', right='off')
    plt.grid(False)
    plt.xlabel( r'$sec ( \theta_\mathrm{det} ) $' )

    if view_pixel is not None : 
        plt.ylabel( r'$E_\mathrm{det}$ (keV)' )

    else : 
        plt.ylabel( r'$\mu$ (channels)' )
        
    titles = [ 'Pu 240', 'Pu 238', 'Cf 249' ]
    
    for i in range(len( source_indices ) ) :

        for j in source_indices[i] :

            axarr[i,j].set_title( titles[i] )
            
            # energies = np.empty( ( len(dbs), 32, 32 ) )
            if view_pixel is None :
                
                if model_params.vary_det_deadlayer:
                    det_constant = lmfit_result.params[ 'source_constant_%d_%d' % (i,j) ].value.item()
                    print( 'effective dl: ' + str ( det_constant / si_stopping_powers[i][j] ) )

                for d in range( len( dbs ) ) :

                    db = dbs[ d ] 

                    if i == 0 :
                        source = 'pu_240'

                    elif i == 1 :
                        source = 'pu_238_' + db.name

                    elif i == 2 :
                        source = 'cf_249'


                    # else:
                    #     dl = lmfit_result.params[ 'det_deadlayer' ].value.item()
                    #     det_constant = dl * si_stopping_powers[i][j]

                    # source_constant = lmfit_result.params[ 'source_constant_%d_%d' % (i,j) ].value

                    for row in model_params.fstrips :

                        x = secant_matrices[source][0][row][ model_params.bstrips ]
                        y = secant_matrices[source][1][row][ model_params.bstrips ]
                        z = mu_matrices[ db.name ][i][j][row][ model_params.bstrips ]

                        test_id = '_' + db.name + '_%d' % ( row, )

                        a = lmfit_result.params[ 'a' + test_id ].value.item()
                        b = lmfit_result.params[ 'b' + test_id ].value.item()

                        E = a * z + b

                        # energies[ d, row, : ] = E.x

                        calibrated_E = energy_from_mu_lmfit( lmfit_result.params,
                                                             z, x, y,
                                                             db.name, row, i, j,
                                                             model_params )

                        Efit = peak_energies[i][j] - ( calibrated_E.x - E.x )

                        if not annotate: 
                            axarr[i,j].errorbar( x.x, E.x, xerr = x.dx, yerr = E.dx,
                                                 color='b', zorder = 1,
                                                 ls = 'none' )
                            
                        else:
                            for a in range( len( model_params.bstrips ) ) :

                                if not meas.isnan( z[a] ) :

                                    label = '(%s,%d,%d,%.1f)' % (db.name, row, model_params.bstrips[a],
                                                               z[a].x ) 

                                    axarr[i,j].scatter( [ x.x[a] ] , [ E.x[a] ], color='b',
                                                        zorder = 1, label = label )


                        mask = ~ meas.isnan( z )

                        axarr[i,j].plot( x[ mask ].x, Efit[ mask ], c = 'r', zorder = 2 ) 


                # # remove outliers (not from fit, just the plot )

                # flattened_energies = energies.flatten()
                # mask = ~ np.isnan( flattened_energies )            
                # sorted_energies = sorted( flattened_energies[ mask ] )

                # newmin = sorted_energies[10]
                # newmax = sorted_energies[ len(sorted_energies) - 5 ]

                # # axarr[i][j].set_ylim( newmin - 10, newmax + 10 )

            else :

                db = view_pixel[0]
                row = view_pixel[1]

                if i == 0 :
                    source = 'pu_240'

                elif i == 1 :
                    source = 'pu_238_' + db.name
                        
                elif i == 2 :
                    source = 'cf_249'
                
                x = secant_matrices[source][0][row][ model_params.bstrips ]
                y = secant_matrices[source][1][row][ model_params.bstrips ]
                z = mu_matrices[ db.name ][i][j][row][ model_params.bstrips ]

                print( (i,j) )
                print( 'mu / peak values: ' )
                print( z )
                print('\n' )
                
                axarr[i,j].errorbar( x.x, z.x, yerr = z.dx,
                                     color='b', zorder = 1, ls = 'none' ) #, label = 'test' )

                                

    if annotate :
        datacursor( formatter = '{label}'.format )
    
    plt.show()

    



# def ignore_coordinates_predicate( i, j ) :
#     return i>1 and i<30 and j>1 and j<30


    
  

# dbs = [ dbmgr.centered, dbmgr.angled ]
# dbs = [ dbmgr.angled, dbmgr.centered, dbmgr.flat ]
# dbs = dbmgr.all_dbs
# dbs = [ dbmgr.angled, dbmgr.centered ]
# dbs = [ dbmgr.angled ] 
dbs = [ dbmgr.centered, dbmgr.angled, dbmgr.moved ] # , dbmgr.flat ] 



# source_indices = [ [0,1], [0,1], [1] ]
source_indices = [ [0,1], [0,1], [1] ]



# subtitle = 'Angled Dataset'
subtitle = 'Centered, Angled, Moved'


fstrips = np.arange( 3, 30 )

# a = 2
# b = 31
# c = 27

# fstrips = np.zeros( b-a - 1, dtype='int' )

# fstrips[ 0 : c-a ] = np.arange( a,c)
# fstrips[ c-a : b-a-1 ] = np.arange( c+1, b )
# print( fstrips )


model_params = dead_layer_model_params( vary_det_deadlayer = 0,
                                        quadratic_source = 0,
                                        quadratic_det = 0,
                                        calibrate_each_pixel = 0,
                                        interp_stopping_power = 1,
                                        mu = 0,
                                        average_over_source = 0,
                                        pulse_height_defect = 0,
                                        fstrips = fstrips,
                                        bstrips = np.arange( 2, 30 ) )
                                        # ignore_outer_pixels = 0 )

                                        

linear_calibration_on_each_x_strip( dbs, source_indices,
                                    model_params,
                                    annotate = 0,
                                    # view_pixel = [ dbmgr.centered, 13 ],
                                    subtitle = subtitle )



# secant_matrices = geom.get_secant_matrices( compute_source_sectheta = 1,
#                                                 average_over_source = 0 )

# print( secant_matrices[ 'cf_249' ][1].x ) 

# secant_matrices = geom.get_secant_matrices( compute_source_sectheta = 1,
#                                                 average_over_source = 1 )
# print( secant_matrices[ 'cf_249' ][1].x ) 



