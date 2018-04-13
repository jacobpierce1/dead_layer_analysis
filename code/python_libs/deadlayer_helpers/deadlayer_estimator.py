# This script reads in dataframes containing: 1. the penetration angle
# of each pixel ( constructed in deadlayer_geometry.py ), 2. the mu
# values of each peak fit ( originally detected in parse_all_data.py
# and read into a matrix in deadlayer_analysis.py ) it then reads in a
# file of stopping power data and interpolates it. using all this data
# we make estimates of both the source dead layer depths and
# the pixel dead layer depth.

import os 
import sys 
import time

import libjacob.jpyplot as jplt
import libjacob.jmeas as meas
import libjacob.jutils as jutils
import libjacob.jmath as jmath
# import libjacob.jstats as jstats


# import deadlayer_helpers.stopping_power_interpolation as stop
# import deadlayer_helpers.stopping_power_interpolation as stop
# import deadlayer_helpers.sql_db_manager as dbmgr
# import deadlayer_helpers.analysis as anal
# import deadlayer_helpers.data_handler as data

import jspectroscopy as spec


from scipy.interpolate import interp1d

import heapq

from mpldatacursor import datacursor 

from matplotlib.colors import LinearSegmentedColormap
from matplotlib.colors import SymLogNorm
import matplotlib.pyplot as plt

plt.rcParams["font.family"] = "Times New Roman"
plt.rcParams[ 'mathtext.default' ] = 'regular' 

import numpy as np
import scipy.optimize
from scipy.stats import chi2

import lmfit

import colorcet




plt.rc('text', usetex=True)





epsilon_0 = 3.67 / 1000  # keV / electron-hole pair 



# peak_energies = np.array( [ [ 5123.68, 5168.17 ],
#                             [ 5456.3, 5499.03 ],
#                             [ 5759.5, 5813.10,  ] ] )   # all in keV

# flattened_peak_calibration_energies = peak_energies.flatten()



# # these are the stopping powers at the above peak energies.

# si_stopping_powers = np.array( [ [ 6.080E+02, 6.046E+02 ],
#                                  [ 5.832E+02, 5.802E+02 ],
#                                  [ 5.626E+02, 5.591E+02 ] ] )

# si_dioxide_stopping_powers = np.array( [ [ 6.489E+02, 6.452E+02 ],
#                                          [ 6.222E+02, 6.190E+02 ],
#                                          [ 6.001E+02, 5.963E+02   ] ] )




def nth_largest(n, iter):
    return heapq.nlargest(n, iter)[-1]






    


# this object stores all the information required to state the model
# being used to describe the relationship between energy entering
# detector and the secant of penetration angle

class deadlayer_model_params( object ) :

    def __init__( self,
                  vary_det_deadlayer = 0,
                  quadratic_source = 0,
                  quadratic_det = 0,
                  calibrate_each_pixel = 0,
                  interp_det_stopping_power = 1,
                  interp_source_stopping_powers = 0,
                  average_over_source = 1,
                  pulse_height_defect = 0,
                  # ignore_outer_pixels = 0,
                  fstrips_requested = None,
                  bstrips = None,
                  different_fstrip_deadlayers = 0,
                  different_pixel_deadlayers = 0,
                  fix_source_deadlayers = None,
                  one_source_constant = 0,
                  det_thickness = None ) :

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
        self.interp_det_stopping_power = interp_det_stopping_power
        self.interp_source_stopping_powers = interp_source_stopping_powers

        
        self.det_stopping_power_interp = None

        # self.ignore_outer_pixels = ignore_outer_pixels

        # these are the rows that will be considered when
        # performing the optimization. others will be ignored.
        
        if fstrips_requested is None:
            self.fstrips_requested = np.arange(32)
        else:
            self.fstrips_requested = fstrips_requested

        self.fstrips = None # to be set later
            
        if bstrips is None:
            self.bstrips = np.arange(32)
        else:
            self.bstrips = bstrips

        self.different_fstrip_deadlayers = different_fstrip_deadlayers
        
        self.different_pixel_deadlayers = different_pixel_deadlayers

        # dict of source deadlayers is supplied and fixed in the fit
        # so that they are not varied in the fit.
        self.fix_source_deadlayers = fix_source_deadlayers 

        self.one_source_constant = one_source_constant

        self.account_for_det = 0
        self.det_thickness = det_thickness
        
        return None

        



    


class source_geom_data( object ) :

    def __init__( self, det_sectheta, source_sectheta, is_angled,
                  det_sectheta_errors = None, source_sectheta_errors = None ) :

        self.det_sectheta = det_sectheta
        self.source_sectheta = source_sectheta
        self.is_angled = is_angled
        self.det_sectheta_errors = det_sectheta_errors
        self.source_sectheta_errors = source_sectheta_errors


        



    






# the energy as measured by the detector, energy_det = m * mu + b,
# does not equal the initial energy of the alpha because of several
# possible losses. the rest of the function adds possible losses and
# returns a prediction for the actual initial energy.  params are the
# parameters of the fit which change as the regression is run,
# model_params are constant parameters saying what type of fit we are
# using (e.g. with a quadratic term in sec theta )

def energy_from_mu_lmfit( params, model_params,
                          channels, det_sectheta, source_sectheta, actual_energies,
                          db_num, x, i, j,
                          compute_weight = 0 ) :

    # print( 'test' ) 

    a = params[ 'a_%d_%d' % ( db_num, x ) ].value.item()
    b = params[ 'b_%d_%d' % ( db_num, x ) ].value.item()

    if model_params.one_source_constant : 
        source_constant = params[ 'source_constant_%d' % (i) ].value
        
    else :
        source_constant = params[ 'source_constant_%d_%d' % (i,j) ].value
    
    energy_det = a * channels + b 

    # in this case, we are using the actual depth of the
    # detector dead layer as a parameter.
    
    if not model_params.vary_det_deadlayer:

        if model_params.different_fstrip_deadlayers : 
            det_deadlayer = params[ 'det_deadlayer_%d' % (x,) ].value.item()

        elif model_params.different_pixel_deadlayers :
            det_deadlayer = np.array( [ params[ 'det_deadlayer_%d_%d' % ( x,y ) ]
                                        for y in model_params.bstrips ] )

        else :
            det_deadlayer = params[ 'det_deadlayer' ].value #.item()
            
        if model_params.interp_det_stopping_power is not None :

            # print( source_sectheta.x ) 
            # print( source_constant )
            # print( actual_energies[i][j] )
            
            S = model_params.det_stopping_power_interp( actual_energies[i][j]
                                                        - source_constant * source_sectheta.x )
            
            det_constant = det_deadlayer * S 
            
        else:
            det_constant = det_deadlayer * si_stopping_powers[ i, j ]

    else:
        det_constant = 0


    if model_params.det_thickness :
        
        # det_thickness = params[ 'det_thickness' ].value
        det_thickness = fit_params[ 'det_thickness' ]
        E1 = actual_energies[i][j] - source_constant * source_sectheta.x
        E2 = E1 - det_deadlayer * det_sectheta.x * model_params.det_stopping_power_interp( E1 )
        Edet = det_thickness * det_sectheta.x * model_params.det_stopping_power_interp( E2 )
        resid = Edet - ( a * channels.x + b ) 

        if compute_weight : 
            weight = 1 / ( a * channels.dx )
            return ( resid, weight )

        return Edet

        
        
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

        energy = np.empty( len( channels ) )

        integrand = lambda E : 1 / ( epsilon_0 - k * si_stopping_power_interpolation( E ) )
                
        for k in range( len( channels ) ) :

            # print( combined_deadlayer_losses ) 
            
            energy[k] = scipy.integrate.quad( integrand,
                                              model_params.si_stopping_power_interpolation_emin,
                                              actual_energies[i][j] - combined_deadlayer_loss )[0]

        energy *= epsilon_0 


    # the weight is the uncertainty of the total computed energy
    # under the model. under all circumstances, this is dominated
    # by the uncertainty in mu.
        
    if compute_weight : 
        weight = 1 / ( a * channels.dx )
        return ( energy, weight ) 

    return energy 










def objective( params, model_params, channels, source_geometries, actual_energies ) :

    # start_time = time.time()
    
    resid = np.empty( model_params.num_data_points ) 

    # keep track of where the 1D residual array has been
    # filled up to.
    resid_idx = 0
    num_bstrips = model_params.num_bstrips 

    for db_num in range( model_params.num_dbs ) :
        
        for i in range( model_params.num_sources ) :

            det_sectheta = source_geometries[ db_num ][i].det_sectheta
            source_sectheta = source_geometries[ db_num ][i].source_sectheta

            for j in range( model_params.num_peaks_per_source[i] ) :
                
                for x in model_params.fstrips[ db_num ] :

                    computed_energy, weight = energy_from_mu_lmfit(
                        params, model_params,
                        channels[ db_num ][i][j][x][ model_params.bstrips ],
                        det_sectheta[x][ model_params.bstrips ],
                        source_sectheta[x][ model_params.bstrips ],
                        actual_energies,
                        db_num, x, i, j, compute_weight = 1  )

                    
                    if not model_params.det_thickness : 

                        # strange bug: when adding meas and scalar, scalar
                        # must be on the right. must fix asap.
                        residual = - computed_energy + actual_energies[i][j]

                    else :
                        residual = computed_energy
                    
                        
                    # resid[ db_num, i, j, x ] = residual.x * weight
                    resid[ resid_idx : resid_idx + num_bstrips ] = residual * weight
                    resid_idx += num_bstrips 
                    
    # ret = resid.flatten()
    
    # print( 'objective: %f' % ( time.time() - start_time, ) )
      
    return resid










# do a fit of the form E = A * mu + b + s * sec(phi) +
# deadlayer_distance * si_stopping_power * sec(theta)

def estimate_deadlayers( model_params, channels, actual_energies,
                         source_geometries, source_stopping_power_interps, source_deadlayer_guesses,
                         det_stopping_power_interp, det_deadlayer_guess,
                         calibration_coefs_guess,
                         gen_plots = 0, annotate = 0,
                         strip_coords = None,
                         names = None ) :

    num_dbs = len( channels )
    num_sources = len( channels[0] ) 
    num_peaks_per_source = [ len( channels[0][i] ) for i in range( num_sources ) ]
    total_num_peaks = np.sum( num_peaks_per_source )

    model_params.num_dbs = num_dbs
    model_params.num_sources = num_sources
    model_params.num_peaks_per_source  = num_peaks_per_source
    model_params.total_num_peaks = total_num_peaks 

    # print( channels[0][0] ) 
    
    print( 'num_dbs: ' + str( num_dbs ) )
    print( 'num_sources: ' + str( num_sources ) )
    print( 'num_peaks_per_source: ' + str( num_peaks_per_source ) ) 
    
    if strip_coords is not None :
        
        if strip_coords == 'all' :
            for db_num in range( num_dbs ) : 
                for i in model_params.fstrips_requested :

                    print( i )
                    
                    ax = plot_one_strip( model_params, channels, source_geometries, [db_num,i] )

                    figpath = '../../images/current_mu_vs_sectheta/'

                    if names is not None :
                        figpath += names[ db_num ] + '/' 

                    if not os.path.exists( figpath ) :
                        os.mkdir( figpath ) 

                    if names is not None :
                        figpath += names[ db_num ]
                        
                    figpath += str(i)
                    
                    plt.savefig( figpath )
                    plt.clf()
                    plt.close()
                
        else :

            db_num, x = strip_coords

            for i in range( num_sources ) :
                for j in range( num_peaks_per_source[i] ) :
            
                    print( (i,j) )

                    # print( source_geometries[ db_num ][i].det_sectheta[x,0] )
                    
                    for k in np.arange( len( channels[0][0][0][0].x ) ) :
                        print( '%d \t %.1f +- %.1f \t %.2f'
                               % ( k, channels[ db_num ][i][j][x,k].x,
                                   channels[ db_num ][i][j][x,k].dx,
                                   source_geometries[ db_num ][i].det_sectheta[x,k].x  ) )
                    print( '\n\n' )                        

            ax = plot_one_strip( model_params, channels, source_geometries, strip_coords )
            plt.show() 
        
        # plot_energy_vs_sectheta( None, secant_matrices, peak_matrices,
        #                          dbs, source_indices, model_params,
        #                          subtitle = subtitle,
        #                          view_pixel = view_pixel )
        return 1

    
        
    if model_params.pulse_height_defect :
        raise NotImplemented() 

    
    # prepare a giant lmfit.Parameters() for the fit    
    fit_params = lmfit.Parameters()

    
    # this says whether we will put a lower bound on the following fit
    # parameters.

    # set_lower_bound = model_params.quadratic_source or model_params.quadratic_det
    set_lower_bound = 1

    # if model_params.det_thickness is not None :
    #     fit_params.add( 'det_thickness', value = 1000 )

    if model_params.det_thickness : 
        fit_params.add( 'det_thickness', value = 1000, min = 0 ) 

    if not model_params.vary_det_deadlayer :

        # add 1 dead layer thickness for entire 
        if  model_params.different_fstrip_deadlayers : 
            for x in model_params.fstrips_requested :
                fit_params.add( 'det_deadlayer_%d' % (x,), value = 100 ) # , min=70, max=130 )
            
                    
        # add 32^2 different dead layers
        elif model_params.different_pixel_deadlayers :
            for x in model_params.fstrips_requested :
                for y in model_params.bstrips : 
                    fit_params.add( 'det_deadlayer_%d_%d' % (x,y), value = 100 ) # , min=70, max=130 )
            
                
        # final option: just one thickness for the entire grid.
        else :
            fit_params.add( 'det_deadlayer', value = 100, vary = 1, min=0, max=200 ) #, min = 0.0 )
            if set_lower_bound :
                fit_params[ 'det_deadlayer' ].min = 0


                
    for i in range( num_sources ) :

        if model_params.one_source_constant :

            source_name = 'source_constant_%d' % (i,)
            
            if model_params.fix_source_deadlayers is None :
                fit_params.add( source_name, value = source_deadlayer_guesses[i] )

            else :
                fit_params.add( source_name, vary = 0,
                                value = model_params.fix_source_deadlayers[ source_name ] )
                    
        else : 
            for j in range( num_peaks_per_source[i] ) :
            
                source_name = 'source_constant_%d_%d' % (i,j)
                
                if model_params.fix_source_deadlayers is None :
                    fit_params.add( source_name, value = source_deadlayer_guesses[i][j],
                                    min = 0, max = 30 )

                else :
                    fit_params.add( source_name, vary = 0,
                                    value = model_params.fix_source_deadlayers[ source_name ] )

                
                if set_lower_bound :
                    fit_params[ 'source_constant_%d_%d' % (i,j) ].min = 0
            
                if model_params.quadratic_source :
                    fit_params.add( 'source_constant2_%d_%d' % (i,j), value = 0.0, vary = 1  )

                if model_params.quadratic_det :
                    fit_params.add( 'det_constant2_%d_%d' % (i,j), value = 0.0, vary = 1  )
                
    

    model_params.fstrips = [0] * num_dbs
    model_params.num_bstrips = len( model_params.bstrips ) 
                    
    for i in range( num_dbs ) :

        frontstrips_to_remove = []

        for x in model_params.fstrips_requested :
                              
            num_data = np.count_nonzero( ~ np.isnan( np.array( [ channels[i][t][v][x][ model_params.bstrips ].x
                                                                 for t in range( num_sources )
                                                                 for v in range( num_peaks_per_source[t] ) ] ) ) )

            if num_data < total_num_peaks + 1 + 2 :
                frontstrips_to_remove.append( x )
                # print( frontstrips_to_remove ) 
                continue
            
            fit_params.add( 'a_' + str(i) + '_%d' % ( x, ), value = 1.99 )
            
            fit_params.add( 'b_' + str(i) + '_%d' % ( x, ), value = -260 )    

        model_params.fstrips[ i ] = model_params.fstrips_requested[ ~ np.in1d( model_params.fstrips_requested,
                                                                                     frontstrips_to_remove ) ]

        print( 'INFO: removed frontstrips: ' + str( frontstrips_to_remove ) )
        # print( model_params.fstrips[ db.name ] )

    model_params.num_data_points = ( total_num_peaks
                                     * sum( [ len( model_params.fstrips[i] ) for i in range( num_dbs ) ] ) 
                                     * len( model_params.bstrips ) )



    if model_params.interp_det_stopping_power :
        model_params.det_stopping_power_interp  = det_stopping_power_interp 
    
    # this is the k parameter of the lennard model. rather
    # than assuming their value for alphas, we let it be a free
    # parameter and will compare to theirs.
    
    if model_params.pulse_height_defect :
        fit_params.add( 'k', value = 0 )
        
        
    mini = lmfit.Minimizer( objective,
                            fit_params,
                            nan_policy = 'omit',
                            fcn_args = ( model_params, channels, source_geometries, actual_energies ) )

    result = mini.minimize()
                            
    
    lmfit.report_fit( result.params, show_correl = 0 )

    print( 'reduced chisq: ' + str( result.redchi ) )
    print( 'chisq: ' + str( result.chisqr ) )
    print( 'ndata: ' + str( result.ndata ) )
    print( 'nfree: ' + str( result.nfree ) )
    print( 'pvalue: ' + str( 1 - chi2.cdf( result.chisqr, result.nfree ) ) )

    # print( '' )

    # ci = lmfit.conf_interval( mini, result, sigmas = ( 1, 2 ),
    #                           p_names = [ 'det_deadlayer', 'source_constant_0_0', 'source_constant_0_1' ] )

    # lmfit.printfuncs.report_ci( ci ) 
    

    # determine whether to use 3d plot or 2d plot
    
    angled_3d_plot = ( plot_3d and ( not model_params.vary_det_deadlayer )
                       and dbmgr.angled in dbs )


    if angled_3d_plot :
        plot_results_3d( result, secant_matrices, peak_matrices,
                         dbs, source_indices,
                         model_params,
                         savefig_dir )

        
    else :
        plot_energy_vs_sectheta( result, secant_matrices, peak_matrices,
                                 dbs, source_indices,
                                 model_params,
                                 annotate = annotate,
                                 subtitle = subtitle,
                                 view_pixel = view_pixel,
                                 residual_scatter_plot = residual_scatter_plot,
                                 angled_3d_plot = angled_3d_plot,
                                 savefig_dir = savefig_dir ) 
    






def plot_one_strip( model_params, channels, source_geometries, coords ) :


    
    db_num, x = coords

    maxpeaks = max( model_params.num_peaks_per_source ) 

    f, axarr = plt.subplots( model_params.num_sources, maxpeaks, figsize = (10,8) ) 

    for i in range( model_params.num_sources ) :
        for j in range( model_params.num_peaks_per_source[i] ) :

            axarr[i,j].errorbar( source_geometries[ db_num ][i].source_sectheta.x[x],
                                 channels[ db_num ][i][j][x].x,
                                 yerr = channels[ db_num ][i][j][x].dx,
                                 ls = 'none' )

    # plt.show() 
    return axarr

        





    

def plot_results_3d( lmfit_result, secant_matrices, mu_matrices,
                     dbs, source_indices,
                     model_params, savefig_dir ) :

    
    f = plt.figure( figsize = ( 10, 8 ) )
    

    
    axarr_2d = [ f.add_subplot( 3,2, 2*i + 1 ) for i in range(3) ]
    for i in range(2) :
        axarr_2d[i].set_xticklabels([])

    axarr_2d[0].set_title( 'Absolute Calibration of Coupled Peaks' )

    # ax1 = f.add_subplot(222, projection='3d' )
    # ax2 = f.add_subplot(224, projection='3d' )

    ax1 = f.add_subplot(222 )
    ax1.set_xticklabels([])
    ax2 = f.add_subplot(224 )

    ax1.set_title( 'Residuals of Decoupled Peaks' )
    
    axarr_3d = [ ax1, ax2 ]

    axarr = [ axarr_2d, axarr_3d ]

    # add some labels
    
    # \mathrm instead of \text

        
    f.suptitle( r'Calibration of Entire Detector Under Model: '
                + r'$ \tilde{\chi}^2 = '
                + '%.2f' % ( lmfit_result.redchi , ) + '$',
                fontsize = 20 ) 
            
    
    # # create common x and y axis labels
    # f.add_subplot(121, frameon=False)
    # plt.tick_params(labelcolor='none', top='off', bottom='off', left='off', right='off')
    # plt.grid(False)
    # plt.xlabel( r'$sec ( \theta_\mathrm{det} ) $' )
    # plt.ylabel( r'$E_\mathrm{det}$ (keV)' )

    # f.add_subplot(221, frameon=False)
    # plt.tick_params(labelcolor='none', top='off', bottom='off', left='off', right='off')
    # plt.grid(False)
    # plt.xlabel( r'$sec ( \theta_\mathrm{det} ) $' )
    # plt.ylabel( r'$sec ( \theta_\mathrm{source} ) $' )
    
    label_fontsize = 16
    
    f.text(0.04, 0.5, r'$E_\mathrm{det}$ (keV)',
                fontsize = 16, va='center', rotation='vertical')

    f.text(0.5, 0.5, r'$ \sec ( \theta_\mathrm{S} ) $',
           fontsize = 16, va='center', ha='center', rotation='vertical')

    # f.text( 0.5, 0.05, r'$ \sec ( \theta_\mathrm{det} ) $',
    #         fontsize = 16, va = 'center', ha='center' )

    f.subplots_adjust( hspace = 0.0, wspace = 0.4 )

    axarr[0][2].set_xlabel( r'$ \sec ( \theta_\mathrm{D} ) $', fontsize = 16 )
    
    axarr[1][1].set_xlabel( r'$ \sec ( \theta_\mathrm{D} ) $', fontsize = 16 )

    
    # plt.ylabel( r'$(E - E_\mathrm{cal}) / \Delta( E_\mathrm{cal} )' )

    
    titles = [ 'Pu 240', 'Pu 238', 'Cf 249' ]
    
    ax_loc2 = -1

    for i in range(len( source_indices ) ) :

        for j in source_indices[i] :

            
            # get the location to put this source in the plot
            # ax_loc1 is the first index of the "grid" (0 -> 2D plots, 1 -> 3D plots )
            # and ax_loc2 is the second index.
            # i == 1 corresponds to the pu_238 data.

            if i == 1 :
                ax_loc1 = 1
                ax_loc2 = j
            else:
                ax_loc1 = 0
                ax_loc2 += 1
            
            ax = axarr[ ax_loc1 ][ ax_loc2 ]

            ax.set_ylabel( titles[i] )
            

            if model_params.vary_det_deadlayer:
                det_constant = lmfit_result.params[ 'source_constant_%d_%d' % (i,j) ].value #.item()
                print( 'effective dl: ' + str ( det_constant / si_stopping_powers[i][j] ) )

 
            # aggragate data in these arrays so that we can add everything to plots
            # at the same time. 

            data_dimensions = ( len(dbs), len( model_params.fstrips_requested ), len( model_params.bstrips ) )

            det_sectheta = meas.meas.empty( data_dimensions )
            source_sectheta = meas.meas.empty( data_dimensions ) 
            Edet = meas.meas.empty( data_dimensions )
            E_residual = np.empty( data_dimensions )
            Efit = np.empty( data_dimensions )

                                
            for db_num in range(len(dbs)) : 

                db = dbs[ db_num ] 
                
                source = db.sources[i]

                for row_num in range( len( model_params.fstrips_requested ) ) :
                    
                    row = model_params.fstrips_requested[ row_num ]

                    idx = ( db_num, row_num )
                    
                    if row not in model_params.fstrips[ db.name ] :
                        det_sectheta[ idx ] = meas.nan
                        continue

                    x = secant_matrices[source][0][row][ model_params.bstrips ]
                    y = secant_matrices[source][1][row][ model_params.bstrips ]
                    z = mu_matrices[ db.name ][i][j][row][ model_params.bstrips ]

                    test_id = '_' + db.name + '_%d' % ( row, )

                    a = lmfit_result.params[ 'a' + test_id ].value.item()
                    b = lmfit_result.params[ 'b' + test_id ].value.item()

                    
                    det_sectheta[ idx ] = x
                    source_sectheta[ idx ] = y

                    E = a * z + b
                    if ( 0 in E.dx ) :
                        print('warning: E.dx == 0' )
                        print(z)
                    
                    Edet[ idx ] = E

                    calibrated_E = energy_from_mu_lmfit( lmfit_result.params,
                                                         z, x, y,
                                                         db.name, row, i, j,
                                                         model_params )

                    E_residual[ idx ] = peak_energies[i][j] - calibrated_E.x
                    Efit[ idx ] = peak_energies[i][j] - ( calibrated_E.x - E.x )

                    
            # add aggragated data to 3d or 2d plot

            det_sectheta = det_sectheta.flatten()
            source_sectheta = source_sectheta.flatten()
            Edet = Edet.flatten()
            E_residual = E_residual.flatten()
            Efit = Efit.flatten()
            
            if ax_loc1 == 1 :
                
                normalized_residuals = E_residual / Edet.dx

                vmax = np.nanmax( np.abs( normalized_residuals[ ( normalized_residuals < 5 ) ] ) )
                # linthresh = np.nanmin( np.abs( normalized_residuals ) ) 
                linthresh = 1.0
                
                im = ax.scatter( det_sectheta.x, source_sectheta.x,
                                 c = normalized_residuals,
                                 s = 5,
                                 # norm = SymLogNorm( linthresh,
                                 #                   vmin = -vmax, vmax = vmax ),
                                 vmax = vmax, vmin = -vmax,
                                 cmap = colorcet.m_diverging_bkr_55_10_c35 )
                plt.colorbar( im, ax = ax )

            else :
                ax.errorbar( det_sectheta.x, Edet.x, xerr = det_sectheta.dx, yerr = Edet.dx,
                             c='b', zorder = 1,
                             ls = 'none' )                        

                mask = ~ np.isnan( Efit )

                ax.plot( det_sectheta.x[ mask ] , Efit[ mask ], c = 'r', zorder = 2 ) 


    if savefig_dir :
        plt.savefig( savefig_dir, format='eps', dpi=2000 )
    
    plt.show()


    
    






    



    
    

# def plot_energy_vs_sectheta( lmfit_result, secant_matrices, mu_matrices,
#                              dbs, source_indices,
#                              model_params,
#                              annotate = 0,
#                              subtitle = '',
#                              view_pixel = None,
#                              residual_scatter_plot = 0,
#                              angled_3d_plot = 0,
#                              savefig_dir = None ) :

    
#     # in this case the angled plots are made in 3d.
#     if angled_3d_plot :
#         f, axarr = plt.subplots( 3, 2, figsize = ( 10, 8 ) )
#         for i in range(2) :
#             axarr[3+i].get_xaxis().set_visible( False )
#             axarr[3+i].get_yaxis().set_visible( False )

#         axarr[3] = f.add_subplot(514, projection='3d' )
#         axarr[4] = f.add_subplot(515, projection='3d' )
        

#     else :
#         f, axarr = plt.subplots( 3, 2, figsize = ( 10, 7.5 ) )
#         f.subplots_adjust( hspace = 0.5 ) 


        
#     # \mathrm instead of \text
#     if view_pixel is None :

#         if not residual_scatter_plot : 
        
#             f.suptitle( subtitle + r'$ \tilde{\chi}^2 = '
#                         + '%.2f' % ( lmfit_result.redchi , ) + '$', fontsize = 20  )

#         else :
#             f.suptitle( r'Model Residuals vs. $\sec ( \theta_\mathrm{det} ) $ For Each Peak'
#                         + '\n' + subtitle + ', ' + r'$ \tilde{\chi}^2 = '
#                         + '%.2f' % ( lmfit_result.redchi , ) + '$' )
            
    
#     # create common x and y axis labels
    
#     f.add_subplot(111, frameon=False)

#     # hide tick and tick label of the big axes

#     plt.tick_params(labelcolor='none', top='off', bottom='off', left='off', right='off')
#     plt.grid(False)
#     plt.xlabel( r'$\sec ( \theta_\mathrm{det} ) $', fontsize = 20  )

#     if view_pixel is None :
        
#         if residual_scatter_plot :
#             plt.ylabel( r'$(E - E_\mathrm{cal}) / \Delta( E_\mathrm{cal} )', fontsize = 20 )
            
#         else :
#             f.text(0.04, 0.5, r'$E_\mathrm{det}$ (keV)',
#                    fontsize = 16, va='center', rotation='vertical')

#     else : 
#         plt.ylabel( r'$\mu$ (channels)' )

        
#     titles = [ 'Pu 240', 'Pu 238', 'Cf 249' ]
    
#     for i in range(len( source_indices ) ) :

#         for j in source_indices[i] :
            
#             axarr[i,j].set_title( titles[i] )
            
#             # energies = np.empty( ( len(dbs), 32, 32 ) )
#             if view_pixel is None :
                
#                 # if model_params.vary_det_deadlayer:
#                 #     det_constant = lmfit_result.params[ 'source_constant_%d_%d' % (i,j) ].value #.item()
#                     # print( 'effective dl: ' + str ( det_constant / si_stopping_powers[i][j] ) )

#                 for db in dbs : 
#                 # for d in range( len( dbs ) ) :

#                     source = db.sources[i]

#                     # else:
#                     #     dl = lmfit_result.params[ 'det_deadlayer' ].value.item()
#                     #     det_constant = dl * si_stopping_powers[i][j]

#                     # source_constant = lmfit_result.params[ 'source_constant_%d_%d' % (i,j) ].value

#                     for row_num in range( len( model_params.fstrips_requested ) ) :
                    
#                         row = model_params.fstrips_requested[ row_num ]

#                         x = secant_matrices[source][0][row][ model_params.bstrips ]
#                         y = secant_matrices[source][1][row][ model_params.bstrips ]
#                         z = mu_matrices[ db.name ][i][j][row][ model_params.bstrips ]

#                         test_id = '_' + db.name + '_%d' % ( row, )

#                         a = lmfit_result.params[ 'a' + test_id ].value.item()
#                         b = lmfit_result.params[ 'b' + test_id ].value.item()

#                         E = a * z + b

#                         # energies[ d, row, : ] = E.x
                            
#                         calibrated_E = energy_from_mu_lmfit( lmfit_result.params,
#                                                              z, x, y,
#                                                              db.name, row, i, j,
#                                                              model_params )

#                         E_residual = peak_energies[i][j] - calibrated_E.x
#                         Efit = peak_energies[i][j] - ( calibrated_E.x - E.x )

#                         if angled_3d_plot :
#                             # add_data_for_angled_3d_plot( 
#                             pass
                        
#                         else :
                        
#                             if not annotate:                                                                 

#                                 # else do a regular 2D plot, either scatter or errorbar.
#                                 if residual_scatter_plot :
#                                     axarr[i,j].scatter( x.x, E_residual / E.dx, color='b' )

#                                 else :
#                                     axarr[i,j].errorbar( x.x, E.x, xerr = x.dx, yerr = E.dx,
#                                                          color='b', zorder = 1,
#                                                          ls = 'none' )

#                             else:
#                                 for a in range( len( model_params.bstrips ) ) :

#                                     if not meas.isnan( z[a] ) :

#                                         label = '(%s,%d,%d,%.1f)' % (db.name, row, model_params.bstrips[a],
#                                                                    z[a].x ) 

#                                         axarr[i,j].scatter( [ x.x[a] ] , [ E.x[a] ], color='b',
#                                                             zorder = 1, label = label )


#                             mask = ~ meas.isnan( z )

#                             if not residual_scatter_plot : 
#                                 axarr[i,j].plot( x[ mask ].x, Efit[ mask ], c = 'r', zorder = 2 ) 


#             else :

#                 db = view_pixel[0]
#                 row = view_pixel[1]

#                 if i == 0 :
#                     source = 'pu_240'

#                 elif i == 1 :
#                     source = 'pu_238_' + db.name
                        
#                 elif i == 2 :
#                     source = 'cf_249'
                
#                 x = secant_matrices[source][0][row][ model_params.bstrips ]
#                 y = secant_matrices[source][1][row][ model_params.bstrips ]
#                 z = mu_matrices[ db.name ][i][j][row][ model_params.bstrips ]

#                 print( (i,j) )
#                 print( 'mu / peak values: ' )
#                 print( z )
#                 print('\n' )
                
#                 axarr[i,j].errorbar( x.x, z.x, yerr = z.dx,
#                                      color='b', zorder = 1, ls = 'none' ) #, label = 'test' )

#     if annotate :
#         datacursor( formatter = '{label}'.format )

    
#     if savefig_dir :
#         plt.savefig( savefig_dir, format='eps', dpi=2000 )

    
#     plt.show()

    




