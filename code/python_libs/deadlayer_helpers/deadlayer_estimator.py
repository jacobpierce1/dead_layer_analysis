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

from matplotlib.colors import LinearSegmentedColormap
from matplotlib.colors import SymLogNorm
import matplotlib.pyplot as plt
plt.rcParams["font.family"] = "Times New Roman"
plt.rcParams[ 'mathtext.default' ] = 'regular' 

import pandas as pd
import numpy as np
import scipy.optimize
from scipy.stats import chi2

from sklearn import linear_model
import lmfit

import colorcet




plt.rc('text', usetex=True)





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
                  fstrips_requested = None,
                  bstrips = None,
                  different_fstrip_deadlayers = 0,
                  different_pixel_deadlayers = 1,
                  fix_source_deadlayers = None,
                  one_source_constant = 0 ) :

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
        
        if fstrips_requested is None:
            self.fstrips_requested = np.arange(32)
        else:
            self.fstrips_requested = fstrips_requested

        self.fstrips = dict()
            
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

    if model_params.one_source_constant : 
        source_constant = params[ 'source_constant_%d' % (i) ].value
        
    else :
        source_constant = params[ 'source_constant_%d_%d' % (i,j) ].value
    
    energy_det = a * mu + b 

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
            
        if model_params.si_stopping_power_interpolation is not None :

            # print( source_constant ) 
            
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
    
    resid = np.empty( ( # len(dbs)
                       len( jutils.flatten_list( source_indices) )
                       * sum( [ len( model_params.fstrips[ db.name ] ) for db in dbs ] ) 
                       * len( model_params.bstrips ) , ) ) # np.zeros( ( len(dbs), 3, 2, 32, 32 ) )

    # print( len( resid ) )
    
    # for i
    db_names = [ db.name for db in dbs ]

    # keep track of where the 1D residual array has been
    # filled up to.
    resid_idx = 0
    num_bstrips = len( model_params.bstrips )

    for db in dbs :

        det_sectheta, source_sectheta = [ [ secant_matrices[ source ][k] 
                                            for source in db.sources ]
                                          for k in range(2) ]
        
        for i in range( len( source_indices ) ) :
            for j in source_indices[i] :
                for x in model_params.fstrips[ db.name ] :

                    computed_energy, weight = energy_from_mu_lmfit(
                        params,
                        mu_matrices[ db.name ][i][j][x][ model_params.bstrips ],
                        det_sectheta[i][x][ model_params.bstrips ],
                        source_sectheta[i][x][ model_params.bstrips ],
                        db.name, x, i, j,
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
      
    return resid










# do a fit of the form E = A * mu + b + s * sec(phi) +
# deadlayer_distance * si_stopping_power * sec(theta)

def linear_calibration_on_each_x_strip( dbs,
                                        source_indices,
                                        model_params,
                                        annotate = 0,
                                        cut_high_sectheta = 0,
                                        subtitle = '',
                                        view_pixel = None,
                                        reset_angles = None,
                                        residual_scatter_plot = 0,
                                        plot_3d = 0,
                                        savefig_dir = None ) : 


    
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
                mask = ( ( peak_matrices[db_name][i][j].dx > 2 ) | 
                        ( peak_matrices[db_name][i][j].dx < 0.01 ) )
                peak_matrices[db_name][i][j][ mask ] = meas.nan
                
    if reset_angles is None :
        reset_angles = not model_params.average_over_source
                
    secant_matrices = geom.get_secant_matrices( compute_source_sectheta = 1,
                                                average_over_source = model_params.average_over_source,
                                                reset = reset_angles )

    
    # cut out large sec theta 
    if cut_high_sectheta : 
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
            fit_params.add( 'det_deadlayer', value = 100, vary = 1 ) #, min=0 ) #, min = 0.0 )
            if set_lower_bound :
                fit_params[ 'det_deadlayer' ].min = 0

                
    for i in range( len( source_indices ) ) :

        if model_params.one_source_constant :

            source_name = 'source_constant_%d' % (i,)
            
            if model_params.fix_source_deadlayers is None :
                fit_params.add( source_name, value = 6.0 )

            else :
                fit_params.add( source_name, vary = 0,
                                value = model_params.fix_source_deadlayers[ source_name ] )
                    
        else : 
            for j in source_indices[i] :
            
                source_name = 'source_constant_%d_%d' % (i,j)
                
                if model_params.fix_source_deadlayers is None :
                    fit_params.add( source_name, value = 6.0 )

                else :
                    fit_params.add( source_name, vary = 0,
                                    value = model_params.fix_source_deadlayers[ source_name ] )

                
                if set_lower_bound :
                    fit_params[ 'source_constant_%d_%d' % (i,j) ].min = 0
            
                if model_params.quadratic_source :
                    fit_params.add( 'source_constant2_%d_%d' % (i,j), value = 0.0, vary = 1  )

                if model_params.quadratic_det :
                    fit_params.add( 'det_constant2_%d_%d' % (i,j), value = 0.0, vary = 1  )
                
    
                
    for db in dbs :
        print( db.name )

        frontstrips_to_remove = []

        for x in model_params.fstrips_requested :
                        
            # check if no fits converged for this strip.
            # print(x)
            # print( np.array( [ peak_matrices[ db.name ][t][v][x][ model_params.bstrips ].x
            #                                  for t in range( len( source_indices ) )
            #                                  for v in source_indices[t] ] ) )
            
            num_data = np.count_nonzero( ~ np.isnan( np.array( [ peak_matrices[ db.name ][t][v][x][ model_params.bstrips ].x
                                                                 for t in range( len( source_indices ) )
                                                                 for v in source_indices[t] ] ) ) )

            if num_data < len( jutils.flatten_list( source_indices ) ) + 1 + 2 :
                frontstrips_to_remove.append( x )
                # print( frontstrips_to_remove ) 
                continue
            
            fit_params.add( 'a_' + db.name + '_%d' % ( x, ), value = 1.99 )
            # if model_params.different_fstrip_deadlayers :
            #     fit_params[ 'a_' + db.name + '_%d' % ( x, ) ].min = 1.85
            #     fit_params[ 'a_' + db.name + '_%d' % ( x, ) ].max = 2.20
            
            fit_params.add( 'b_' + db.name + '_%d' % ( x, ), value = -260 )    

        model_params.fstrips[ db.name ] = model_params.fstrips_requested[ ~ np.in1d( model_params.fstrips_requested,
                                                                                     frontstrips_to_remove ) ]

        print( 'INFO: removed frontstrips: ' + str( frontstrips_to_remove ) )
        # print( model_params.fstrips[ db.name ] )

        
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


        
    mini = lmfit.Minimizer( objective,
                            fit_params,
                            nan_policy = 'omit',
                            fcn_args = ( peak_matrices, secant_matrices, peak_energies,
                                         dbs, source_indices, model_params ) )

    result = mini.minimize()
                            
        
        
    # result = lmfit.minimize( objective,
    #                          fit_params,
    #                          args = ( peak_matrices, secant_matrices, peak_energies,
    #                                   dbs, source_indices, model_params ),
    #                          nan_policy = 'omit' )

    
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


    
    






    



    
    

def plot_energy_vs_sectheta( lmfit_result, secant_matrices, mu_matrices,
                             dbs, source_indices,
                             model_params,
                             annotate = 0,
                             subtitle = '',
                             view_pixel = None,
                             residual_scatter_plot = 0,
                             angled_3d_plot = 0,
                             savefig_dir = None ) :

    
    # in this case the angled plots are made in 3d.
    if angled_3d_plot :
        f, axarr = plt.subplots( 3, 2, figsize = ( 10, 8 ) )
        for i in range(2) :
            axarr[3+i].get_xaxis().set_visible( False )
            axarr[3+i].get_yaxis().set_visible( False )

        axarr[3] = f.add_subplot(514, projection='3d' )
        axarr[4] = f.add_subplot(515, projection='3d' )
        

    else :
        f, axarr = plt.subplots( 3, 2, figsize = ( 10, 7.5 ) )
        f.subplots_adjust( hspace = 0.5 ) 


        
    # \mathrm instead of \text
    if view_pixel is None :

        if not residual_scatter_plot : 
        
            f.suptitle( subtitle + r'$ \tilde{\chi}^2 = '
                        + '%.2f' % ( lmfit_result.redchi , ) + '$', fontsize = 20  )

        else :
            f.suptitle( r'Model Residuals vs. $\sec ( \theta_\mathrm{det} ) $ For Each Peak'
                        + '\n' + subtitle + ', ' + r'$ \tilde{\chi}^2 = '
                        + '%.2f' % ( lmfit_result.redchi , ) + '$' )
            
    
    # create common x and y axis labels
    
    f.add_subplot(111, frameon=False)

    # hide tick and tick label of the big axes

    plt.tick_params(labelcolor='none', top='off', bottom='off', left='off', right='off')
    plt.grid(False)
    plt.xlabel( r'$\sec ( \theta_\mathrm{det} ) $', fontsize = 20  )

    if view_pixel is None :
        
        if residual_scatter_plot :
            plt.ylabel( r'$(E - E_\mathrm{cal}) / \Delta( E_\mathrm{cal} )', fontsize = 20 )
            
        else :
            f.text(0.04, 0.5, r'$E_\mathrm{det}$ (keV)',
                   fontsize = 16, va='center', rotation='vertical')

    else : 
        plt.ylabel( r'$\mu$ (channels)' )

        
    titles = [ 'Pu 240', 'Pu 238', 'Cf 249' ]
    
    for i in range(len( source_indices ) ) :

        for j in source_indices[i] :
            
            axarr[i,j].set_title( titles[i] )
            
            # energies = np.empty( ( len(dbs), 32, 32 ) )
            if view_pixel is None :
                
                # if model_params.vary_det_deadlayer:
                #     det_constant = lmfit_result.params[ 'source_constant_%d_%d' % (i,j) ].value #.item()
                    # print( 'effective dl: ' + str ( det_constant / si_stopping_powers[i][j] ) )

                for db in dbs : 
                # for d in range( len( dbs ) ) :

                    source = db.sources[i]

                    # else:
                    #     dl = lmfit_result.params[ 'det_deadlayer' ].value.item()
                    #     det_constant = dl * si_stopping_powers[i][j]

                    # source_constant = lmfit_result.params[ 'source_constant_%d_%d' % (i,j) ].value

                    for row_num in range( len( model_params.fstrips_requested ) ) :
                    
                        row = model_params.fstrips_requested[ row_num ]

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

                        E_residual = peak_energies[i][j] - calibrated_E.x
                        Efit = peak_energies[i][j] - ( calibrated_E.x - E.x )

                        if angled_3d_plot :
                            # add_data_for_angled_3d_plot( 
                            pass
                        
                        else :
                        
                            if not annotate:                                                                 

                                # else do a regular 2D plot, either scatter or errorbar.
                                if residual_scatter_plot :
                                    axarr[i,j].scatter( x.x, E_residual / E.dx, color='b' )

                                else :
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

                            if not residual_scatter_plot : 
                                axarr[i,j].plot( x[ mask ].x, Efit[ mask ], c = 'r', zorder = 2 ) 


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

    
    if savefig_dir :
        plt.savefig( savefig_dir, format='eps', dpi=2000 )

    
    plt.show()

    




