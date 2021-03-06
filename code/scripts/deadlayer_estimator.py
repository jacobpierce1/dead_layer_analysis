# This script reads in dataframes containing: 1. the penetration angle
# of each pixel ( constructed in deadlayer_geometry.py ), 2. the mu
# values of each peak fit ( originally detected in parse_all_data.py
# and read into a matrix in deadlayer_analysis.py ) it then reads in a
# file of stopping power data and interpolates it. using all this data
# we make estimates of both the source dead layer depths and
# the pixel dead layer depth.

import matplotlib 
# matplotlib.use('Agg')

import os 
import sys 
import time

import jutils
import jutils.meas as meas 

import jspectroscopy as spec


from scipy.interpolate import interp1d
import scipy.integrate


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

import scipy.integrate



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
                  disable_sources = 0,
                  fix_source_deadlayers = None,
                  one_source_constant = 0,
                  det_thickness = None,
                  vary_source_thickness = 0,
                  source_radii = None,
                  use_different_intercepts = 0,
                  constant_energy_offset = 0 ) :

        self.dimx = 32
        self.dimy = 32 
        
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

        self.disable_sources = disable_sources

        self.one_source_constant = one_source_constant

        self.account_for_det = 0
        self.det_thickness = det_thickness

        self.vary_source_thickness = vary_source_thickness

        # not options yet 

        self.num_source_basis_functions = 8
        self.basis_function_integrals = None
        self.source_radii = source_radii
        
        self.estimate_source_stopping_power = 0
        self.use_different_intercepts = use_different_intercepts
        self.constant_energy_offset = constant_energy_offset 
        
        return None

        



    




        

# multiply by stopping power at emission energy and take dot product
# with (h0, ..., h_N-1 ) to get energy loss due to source distribution
# with height given by P0 finite elements

def compute_basis_function_integrals( displacement, nhat, mesh ) :

    print( 'INFO: computing integrals...' )
    
    R = mesh[-1]

    def normalizer_integrand( r, phi, d ) :
        rprime = np.array( [ r * np.cos( phi ), r * np.sin( phi ), 0 ]  )
        tmp = ( d - rprime ) 
        return r * np.cos( phi ) * ( tmp ).dot( nhat ) / np.linalg.norm( tmp ) ** 3 

    def mesh_integrand( r, phi, d, i ) :
        rprime = np.array( [ r * np.cos( phi ), r * np.sin( phi ), 0 ]  )
        tmp = ( d - rprime ) 
        return ( ( 1 / d[2] ) * r * np.cos( phi ) * ( tmp ).dot( nhat )
                 / np.linalg.norm( tmp ) ** 2 )
        
    normalizers = np.zeros( (32,32) )
    
    ret = np.zeros( (32,32, len(mesh)-1 ) )

    for i in range(32) :
        for j in range(32) :

            # this depends on nhat 
            d = displacement + np.array( 2 * i, 2 * j, 0 ) 

            normalizers[i,j] = scipy.integrate.dblquad(
                normalizer_integrand, (0, 2*np.pi ), 0, R,
                args = (d,) )
                    
            # compute integrals along each interval in the mesh (fencepost ) 
                    
            for k in range( len( mesh ) - 1 ) : 

                ret[i,j,k] = scipy.integrate.dblquad(
                    normalizer_integrand, (0, 2*np.pi ), mesh[k], mesh[k+1],
                    args = (d,) )

    ret /= normalizers

    print( 'Done.' )
    
    return ret 






def compute_all_basis_function_integrals( displacements, nhats, meshes, savepath = None ) :

    if savepath is not None :

        return dill.loads( savepath )

    else :
        num_sources = len( displacements ) 

        ret = [ 0 ] * num_sources

        for i in range( num_sources ) : 

            ret[i] = compute_basis_function_integrals( displacements[i], nhats[i], meshes[i] )
            
        dill.dump( ret, savepath )

        return ret 


                                          


                
def get_stopping_power_interpolation( Z ) :
    pass


                
        
        
def create_empty_data_container( model_params, _meas = 0 ) :
    
    ret =  [ [ [ 0 for j in range( model_params.num_peaks_per_source[i] ) ]
               for i in range( model_params.num_sources ) ]
             for dbnum in range( model_params.num_dbs ) ]
    
    for d in range( model_params.num_dbs ) :
        
        for i in range( model_params.num_sources ) :
            
            for j in range( model_params.num_peaks_per_source[i] ) :

                if _meas :
                    data_type = meas

                else :
                    data_type = np 
                
                # tmp = constructor( ( len( model_params.fstrips[d] ),
                #                      len( model_params.bstrips ) ) )
                                
                tmp = data_type.zeros( ( 32, 32 ) )

                tmp[:] = data_type.nan

                ret[d][i][j] = tmp
                
    return ret 
                                
        


    
def compute_all_energies( params, model_params, channels,
                          det_sectheta, source_sectheta, actual_energies ) :

    e = create_empty_data_container( model_params, 1 )
    edet = create_empty_data_container( model_params, 1 )
    loss = create_empty_data_container( model_params, 0 )
    resid = create_empty_data_container( model_params, 0 ) 
    
    for d in range( model_params.num_dbs ) :
        
        for i in range( model_params.num_sources ) :

            tmp_det_sectheta = det_sectheta[ d ][ i ]
            tmp_source_sectheta = source_sectheta[ d ][ i ]
            
            for j in range( model_params.num_peaks_per_source[i] ) :
                
                # for x in range( len( model_params.fstrips[ d ] ) ) :
                for x in range( len( model_params.fstrips[d] ) ) : 
                
                    fstrip = model_params.fstrips[d][ x ]

                    e_, edet_, loss_ = energy_from_mu_lmfit(
                        params, model_params,
                        channels[ d ][i][j][ fstrip ],
                        tmp_det_sectheta[ fstrip ],
                        tmp_source_sectheta[ fstrip ],
                        actual_energies,
                        d, fstrip, i, j, compute_weight = 1  )
                    
                    e[d][i][j][fstrip] = e_
                    edet[d][i][j][fstrip] = edet_
                    loss[d][i][j][fstrip] = loss_
                    resid[d][i][j][fstrip] = ( actual_energies[i][j] - e_.x ) / e_.dx 
                    
    return e, edet, loss, resid 



# def compute_residuals( params, model_params, channels,
#                        source_geometries, actual_energies ) :

#     pass








# the energy as measured by the detector, energy_det = m * mu + b,
# does not equal the initial energy of the alpha because of several
# possible losses. the rest of the function adds possible losses and
# returns a prediction for the actual initial energy.  params are the
# parameters of the fit which change as the regression is run,
# model_params are constant parameters saying what type of fit we are
# using (e.g. with a quadratic term in sec theta )

def energy_from_mu_lmfit( params, model_params,
                          channels, det_sectheta, source_sectheta, actual_energies,
                          db_num, fstrip, i, j,
                          compute_weight = 0 ) : 

    # fstrip = model_params.fstrips[ db_num ][x]
    
    a = params[ 'a_%d_%d' % ( db_num, fstrip ) ].value

    if not model_params.use_different_intercepts : 
        b = params[ 'b_%d_%d' % ( db_num, fstrip ) ].value
    else :
        b = params[ 'b_%d_%d_%d' % ( db_num, fstrip, i ) ].value

        
    if not model_params.disable_sources :
        
        if model_params.one_source_constant : 
            source_constant = params[ 'source_constant_%d' % (i) ].value

        else :
            source_constant = params[ 'source_constant_%d_%d' % (i,j) ].value

    else :
        source_constant = 0

    
    energy_det = (float(a) * channels ) + b 
    
    
    # in this case, we are using the actual depth of the
    # detector dead layer as a parameter.
    
    if not model_params.vary_det_deadlayer:

        # if model_params.different_fstrip_deadlayers : 
        #     det_deadlayer = params[ 'det_deadlayer_%d' % (fstrip,) ].value.item()

        # elif model_params.different_pixel_deadlayers :
        #     det_deadlayer = np.array( [ params[ 'det_deadlayer_%d_%d' % ( x,y ) ]
        #                                 for y in model_params.bstrips ] )

        det_deadlayer = params[ 'det_deadlayer' ].value #.item()
        
        if model_params.interp_det_stopping_power is not None :
            
            S = model_params.det_stopping_power_interp( actual_energies[i][j]
                                                        - source_constant * source_sectheta )
            
            det_constant = det_deadlayer * S 
            
        else:
            det_constant = det_deadlayer * si_stopping_powers[ i, j ]

    else:
        det_constant = 0

    # print( det_constant ) 


    # if model_params.det_thickness :
        
    #     # det_thickness = params[ 'det_thickness' ].value
    #     det_thickness = params[ 'det_thickness' ]
    #     E1 = actual_energies[i][j] - source_constant * source_sectheta.x
    #     E2 = E1 - det_deadlayer * det_sectheta.x * model_params.det_stopping_power_interp( E1 )
    #     Edet = det_thickness * det_sectheta.x * model_params.det_stopping_power_interp( E2 )
    #     resid = Edet - ( a * channels.x + b ) 

    #     if compute_weight : 
    #         weight = 1 / ( a * channels.dx )
    #         return ( resid, weight )

    #     return Edet

    
    combined_deadlayer_losses = det_constant * det_sectheta
        
        
    if not model_params.vary_source_thickness : 
        source_loss =  source_constant * source_sectheta

    else :

        if not model_params.disable_sources : 
        
            thicknesses = np.array( [ params[ 'source_thickness_%d_%d' % ( i, k ) ]
                                      for k in range(
                                              model_params.num_source_basis_functions ) ] )
            
            source_loss = source_constant * thicknesses.dot(
                model_params.basis_function_integrals ) 

    combined_deadlayer_losses += source_loss

    
        
        
    # if model_params.quadratic_source :
    #     source_constant2 = params[ 'source_constant2_%d_%d' % (i,j) ].value.item()

    #     combined_deadlayer_losses += source_constant2 * (source_sectheta.x ** 2)

    # if model_params.quadratic_det :
    #     det_constant2 = params[ 'det_constant2_%d_%d' % (i,j) ].value.item()

    #     combined_deadlayer_losses += det_constant2 * ( det_sectheta.x - 1 ) ** 2
        

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


    
    if model_params.constant_energy_offset :
        if i > 0 : 
            constant_offset = params[ 'constant_offset_%d' % (i) ]
            energy += constant_offset
            energy_det += constant_offset 

        
    return ( energy, energy_det, combined_deadlayer_losses ) 









def objective( params, model_params, channels,
               det_sectheta, source_sectheta, actual_energies ) :

    e, edet, loss, resid = compute_all_energies( params, model_params,
                                                 channels, det_sectheta, source_sectheta,
                                                 actual_energies )
    
    flattened = np.array( [ resid[d][i][j][x,y]
                            for d in np.arange( model_params.num_dbs )
                            for i in np.arange( model_params.num_sources )
                            for j in np.arange( model_params.num_peaks_per_source[i] )
                            for x in range(32)
                            for y in range(32) ] )
                  # for x in np.arange( model_params.num_fstrips[d] )
                  # for y in np.arange( model_params.num_bstrips ) ]

    # print( np.sum( ~ np.isnan( flattened ) ) )
    # print( 'flattened' ) 
    # print( flattened ) 
    
    # ret = np.nansum( flattened ** 2 )  
    ret = flattened
    # print( ret )
    return ret 











# do a fit of the form E = A * mu + b + s * sec(phi) +
# deadlayer_distance * si_stopping_power * sec(theta)

def estimate_deadlayers( model_params, channels_in, actual_energies,
                         det_sectheta, source_sectheta,
                         source_stopping_power_interps,
                         source_deadlayer_guesses,
                         det_stopping_power_interp, det_deadlayer_guess,
                         calibration_coefs_guess,
                         gen_plots = 0, annotate = 0,
                         strip_coords = None,
                         names = None,
                         figpath = None,
                         savepath = None ) :
    
    num_dbs = len( channels_in )
    num_sources = len( channels_in[0] ) 
    num_peaks_per_source = [ len( channels_in[0][i] ) for i in range( num_sources ) ]
    total_num_peaks = np.sum( num_peaks_per_source )

    model_params.num_dbs = num_dbs
    model_params.num_sources = num_sources
    model_params.num_peaks_per_source  = num_peaks_per_source
    model_params.total_num_peaks = total_num_peaks

    
    channels = [ [ [ channels_in[db][i][j].copy()
                     for j in range( num_peaks_per_source[i] ) ]
                   for i in range( num_sources ) ]
                 for db in range( num_dbs ) ]

        
    for db in range( num_dbs ) :
        for i in range( num_sources ) :
            for j in range( num_peaks_per_source[i] ) :
                mask = np.isnan( det_sectheta[db][i] )
                channels[db][i][j][ mask ] = meas.nan

                    
    print( 'num_dbs: ' + str( num_dbs ) )
    print( 'num_sources: ' + str( num_sources ) )
    print( 'num_peaks_per_source: ' + str( num_peaks_per_source ) ) 

    for db in range( num_dbs ) : 
        for i in range( num_sources ) :
            for j in range( num_peaks_per_source[i] ) : 
                for y in range( 32 ) :
                    if y not in model_params.bstrips : 
                        channels[db][i][j][:,y] = meas.nan

    if model_params.pulse_height_defect :
        raise NotImplemented() 

    
    # prepare a giant lmfit.Parameters() for the fit    
    fit_params = lmfit.Parameters()

    
    # this says whether we will put a lower bound on the following fit
    # parameters.

    # set_lower_bound = model_params.quadratic_source or model_params.quadratic_det
    set_lower_bound = 0

    # if model_params.det_thickness is not None :
    #     fit_params.add( 'det_thickness', value = 1000 )

    # if model_params.det_thickness : 
    #     fit_params.add( 'det_thickness', value = 1000, min = 0 ) 

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
            
                
        # final option: just one thickness for the entire grid (nanometers) 
        else :
            fit_params.add( 'det_deadlayer', value = 100, vary = 1 )#, min=0, max=200 )
            if set_lower_bound :
                fit_params[ 'det_deadlayer' ].min = 0


    if not model_params.disable_sources :
        
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
                                        min = 0, max = 100 )

                    else :
                        fit_params.add( source_name, vary = 0,
                                        value = model_params.fix_source_deadlayers[ source_name ] )


                    if set_lower_bound :
                        fit_params[ 'source_constant_%d_%d' % (i,j) ].min = 0

                    if model_params.quadratic_source :
                        fit_params.add( 'source_constant2_%d_%d' % (i,j), value = 0.0, vary = 1  )

                    if model_params.quadratic_det :
                        fit_params.add( 'det_constant2_%d_%d' % (i,j), value = 0.0, vary = 1  )

            if model_params.vary_source_thickness :

                # add coefficients of the radial basis functions
                for k in range( model_params.num_source_basis_functions ) : 

                    fit_params.add( 'source_thickness_%d_%d' % (i,k), value = 0.5,
                                    min = 0, vary = 1  )

                    if not model_params.estimate_source_stopping_power :
                        fit_params[ 'source_thickness_%d_%d' % (i,k) ].max = 1 

                # self.num_source_basis_functions = 8


        if model_params.vary_source_thickness :

            model_params.basis_function_integrals = compute_basis_function_integrals(
                displacements, nhats, savepath ) 

        if model_params.constant_energy_offset :
            for i in range( 1, num_sources ) :
                fit_params.add( 'constant_offset_%d' % (i), value = 0.0, min = -10, max = 10 )
            

    model_params.fstrips = [0] * num_dbs
    model_params.num_fstrips = [0] * num_dbs 
    model_params.num_bstrips = len( model_params.bstrips ) 

    total_num_data = 0 
    
    for d in range( num_dbs ) :

        frontstrips_to_remove = []

        for x in model_params.fstrips_requested :
            skip = [ 0 for i in range( num_sources ) ] 
            total_strip_data = 0 
            for i in range( num_sources ) :
                num_source_data = 0
                for j in range( num_peaks_per_source[i] ) :
                                        
                    num_source_data += np.sum( ~ np.isnan( channels[d][i][j][x]
                                                    [ model_params.bstrips ].x ) ) 
                    
                    # print( channels[d][i][j][x] )
                    # print( num_source_data ) # print( "num_source_data: ", num_source_data )
                    
                total_strip_data += num_source_data
                
                if num_source_data < 3 :
                    skip[i] = 1 

            if any( skip ) :
                frontstrips_to_remove.append( x )
                continue

            else : 
                total_num_data += total_strip_data 
                

            
            fit_params.add( 'a_' + str(d) + '_%d' % ( x, ),
                            value = calibration_coefs_guess[0] )

            if not model_params.use_different_intercepts : 
                fit_params.add( 'b_' + str(d) + '_%d' % ( x, ),
                                value = calibration_coefs_guess[1] )    

            else :
                for i in range( num_sources ) :
                    if not skip[i] : 
                        fit_params.add( 'b_' + str(d) + '_%d_%d' % ( x, i),
                                        value = calibration_coefs_guess[1] )    
                
                
        # end for x in ...
        model_params.fstrips[ d ] = model_params.fstrips_requested[
            ~ np.in1d( model_params.fstrips_requested,
                       frontstrips_to_remove ) ]

        model_params.num_fstrips[d] = len( model_params.fstrips[d] )
        print( 'INFO: removed frontstrips: ' + str( frontstrips_to_remove ) )
        # print( model_params.fstrips[ db.name ] )

    model_params.num_data_points = ( total_num_peaks
                                     * sum( [ len( model_params.fstrips[d] )
                                              for d in range( num_dbs ) ] ) 
                                     * len( model_params.bstrips ) )
    
    print( 'total_num_data: ', total_num_data ) 
    
    if model_params.interp_det_stopping_power :
        model_params.det_stopping_power_interp  = det_stopping_power_interp 
    
    # this is the k parameter of the lennard model. rather
    # than assuming their value for alphas, we let it be a free
    # parameter and will compare to theirs.
    
    if model_params.pulse_height_defect :
        fit_params.add( 'k', value = 0 )
        

    num_params = ( 2 * np.sum( [ len( model_params.fstrips[d] )
                                 for d in range( model_params.num_dbs ) ] )
                   + np.sum( num_peaks_per_source ) * num_dbs )

    print( 'num_params: ', num_params )

    
    mini = lmfit.Minimizer( objective,
                            fit_params,
                            nan_policy = 'omit',
                            fcn_args = ( model_params, channels,
                                         det_sectheta, source_sectheta,
                                         actual_energies ),
                            # options = { 'maxiter' : int(1e7),
                            #             'xatol' : 1e-6,
                            #             'fatol' : 1e-6 }
                            maxfev = int(1e7),
                            xtol = 1e-12,
                            ftol = 1e-12
                            # max_nfev = int(1e7),
                            # xtol = 1e-12,
                            # ftol = 1e-12
    )

    # print( model_params.fstrips ) 
    # print( channels[0][0][0][1] ) 
    # print( '\n\n' ) 
        
    result = mini.minimize( method = 'leastsq' )
                            
    
    lmfit.report_fit( result.params, show_correl = 0 )

    print( 'reduced chisq: ' + str( result.redchi ) )
    print( 'chisq: ' + str( result.chisqr ) )
    print( 'ndata: ' + str( result.ndata ) )
    print( 'nfree: ' + str( result.nfree ) )
    print( 'pvalue: ' + str( 1 - chi2.cdf( result.chisqr, result.nfree ) ) )
    
    e, edet, loss, resid = compute_all_energies( result.params, model_params,
                                                 channels, det_sectheta, source_sectheta,
                                                 actual_energies )
    
    if savepath is not None :
        plot_vs_sectheta( channels, det_sectheta,
                          actual_energies, result.params, model_params, savepath ) 

    for db in range( num_dbs ) :
        for i in range( num_sources ) :
            for j in range( num_peaks_per_source[i] ) :
                resid[db][i][j] = np.expand_dims( resid[db][i][j], axis = 0 ) 
        resid = spec.dssd_data_container( resid[db]  )
        print( resid.shape )
        heatmap_savepath = savepath + 'cal_heatmap.png'
        resid.plot_heatmap( cmap = colorcet.m_diverging_bkr_55_10_c35,
                            show = 1, savepath = heatmap_savepath ) 


        
    return result
    # ci = lmfit.conf_interval(mini, result)
    # lmfit.printfuncs.report_ci(ci)

    
    
            
def plot_vs_sectheta( channels, secant_matrices, actual_energies,
                      params, model_params, savepath ) : 

    # channel_fit = edet_to_channel_fit( edet, params ) 
    # num_data_per_group = [ len( data[i] ) for i in range( len( data ) ) ]

    num_groups = len( channels[0] )
    num_data_per_group = [ len(t) for t in channels[0] ] 
    max_peaks = max( num_data_per_group ) 

    for db in range( len( channels ) ) :
        for x in model_params.fstrips[db] : 

            no_data = 1 

            f, axarr = plt.subplots( max_peaks, num_groups,
                                             figsize = ( 12,8 ), squeeze = 0 )

            f.subplots_adjust( wspace = 0.5, hspace = 0.5 )

            for i in range( num_groups ) :
                for j in range( max_peaks ) :
                    
                    if j < num_data_per_group[i] :
                        
                        a = params[ 'a_%d_%d' % (db, x ) ]
                        b = params[ 'b_%d_%d' % (db, x ) ]
                        source_constant = params[ 'source_constant_%d_%d'
                                                  % ( i, j ) ]

                        tmp1 = channels[db][i][j][x].x
                        tmp2 = channels[db][i][j][x].dx
                        axarr[j,i].errorbar( secant_matrices[db][i][x],
                                             tmp1, tmp2,
                                             ls = 'none' ) 
                        sec = secant_matrices[db][i][x]
                        min_sectheta = min( sec[ ~ np.isnan( tmp1 ) ] ) 
                        max_sectheta = max( sec[ ~ np.isnan( tmp1 ) ] ) 
                        min_chan = ( ( actual_energies[i][j] - source_constant * min_sectheta ) - b ) / a
                        max_chan = ( ( actual_energies[i][j] - source_constant * max_sectheta ) - b ) / a
                        axarr[j,i].plot( [ min_sectheta, max_sectheta ], [ min_chan, max_chan ] ) 
                        
                        nan_count = np.count_nonzero(
                            np.isnan( secant_matrices[db][i][x] )
                            | np.isnan( tmp1 ) ) 
                        
                        if nan_count < len( tmp1 ) :
                            no_data = 0

                    else :
                        axarr[j,i].axis( 'off' )

            outdir = savepath

            os.makedirs( outdir, exist_ok = 1 )

            if not no_data : 
                plt.savefig( outdir +  '%d.png' % x )

            plt.close( 'all' ) 


 




    

# def plot_results_3d( lmfit_result, source_geometries, mu_matrices,
#                      dbs, source_indices,
#                      model_params, savefig_dir ) :

#     # det_sectheta = source_geometries.det_sectheta
#     # source_sectheta = source_geometries.source_sectheta
    
    
#     f = plt.figure( figsize = ( 10, 8 ) )
    

    
#     axarr_2d = [ f.add_subplot( 3,2, 2*i + 1 ) for i in range(3) ]
#     for i in range(2) :
#         axarr_2d[i].set_xticklabels([])

#     axarr_2d[0].set_title( 'Absolute Calibration of Coupled Peaks' )

#     # ax1 = f.add_subplot(222, projection='3d' )
#     # ax2 = f.add_subplot(224, projection='3d' )

#     ax1 = f.add_subplot(222 )
#     ax1.set_xticklabels([])
#     ax2 = f.add_subplot(224 )

#     ax1.set_title( 'Residuals of Decoupled Peaks' )
    
#     axarr_3d = [ ax1, ax2 ]

#     axarr = [ axarr_2d, axarr_3d ]

#     # add some labels
    
#     # \mathrm instead of \text

        
#     f.suptitle( r'Calibration of Entire Detector Under Model: '
#                 + r'$ \tilde{\chi}^2 = '
#                 + '%.2f' % ( lmfit_result.redchi , ) + '$',
#                 fontsize = 20 ) 
            
    
#     # # create common x and y axis labels
#     # f.add_subplot(121, frameon=False)
#     # plt.tick_params(labelcolor='none', top='off', bottom='off', left='off', right='off')
#     # plt.grid(False)
#     # plt.xlabel( r'$sec ( \theta_\mathrm{det} ) $' )
#     # plt.ylabel( r'$E_\mathrm{det}$ (keV)' )

#     # f.add_subplot(221, frameon=False)
#     # plt.tick_params(labelcolor='none', top='off', bottom='off', left='off', right='off')
#     # plt.grid(False)
#     # plt.xlabel( r'$sec ( \theta_\mathrm{det} ) $' )
#     # plt.ylabel( r'$sec ( \theta_\mathrm{source} ) $' )
    
#     label_fontsize = 16
    
#     f.text(0.04, 0.5, r'$E_\mathrm{det}$ (keV)',
#                 fontsize = 16, va='center', rotation='vertical')

#     f.text(0.5, 0.5, r'$ \sec ( \theta_\mathrm{S} ) $',
#            fontsize = 16, va='center', ha='center', rotation='vertical')

#     # f.text( 0.5, 0.05, r'$ \sec ( \theta_\mathrm{det} ) $',
#     #         fontsize = 16, va = 'center', ha='center' )

#     f.subplots_adjust( hspace = 0.0, wspace = 0.4 )

#     axarr[0][2].set_xlabel( r'$ \sec ( \theta_\mathrm{D} ) $', fontsize = 16 )
    
#     axarr[1][1].set_xlabel( r'$ \sec ( \theta_\mathrm{D} ) $', fontsize = 16 )

    
#     # plt.ylabel( r'$(E - E_\mathrm{cal}) / \Delta( E_\mathrm{cal} )' )

    
#     titles = [ 'Pu 240', 'Pu 238', 'Cf 249' ]
    
#     ax_loc2 = -1

#     for i in range(len( source_indices ) ) :

#         for j in source_indices[i] :

            
#             # get the location to put this source in the plot
#             # ax_loc1 is the first index of the "grid" (0 -> 2D plots, 1 -> 3D plots )
#             # and ax_loc2 is the second index.
#             # i == 1 corresponds to the pu_238 data.

#             if i == 1 :
#                 ax_loc1 = 1
#                 ax_loc2 = j
#             else:
#                 ax_loc1 = 0
#                 ax_loc2 += 1
            
#             ax = axarr[ ax_loc1 ][ ax_loc2 ]

#             ax.set_ylabel( titles[i] )
            

#             if model_params.vary_det_deadlayer:
#                 det_constant = lmfit_result.params[ 'source_constant_%d_%d' % (i,j) ].value #.item()
#                 print( 'effective dl: ' + str ( det_constant / si_stopping_powers[i][j] ) )

 
#             # aggragate data in these arrays so that we can add everything to plots
#             # at the same time. 

#             data_dimensions = ( len(dbs), len( model_params.fstrips_requested ), len( model_params.bstrips ) )

#             det_sectheta = meas.meas.empty( data_dimensions )
#             source_sectheta = meas.meas.empty( data_dimensions ) 
#             Edet = meas.meas.empty( data_dimensions )
#             E_residual = np.empty( data_dimensions )
#             Efit = np.empty( data_dimensions )

                                
#             for db_num in range(len(dbs)) : 

#                 db = dbs[ db_num ] 
                
#                 source = db.sources[i]

#                 for row_num in range( len( model_params.fstrips_requested ) ) :
                    
#                     row = model_params.fstrips_requested[ row_num ]

#                     idx = ( db_num, row_num )
                    
#                     if row not in model_params.fstrips[ db.name ] :
#                         det_sectheta[ idx ] = meas.nan
#                         continue

#                     # x = secant_matrices[source][0][row][ model_params.bstrips ]
#                     # y = secant_matrices[source][1][row][ model_params.bstrips ]
#                     z = mu_matrices[ db.name ][i][j][row][ model_params.bstrips ]

#                     test_id = '_' + db.name + '_%d' % ( row, )

#                     a = lmfit_result.params[ 'a' + test_id ].value.item()
#                     b = lmfit_result.params[ 'b' + test_id ].value.item()

                    
#                     det_sectheta[ idx ] = x
#                     source_sectheta[ idx ] = y

#                     E = a * z + b
#                     if ( 0 in E.dx ) :
#                         print('warning: E.dx == 0' )
#                         print(z)
                    
#                     Edet[ idx ] = E

#                     calibrated_E = energy_from_mu_lmfit( lmfit_result.params,
#                                                          z, x, y,
#                                                          db.name, row, i, j,
#                                                          model_params )

#                     E_residual[ idx ] = peak_energies[i][j] - calibrated_E.x
#                     Efit[ idx ] = peak_energies[i][j] - ( calibrated_E.x - E.x )

                    
#             # add aggragated data to 3d or 2d plot

#             det_sectheta = det_sectheta.flatten()
#             source_sectheta = source_sectheta.flatten()
#             Edet = Edet.flatten()
#             E_residual = E_residual.flatten()
#             Efit = Efit.flatten()
            
#             if ax_loc1 == 1 :
                
#                 normalized_residuals = E_residual / Edet.dx

#                 vmax = np.nanmax( np.abs( normalized_residuals[ ( normalized_residuals < 5 ) ] ) )
#                 # linthresh = np.nanmin( np.abs( normalized_residuals ) ) 
#                 linthresh = 1.0
                
#                 im = ax.scatter( det_sectheta.x, source_sectheta.x,
#                                  c = normalized_residuals,
#                                  s = 5,
#                                  # norm = SymLogNorm( linthresh,
#                                  #                   vmin = -vmax, vmax = vmax ),
#                                  vmax = vmax, vmin = -vmax,
#                                  cmap = colorcet.m_diverging_bkr_55_10_c35 )
#                 plt.colorbar( im, ax = ax )

#             else :
#                 ax.errorbar( det_sectheta.x, Edet.x, xerr = det_sectheta.dx, yerr = Edet.dx,
#                              c='b', zorder = 1,
#                              ls = 'none' )                        

#                 mask = ~ np.isnan( Efit )

#                 ax.plot( det_sectheta.x[ mask ] , Efit[ mask ], c = 'r', zorder = 2 ) 


#     if savefig_dir :
#         plt.savefig( savefig_dir, format='eps', dpi=2000 )
    
#     plt.show()


    
    






    



    
    

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

    




