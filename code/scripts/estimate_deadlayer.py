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

import deadlayer_helpers.stopping_power_interpolation as stop
import deadlayer_helpers.geometry as geom
import deadlayer_helpers.sql_db_manager as db
import deadlayer_helpers.analysis as anal

import matplotlib.pyplot as plt
import pandas as pd
import numpy as np
import scipy.optimize


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



# this function takes an array of observed channels x and converts to
# energies of the alphas as they were originally emitted, accounting
# for loss of energy in both the detector and source dead layers and
# assuming a linear calibration of the detector.
#
# format: p[ si_deadlayer_efficiency, a, b, detector_dl, pu238_dl,
# pu240_dl, cf249_dl ]
#
# m and b are the parameters of a linear fit y =
# mx + b identified_energy is something extracted from the peaks that
# is expected to be proportional to the energy of the alpha as it
# entered the detector, such as the extrapolated peak position or the
# mu value obtained from the fit. i am going to try using both
# separately.
#
# the x array must have one entry for each source.
def calibrated_energy( p, x, stop_power_interp_functions ):

    E_current = np.empty( x.shape )
    
    # first use linear calibration to get the energy of the alpha as
    # recorded by the detector:
    for i in np.arange(len(E_current)):

        # this is the energy of the alpha entering the detector, assumed to be
        # independent of the distance traveled. such a loss could be incorporated
        # later.

        E_current[i] = (
            _energy_entering_detector( p[0], p[1], p[2], p[4],
                                      stop_power_interp_funcs[-1],
                                      x[i], cosine_matrices[i] ) )
        
    # todo: account for distance traveled through air by the alpha
    #
    #
    
    # # next: account for energy lost in the source dead layer.
    # for i in np.arange(len(E_current)):
    #     E_current[i] += 


    return E_current



# map from E_0 to a spectrum feature, e.g. mu value or peak value.
# p[0] = linear coefficient, p[1] = const offset, p[2] = distance of the
# deadlayre in mm
def energy_to_feature( p, E0, energy_loss_func ):
    return p[0]*( E0 - energy_loss_func( E0, p[2] ) ) + p[1] 




# the procedure is complicated a bit by the fact that we measure both charge from
# the detector (assuming 100% collection efficiency) along with charge from the
# dead layer which by hypothesis does not have perfect collection efficiency.
# this function uses Newton's method to obtain what the energy entering the
# dead layer of the detector must have been given an interpolation function
# of the stopping power, the channel number, and guesses for m, b, and the
# collection efficiency of the dead layer.
#
# in short, solve ax + b = det_dl_efficiency * ( E0 - interp(E0) ) + interp(E0)
# for E0.
def _energy_entering_detector( det_dl_efficiency, a, b, stop_power_interp_func, channel, cosine_matrix ):

    # guess for solution: calibration ignoring the dead layer
    E0_guess = a*x + b

    
    func = lambda E0: ( det_dl_efficiency * (E0 - stop_power_interp_func( E0 ) ) +
                        stop_power_interp_func(E0) - ( a * channel + b ) )

    return scipy.optimize.newton( func, E0_guess )




# correct for the energy lost in the source deadlayer. use Newton's method. 
def _energy_leaving_source( E_entering_det, ):

    E0_guess = E_entering_det

    



def _populate_stop_power_interp_funcs( stop_power_interp_funcs, all_objects_deadlayer_ids ):
    for i in range( len( stop_power_interp_funcs ) ):

        current_file = (
            _dir_stopping_power_data + 'alpha_stopping_power_' +
            _deadlayer_ids[i] + '.txt'
        )
        
        stop_power_interp_funcs[i] = (
            stop.stopping_power_interpolation.from_nist_file( current_file,
                                                              _deadlayer_densities[i]
            ).interp
        )



# def _estimate_deadlayer( cosine_matrices, alpha_energies, p0 ):
    


def _main():

    # this takes a while to populate. use debug=1 option when developing.
    cosine_matrices = geom.get_cosine_matrices( debug=0 )

    # read in grid of the mu values 
    mu_grid_center = anal.get_mu_grid_where_valid( 2, get_db_contents(
        sql_db_manager.centered_db ) )

    # anal.get_mu_grid_where_valid( 
    
    # guess 100 nm depth for all dead layers, both detector and source.
    # calculation is 100 nm * ( 1 m / 10^9 nm ) * 1000 mm / nm 

    # deadlayer_depth_guess = [ 100.0 / 1e9 * 1000 ] * len( _all_objects )
    det_deadlayer_depth_guess = [ 100.0 / 1e9 * 1000 ]
    
    # # guess for the fraction of energy that is collected in the detector deadlayer
    # deadlayer_efficiency_guess = 0.5

    # actual energies of the 5 peaks we are looking at:
    alpha_energies = jutils.flatten_list( anal.peak_energies )

    # construct the stopping power interpolation functions
    _populate_stop_power_interp_funcs( _stop_power_interp_funcs,
                                       _deadlayer_ids )

    # guess for all params: A, B, det deadlayer guess,
    p0 = [ 2.0, 50.0, det_deadlayer_depth_guess ]
    #    estimate_deadlayer( cosine_matrices, alpha_energies, p0 )

    current_row_cosines = cosine_matrices[:,:,0,:]


    calibration = lambda _p, _x: calibrated_energy( _p, _x,
                                                    stop_power_interp_functions )

    
    ret = jmath.jacob_least_squares( current_row_cosines.x, mu_values, mu_values_delta,
                                     p0, calibrated_energy )


    

# in principle we should be able to see the expected effect for any source
# because of small angular variation throughout pixels.
# this function makes an estimate using this fact. for each strip we
# look at small angular variations and fit to a line to extract the
# parameter identified as the dead layer distance. then we look at the
# distribution of results obtained using this method.
def example_estimate_for_one_source():

    # this takes a while to populate. use debug=1 option when developing.
    cosine_matrices = geom.get_cosine_matrices( debug=0 )

    # read in grid of the mu values as meas class
    mu_grid_center = anal.get_mu_grid_where_valid( 2, anal.read_all( db.centered_db ) )
    
    # guess 100 nm depth for all dead layers, both detector and source.
    # calculation is 100 nm * ( 1 m / 10^9 nm ) * 1000 mm / nm 

    det_deadlayer_depth_guess = 100.0 / 1e9 * 1000 
    
    # # guess for the fraction of energy that is collected in the detector deadlayer
    # deadlayer_efficiency_guess = 0.5

    # actual energies of the 5 peaks we are looking at:
    alpha_energies = jutils.flatten_list( anal.peak_energies )

    # construct the stopping power interpolation functions
    _populate_stop_power_interp_funcs( _stop_power_interp_funcs,
                                       _deadlayer_ids )

    # guess for all params: A, B, det deadlayer guess,
    p0 = [ 2.0, 50.0, det_deadlayer_depth_guess ]
    #    estimate_deadlayer( cosine_matrices, alpha_energies, p0 )

    current_row_cosines = cosine_matrices[:,:,0,:]


    calibration = lambda _p, _x: calibrated_energy( _p, _x,
                                                    stop_power_interp_functions )

    
    ret = jmath.jacob_least_squares( current_row_cosines.x, mu_values, mu_values_delta,
                                     p0, calibrated_energy )



    

    
    # # make the estimate by ignoring the variations in pixel array, only
    # # considering the strip number.
    # if average:
    #     cosine_matrices = cosine_matrices.mean( axis=3 )
    #     estimate_deadlayers_from_average_strip_position()
        
    # else:
    #     print 'AVERAGE=0 NOT IMPLEMENTED'

    # return 1
    



example_estimate_for_one_source()
