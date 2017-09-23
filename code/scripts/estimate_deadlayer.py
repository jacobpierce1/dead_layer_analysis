# this script reads in dataframes containing: 1. the penetration angle
# of each pixel ( constructed in deadlayer_geometry.py ), 2. the mu
# values of each peak fit ( originally detected in parse_all_data.py
# and read into a matrix in deadlayer_analysis.py ) it then reads in a
# file of stopping power data and interpolates it. using all this data
# we make estimates of both the source dead layer depths and
# the pixel dead layer depth.


import libjacob.jacob_pyplot as jplt
import libjacob.error as error
import libjacob.jacob_utils as jutils

import deadlayer_helpers.stopping_power_interpolation as stop
import deadlayer_helpers.geometry as geom
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
def calibrated_energy( p, x, stop_power_interp_functions, cosine_matrices ):

    E_current = np.empty( x.shape )
    
    # first use linear calibration to get the energy of the alpha as
    # recorded by the detector:
    for i in np.arange(len(E_current)):

        # this is the energy of the alpha entering the detector, assumed to be
        # independent of the distance traveled. such a loss could be incorporated
        # later.

        E_current[i] = (
            energy_entering_detector( p[0], p[1], p[2],
                                      stop_power_interp_funcs[-1],
                                      x[i], cosine_matrices[i] ) )
        
    # todo: account for distance traveled through air by the alpha
    #
    #
    
    # # next: account for energy lost in the source dead layer.
    # for i in np.arange(len(E_current)):
    #     E_current[i] += 


    return E_current





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

        


def _main( average=0 ):

    # this takes a while to populate. use debug=1 option when developing.
    cosine_matrices = geom.get_cosine_matrices( debug=1 )

    # guess 100 nm depth for all dead layers, both detector and source.
    # calculation is 100 nm * ( 1 m / 10^9 nm ) * 1000 mm / nm 
    deadlayer_depth_guess = [ 100.0 / 1e9 * 1000 ] * len( _all_objects )

    # guess for the fraction of energy that is collected in the detector deadlayer
    deadlayer_efficiency_guess = 0.5

    # actual energies of the 5 peaks we are looking at:
    alpha_energies = jutils.flatten_list( anal.peak_energies )

    # guess for all params
    p0 = [ deadlayer_efficiency_guess, 2.0, 50.0, 100.0, 100.0, 100.0, 100.0 ]


    # construct the stopping power interpolation functions
    _populate_stop_power_interp_funcs( _stop_power_interp_funcs,
                                       _deadlayer_ids )

    # make the estimate by ignoring the variations in pixel array, only
    # considering the strip number.
    if average:
        cosine_matrices = error.emean( cosine_matrices )
        estimate_deadlayers_from_average_strip_position()
        
    else:
        print 'AVERAGE=0 NOT IMPLEMENTED'

    return 1
    

_main( average=1 )
