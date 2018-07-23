# in this script, we will use integrated hit hit counts from each source
# on each pixel to estimate source distributions from sources placed
# below the detector. given the massive amount of information (32^2 data points )
# for each source, not to mention more data obtained from aggregating datasets,
# the source height distributions can be heavily constrained.

import matplotlib
matplotlib.use('agg')

import numpy as np
import jutils 
import scipy.integrate
import matplotlib.pyplot as plt
import jspectroscopy as spec 
import bpt
import scipy.optimize
from scipy.stats import chi2
import os 
from mpl_toolkits.axes_grid1 import make_axes_locatable
import matplotlib.cm
import colorcet
import sys 

import exp1_geometry

# config
# use_variable_thickness = 1
mesh_size = 4


                   

num_basis_functions = mesh_size

cut_strips = 1

integrals = None 



# a function that is guarenteed to detect the reference peak
# of each group given a list of the peaks detected in the spectrum

def primary_peak_detector( peaks_detected, histo ) :

    if len( peaks_detected ) < 6 :
        return None

    indices = np.empty(3, dtype = int )

    ret = np.empty( 3, dtype = int ) 
    
    # primary peaks 0 and 2 are the largest peaks from the left and right, respectively.

    # we are counting on detection of the first peak. without that
    # all is lost
    
    if abs( peaks_detected[1] - peaks_detected[0] ) > 45 :
        return None 
        
    indices[0] = 1
    
    peak0 = peaks_detected[1]
    ret[0] = peak0
    
    next_peak_offsets = [ 174, 320 ]
    
    for i in range(2) :
        ret[ i+1 ] = ( peak0 + next_peak_offsets[i] - 15
                       + np.argmax( histo[ peak0 + next_peak_offsets[i] - 15 :
                                           peak0 + next_peak_offsets[i] + 15 ] ) )
    
    return ret 






def filter_data( hitmap ) :

    diff = np.diff( hitmap, axis = 0 ) 

    reduced_hitmap = hitmap[ : len( hitmap ) - 1, : ]

    tmp = np.abs( diff / reduced_hitmap ) 
    
    bad_data = ( tmp >= 0.3 ).astype( bool ) 
    
    bad_data = np.append( bad_data, np.zeros((1,32), dtype = bool ), axis = 0 )     

    hitmap[ bad_data ] = np.nan
    
    




def remove_strips( hitmap, fstrips, bstrips ) :

    for f in fstrips :
        hitmap[ f, : ] = np.nan

    for b in bstrips :
        hitmap[ :, b ] = np.nan



group_ranges = [ [ [-60, 20] ], [[-60, 20]], [[-93, 85]] ]


num_peaks_to_detect = 6




exp1_secant_matrices = exp1_geometry.get_secant_matrices()
peak_group_ranges = [ [ -np.inf, 2700 ], [ 2700, 2900 ], [ 2900, np.inf ] ]

thresholds = [ 20, 20, 20 ] 

# for name in [ 'centered', 'moved', 'flat', 'angled' ] : # , 'det3_moved', 'det3_cent' ] : 
for name in [ 'moved' ] :

    db = spec.spectrum_db( name, '../../../storage/' ) 

    
    secant_matrices = exp1_secant_matrices[ name ] 

    # primary_peaks = spec.write_peakdetect( db, primary_peak_detector, data_retriever,
    #                                        num_peaks_to_detect )

    primary_peaks = db.compute_primary_peaks( group_ranges,
                                              primary_peak_detector,
                                              num_peaks_to_detect = num_peaks_to_detect,
                                              primary_peak_detector = primary_peak_detector,
                                              plot = 0, load = 1 )

    # primary_peaks = db.compute_primary_peaks( peak_group_ranges, thresholds,
    #                                           load = 0, plot = 1  )
    
    stds = db.compute_stds( primary_peaks, group_ranges,
                              plot = 0, load = 1 ) 

    means = db.compute_means( primary_peaks, group_ranges,
                              plot = 0, load = 1 ) 

    # averages = spec.compute_averages( db, primary_peaks, data_retriever


    source_names = [ '240 Pu (moved)', '238 Pu (centered)', '249 Cf (moved)' ]


    db.plot_heatmap( 'primary_peaks', source_names )
    db.plot_heatmap( 'stds', source_names )
    db.plot_heatmap( 'means', source_names ) 

    db.plot_vs_sectheta( 'stds', source_names, secant_matrices )
    db.plot_vs_sectheta( 'primary_peaks', source_names, secant_matrices )
    db.plot_vs_sectheta( 'means', source_names, secant_matrices ) 

    db.disconnect()


