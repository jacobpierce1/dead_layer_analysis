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






def data_retriever( d, x, y ) :
    return bpt.data_fetcher( bpt_data_path, 'centered', 1, x, y ) 


def raw_data_retriever( d, x, y ) :
    return bpt.raw_data_fetcher( bpt_data_path, 'centered', 1, x, y ) 




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



group_ranges = [ [ -60, 20 ], [-60, 20], [-93, 85] ]


num_peaks_to_detect = 6


bpt_data_path = '../../../bpt-data/extracted_root_tree_data'


db = spec.spectrum_db( 'centered', '../../../storage/' ) 



# primary_peaks = spec.write_peakdetect( db, primary_peak_detector, data_retriever,
#                                        num_peaks_to_detect )

primary_peaks = db.get_primary_peaks( primary_peak_detector, data_retriever,
                                      num_peaks_to_detect, load = 1 )


widths = db.get_widths( primary_peaks, raw_data_retriever, group_ranges,
                            plot = 0, load = 0 ) 

# averages = spec.compute_averages( db, primary_peaks, data_retriever

savepath = db.storage_path + 'heatmaps/widths.png'

exp1_secant_matrices = exp1_geometry.get_secant_matrices()

secant_matrices = [ [ exp1_secant_matrices[ 'pu_240' ][0].x ],
                    [ exp1_secant_matrices[ 'pu_238_centered' ][0].x ],
                    [ exp1_secant_matrices[ 'cf_249' ][0].x ] ]

source_names = [ '240 Pu (moved)', '238 Pu (centered)', '249 Cf (moved)' ]


db.plot_heatmap( 'primary_peaks', source_names, 1 )
db.plot_heatmap( 'widths', source_names, 1 ) 
db.plot_vs_sectheta( 'widths', source_names, secant_matrices, 1 )
db.plot_vs_sectheta( 'primary_peaks', source_names, secant_matrices, 1 )





# if cut_strips :
#     title = 'Point Source Fits: With Strip Cuts' 
#     savepath += 'with_strip_cuts.eps'
    
# else :
#     title = 'Point Source Fits: No Strip Cuts'
#     savepath += 'no_strip_cuts.eps'



# params_guesses = np.array( [ [ 3.0e7, 87.10, -9.31, 58.35 ],
#                                [ 3.0e7, 32.95, -31.45, 57.88 ],
#                                [ 3.0e7, 85.35, -46.06, 57.74 ] ] )


f, axarr = plt.subplots( 1, 3, figsize = ( 10, 8 ) ) 

# f.suptitle( title ) 

f.subplots_adjust( wspace = 0.5 )

# for i in range( 3 ) :
#     cmap = colorcet.m_fire
#     cmap.set_bad('white',1.)
#     im = axarr[ i ].imshow( widths[i,0], cmap = colorcet.m_rainbow )
#     divider = make_axes_locatable( axarr[i] )
#     cax = divider.append_axes("right", size="5%", pad=0.05)
#     f.colorbar(im, cax=cax)



db.disconnect()

    
