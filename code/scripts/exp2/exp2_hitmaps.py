# in this script, we will use integrated hit hit counts from each source
# on each pixel to estimate source distributions from sources placed
# below the detector. given the massive amount of information (32^2 data points )
# for each source, not to mention more data obtained from aggregating datasets,
# the source height distributions can be heavily constrained.


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

# config
# use_variable_thickness = 1
mesh_size = 4

                   
num_basis_functions = mesh_size

cut_strips = 1

integrals = None 

debug = 0

detnum = 0

z_displacement = 4.155 * 25.4



def primary_peak_detector( peaks_detected, histo ) :
        
    if len( peaks_detected ) < 2 :
        return None

    if( debug ) :
        print( 'peaks_detected', peaks_detected )

    if ( peaks_detected[0] > 2000
         or peaks_detected[1] < 2000
         or peaks_detected[0] < 1000 ) :

        if debug :
            print( "peak detect failed" )
 
        return None
    
    return peaks_detected 






def data_retriever( x, y ) :

    return bpt.data_fetcher( bpt_data_path, 'full_bkgd_tot', 4, x, y ) 





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



def remove_data( hitmap, data_to_remove ) :

    for coords in data_to_remove :
        # print( 'removing coords: ', coords ) 
        # hitmap[ coords[0], coords[1] ] = np.nan
        hitmap[ coords ] = np.nan



def fitfunc( params, cuts ) :
    
    nhat = np.array( [0,0,1]  )

    A = np.zeros(2)
    first_displacements = np.zeros((2,3))
    
    A[0] = params[0]
    first_displacements[0,0:2] = params[1:3]
    first_displacements[0,2] = params[6]
    
    A[1] = params[3]
    first_displacements[1,0:2] = params[4:6]
    first_displacements[1,2] = params[6]
    
    counts = np.zeros( (32,32) )
    
    for x in range(32) :
        for y in range( 32 ) :
            for k in range(2) :

                if ( k==0 and x >= cuts[k] ) or ( k==1 and x <= cuts[k] ) :
                    
                    displacement = 2.0 * np.array( [ y, x, 0 ] ) + first_displacements[k] 
                    costheta = displacement.dot( nhat ) / np.linalg.norm( displacement )
                    counts[x,y] += ( A[k] * np.abs( costheta )
                                     / displacement.dot( displacement ) )

    return counts 





def point_source_resid( params, hitmap, d_hitmap, cuts ) :
    counts = point_source_counts( params, cuts ) 
    resid = ( counts - hitmap ) / d_hitmap 
    resid[ np.isnan(resid) ] = 0
    return resid.flatten()
    



group_ranges = [ [ -140, 32 ], [-110, 50] ]


num_peaks_to_detect = 2


bpt_data_path = '../../../bpt-data/extracted_root_tree_data'


source_names = [ 'Gd 148', 'Cm 244' ]





# extra_strips_to_remove = [ [ [6,18,17,19,20], [7] ],
#                            [ [1,9,10], [7,8] ] ]




cut_strips = 1


params_guesses = np.array( [ [ [ 3.0e7, -10.0, -30.0, 
                                 3.0e7, -5.0, -30.0,
                                 100.0 ] ] * 4,
                             [ [ 3.0e7, -10.0, -30.0,
                                 3.0e7, -5.0, -30.0,
                                 100.0 ] ] * 4 ] )

strips_to_remove = [
    [
        [ [1,9,10,24,25], [8,9,16,31] ],
        [ [2,10,11,13,14,19,21,25,26,31], [8,16,31] ],
        [ [6,14,17,18,19,20,31], [7,31] ],
        [ [9,16,17,18,19,31], [7,8,16,31] ],
    ],
    [
        [ [4,29,10,16,17], [8,16,31] ],
        [ [5,6,7,8,21,22,30,31], [8,12,16,31] ],
        [ [1,9,10,13,31], [7,8,31] ],
        [ [1,5,6,23,24,29,30], [7,8,16,31] ]
    ]
] 



data_to_remove = [
    [
        [ [] ],
        [ [15,24] ],
        [ [] ],
        [ [30,9] ]
    ],
    [
        [ [] ],
        [ [24,15] ],
        [ [] ],
        [ [] ]
    ]
]


cuts = [
    [
        [10,23],
        [11,25],
        [6,20],
        [10,18]
    ],
    [
        [5,16],
        [10,21],
        [9,23],
        [6,23]
    ]
]
                         


db = spec.spectrum_db( 'full_bkgd_tot', '../../../storage/' ) 

primary_peaks = db.compute_primary_peaks( primary_peak_detector, num_peaks_to_detect,
                                          load = 1 )

hitmaps = db.compute_hitmaps( group_ranges, plot = 1, load = 1 )

print( hitmaps.shape ) 

for i in range( db.num_groups ) :
    for d in range( db.num_dets ) :
        filter_data( hitmaps[i][d] )
        if cut_strips :
            remove_strips( hitmaps[i][d], strips_to_remove[i][d][0],
                           strips_to_remove[i][d][1] ) 
            remove_data( hitmaps[i][d], data_to_remove[i][d] )
                

db.fit_hitmaps( hitmaps, fitfunc, params_guesses,
                source_names = source_names, plot = 0,
                extra_args = cuts )# , label = 'no_cut'  )

group_ranges = [ [-30,20], [-50,20] ]

widths = db.compute_widths( primary_peaks, group_ranges,
                            plot = 0, load = 0 ) 

db.plot_heatmap( 'widths', source_names, 1 )


db.disconnect()
