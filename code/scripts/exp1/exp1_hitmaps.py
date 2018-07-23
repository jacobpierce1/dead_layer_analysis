# in this script, we will use integrated hit hit counts from each source
# on each pixel to estimate source distributions from sources placed
# below the detector. given the massive amount of information (32^2 data points )
# for each source, not to mention more data obtained from aggregating datasets,
# the source height distributions can be heavily constrained.
import matplotlib
# matplotlib.use('agg') 

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



        
def remove_data( hitmap, data_to_remove ) :
    for coords in data_to_remove :
        # print( coords ) 
        hitmap[ coords ] = np.nan






def point_source_counts( params ) :

    
    nhat = np.array( [0,0,1]  )
    
    A = params[0]
    first_displacement = params[1:4]

    counts = np.zeros( (32,32) ) 
    
    for x in range(32) :
        for y in range( 32 ) :
            displacement = 2.0 * np.array( [ -y, x, 0 ] ) + first_displacement 
            costheta = displacement.dot( nhat ) / np.linalg.norm( displacement )
            counts[x,y] = A * np.abs( costheta ) / displacement.dot( displacement )

    return counts 





def point_source_resid( params, hitmap, d_hitmap ) :
    counts = point_source_counts( params ) 
    resid = ( counts - hitmap ) / d_hitmap 
    resid[ np.isnan(resid) ] = 0
    return resid.flatten()
    




# multiply by stopping power at emission energy and take dot product
# with (h0, ..., h_N-1 ) to get energy loss due to source distribution
# with height given by P0 finite elements

def compute_basis_function_integrals( first_displacement, nhat, max_radius, mesh_size,
                                      save_path = None  ) :
    
    # if save_path is not None : 
    #     if os.path.exists( save_path ) :
    #         pass
    if np.allclose( first_displacement, compute_basis_function_integrals.first_displacement_prev,
                    1.0e-3 ) : 
        # print( 'INFO: recycling integrals...' ) 
        return compute_basis_function_integrals.integrals_prev

    compute_basis_function_integrals.first_displacement_prev = first_displacement
    
    # print( 'INFO: computing integrals...' )
    
    mesh = np.linspace( 0.0, max_radius, mesh_size ) 

    num_basis_functions = len( mesh )
    
    # consrtuct supports
    
    supports = np.zeros( ( num_basis_functions, 2 ) )

    # print( num_basis_functions )
    
    for i in range( num_basis_functions  ) :

        # print( i ) 

        if i == 0 :
            supports[i] = np.array( [ mesh[0], mesh[1] ] )

        elif i == num_basis_functions - 1 :
            supports[i] = np.array( [ mesh[-2], mesh[-1] ] )

        else :
            supports[i] = np.array( [ mesh[i-1], mesh[i+1] ] )

    # print( supports ) 


    basis_functions = [ lambda x : P1_basis_function( mesh[i], mesh[1], x )
                        for i in range( num_basis_functions ) ]
    
        
    def integrand( r, phi, d, i ) :
        # print( k ) 
        rprime = np.array( [ r * np.cos( phi ), r * np.sin( phi ), 0.0 ]  )
        tmp = ( d - rprime ) 
        return  ( r * np.abs( basis_functions[i]( r )
                              * ( tmp ).dot( nhat )
                              / np.linalg.norm( tmp ) ** 3.0 )  )        
    
    ret = np.zeros( ( 32, 32, num_basis_functions ) )

    for x in range(32) :
        for y in range(32) :

            # this depends on nhat 
            d = first_displacement + 2 * np.array( [ -y, x, 0 ] ) 
                    
            for k in range( num_basis_functions ) : 
                integral = scipy.integrate.dblquad(
                    integrand, 0.0, 2 * np.pi, * supports[k],
                    args = (d,k),
                    epsabs = 1e-11, epsrel = 1e-11 )

                ret[x,y,k] = integral[0]

    compute_basis_function_integrals.integrals_prev = ret 
    
    return ret 
compute_basis_function_integrals.first_displacement_prev = np.zeros( 3 )  
compute_basis_function_integrals.integrals_prev = None


def P1_basis_function( center, mesh_space, x ) :

    if x > center :
        return ( center + mesh_space - x ) / mesh_space

    else :
        return ( x - ( center - mesh_space ) ) / mesh_space 





def variable_thickness_source_counts( params, hitmap, max_radius, mesh_size ) :

    
    nhat = np.array( [ 0, 0, 1.0]  )
    
    first_displacement = params[0:3 ]

    # restrict the source thicknesses to be positive by
    # applying positive injective function 
    A = np.exp( params[ 3: ] )
    
    integrals = compute_basis_function_integrals( first_displacement, nhat,
                                                  max_radius, mesh_size )

    
    counts = integrals.dot( A ) 


    return counts 








def variable_thickness_resid( params, hitmap, d_hitmap, max_radius, mesh_size ) :

    counts = variable_thickness_source_counts( params, hitmap, max_radius, mesh_size ) 

    resid = ( counts - hitmap ) / d_hitmap 
    
    resid[ np.isnan(resid) ] = 0


    ret = resid.flatten()

    return ret




# group_ranges = [ [ -78, 20 ], [-90, 20], [-93, 85] ]


num_peaks_to_detect = 6



source_names = [ '240 Pu', '238 Pu', '249 Cf' ]


# cut_strips = 1 
# strips_to_remove = [0,31,10,12,27,15,3], [0,31,15,16,21]
# params_guesses = np.array( [ [ [ 3.0e7, 87.10, -9.31, 58.35 ] ],
#                              [ [ 3.0e7, 32.95, -31.45, 57.88 ] ],
#                              [ [ 3.0e7, 85.35, -46.06, 57.74 ] ] ] )
# name = 'centered'


# cut_strips = 1
# strips_to_remove = [ [0,31,3,10,12,15,27], [0,31,15,16,21,23,26,29] ]
# params_guesses = np.array( [ [ [ 3.0e7, 87.10, -9.31, 58.35 ] ],
#                              [ [ 3.0e7, 32.95, -31.45, 57.88 ] ],
#                              [ [ 3.0e7, 85.35, -46.06, 57.74 ] ] ] )
# extra_data_to_remove = [ [14,6] ]
# name = 'flat' 


# cut_strips = 1
# strips_to_remove = [ [0,3,10,12,15,27,31], [0,15,16,17,20,21,25,26,31] ]
# params_guesses = np.array( [ [ [ 3.0e7, 87.10, -9.31, 58.35 ] ],
#                              [ [ 3.0e7, -44.59, -28.58, 57.88 ] ],
#                              [ [ 3.0e7, 85.35, -46.06, 57.74 ] ] ] )
# extra_data_to_remove = []
# name = 'moved'



# cut_strips = 0
# strips_to_remove = [ [0,31,3,10,12,15], [0,31,15,16,21,30] ]
# params_guesses = np.array( [ [ [ 3.0e7, 87.10, -9.31, 58.35 ] ],
#                              [ [ 3.0e7, 32.95, -31.45, 58.72 ] ],
#                              [ [ 3.0e7, 85.35, -46.06, 57.74 ] ] ] )
# extra_data_to_remove = []
# name = 'angled'


# cut_strips = 1
# strips_to_remove = [ [0,31,10,27,12,15], [0,31,25,26,16,15,21,27] ]
# params_guesses = np.array( [ [ [ 3.0e7, 87.10, -9.31, 58.35 ] ],
#                              [ [ 3.0e7, 32.95, -31.45, 58.72 ] ],
#                              [ [ 3.0e7, 85.35, -46.06, 57.74 ] ] ] )
# extra_data_to_remove = [ [3,7], [3,7], [11,6] ]
# name = 'det3_moved'


load = 1
cut_strips = 0
strips_to_remove = [ [0,31,10,27,12,15], [0,31,25,26,16,15,21,27] ]
params_guesses = np.array( [ [ [ 3.0e7, 87.10, -9.31, 58.35 ] ],
                             [ [ 3.0e7, 32.95, -31.45, 58.72 ] ],
                             [ [ 3.0e7, 85.35, -46.06, 57.74 ] ] ] )
extra_data_to_remove = [ [3,7], [3,7], [11,6] ]
name = 'det1_moved'



group_ranges = [ [ [-45, 20] ], [[-60, 20]], [[-60, 30]] ]


num_peaks_to_detect = 6



if cut_strips :
    hitmap_label = None
else :
    hitmap_label = 'no_cuts'



exp1_secant_matrices = exp1_geometry.get_secant_matrices()
peak_group_ranges = [ [ -np.inf, 2700 ], [ 2700, 2900 ], [ 2900, np.inf ] ]

thresholds = [ 20, 20, 20 ] 

# for name in [ 'centered', 'moved', 'flat', 'angled' ] : # , 'det3_moved', 'det3_cent' ] : 


analysis_mgr = spec.dssd_analysis_manager( name, '../../../storage/', (32,32),
                                           [1,1,1] ) 

    

# primary_peaks = spec.write_peakdetect( analysis_mgr, primary_peak_detector, data_retriever,
#                                        num_peaks_to_detect )

primary_peaks = analysis_mgr.compute_primary_peaks( group_ranges,
                                          primary_peak_detector,
                                          num_peaks_to_detect = num_peaks_to_detect,
                                          primary_peak_detector = primary_peak_detector,
                                                    plot = 0, load = 1 )



# primary_peaks = analysis_mgr.compute_primary_peaks( peak_group_ranges, thresholds,
#                                           load = 0, plot = 1  )

stds = analysis_mgr.compute_stds( primary_peaks, group_ranges,
                                  plot = 0, load = load ) 

means = analysis_mgr.compute_means( primary_peaks, group_ranges,
                          plot = 0, load = load ) 

# averages = spec.compute_averages( analysis_mgr, primary_peaks, data_retriever



hitmaps = analysis_mgr.compute_hitmaps( group_ranges, plot = 0, load = 1 )

fitfuncs = [ [point_source_counts] for i in range(3) ] 

print( analysis_mgr.num_dets )

for d in range( analysis_mgr.num_dets ) :
    for i in range( analysis_mgr.num_groups ) :
        # filter_data( hitmaps[i][d] )
        if cut_strips : 
            remove_strips( hitmaps[i][0][d], strips_to_remove[0], strips_to_remove[1] ) 
            remove_data( hitmaps[i][0][d], extra_data_to_remove ) 

# analysis_mgr.fit_hitmaps( hitmaps, fitfuncs, params_guesses,
#                 source_names = source_names, plot = 1, label = hitmap_label )

source_names = [ '240 Pu', '238 Pu', '249 Cf' ]


analysis_mgr.plot_heatmap( 'primary_peaks', source_names, show = 1  )
analysis_mgr.plot_heatmap( 'stds', source_names )
analysis_mgr.plot_heatmap( 'means', source_names ) 


# secant_matrices = exp1_secant_matrices[ name ] 
# analysis_mgr.plot_vs_sectheta( 'stds', source_names, secant_matrices )
# analysis_mgr.plot_vs_sectheta( 'primary_peaks', source_names, secant_matrices )
# analysis_mgr.plot_vs_sectheta( 'means', source_names, secant_matrices ) 




# db = spec.spectrum_db( name, '../../../storage/' ) 

# hitmaps = db.compute_hitmaps( group_ranges,
#                               plot = 1, load = 0 )

# for d in range( db.num_dets ) :
#     for i in range( db.num_groups ) :
#         filter_data( hitmaps[i][d] )
#         if cut_strips : 
#             remove_strips( hitmaps[i][d], strips_to_remove[0], strips_to_remove[1] ) 
#             remove_data( hitmaps[i][d], extra_data_to_remove ) 

# db.fit_hitmaps( hitmaps, point_source_counts, params_guesses,
#                 source_names = source_names, plot = 1 )

# # db.plot_heatmap( 'hitmaps', source_names, 1, show = 1 )


