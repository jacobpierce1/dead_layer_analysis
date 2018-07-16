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






def data_retriever( x, y ) :

    return bpt.data_fetcher( bpt_data_path, 'centered', 1, x, y ) 





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



# spec.plot_hitmaps ( hitmaps )






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
    
    # print( params )
    # print( params.dtype ) 
    
    first_displacement = params[0:3 ]

    # restrict the source thicknesses to be positive by
    # applying positive injective function 
    A = np.exp( params[ 3: ] )
    # A = params[ 3: ] 
    
    integrals = compute_basis_function_integrals( first_displacement, nhat,
                                                  max_radius, mesh_size )

    # print( integrals.dtype ) 


    # print( integrals ) 
    # print( np.max(integrals ) )
    # print( np.min( integrals ) ) 
    
    counts = integrals.dot( A ) 

    # print( counts )
    # print( np.max( counts ) )

    # print( counts.shape ) 
    
    return counts 








def variable_thickness_resid( params, hitmap, d_hitmap, max_radius, mesh_size ) :

    # print( params )
    # print( hitmap )
    # print( d_hitmap )
    # print( max_radius )
    # print( mesh_size ) 
    
    counts = variable_thickness_source_counts( params, hitmap, max_radius, mesh_size ) 

    resid = ( counts - hitmap ) / d_hitmap 
    
    resid[ np.isnan(resid) ] = 0

    # print( hitmap )
    # print( counts ) 
    # print( resid ) 

    ret = resid.flatten()

    # print( ret.dot( ret ) )

    return ret




group_ranges = [ [ -78, 20 ], [-90, 20], [-93, 85] ]


num_peaks_to_detect = 6



cut_strips = 1 


params_guesses = np.array( [ [ [ 3.0e7, 87.10, -9.31, 58.35 ] ],
                             [ [ 3.0e7, 32.95, -31.45, 57.88 ] ],
                             [ [ 3.0e7, 85.35, -46.06, 57.74 ] ] ] )

source_names = [ '240 Pu', '238 Pu', '249 Cf' ]



for name in [ 'centered' ] :

    db = spec.spectrum_db( name, '../../../storage/' ) 

    hitmaps = db.compute_hitmaps( group_ranges,
                                  plot = 0, load = 1 )

    for d in range( db.num_dets ) :
        for i in range( db.num_groups ) :
            filter_data( hitmaps[i][d] )
            if cut_strips : 
                remove_strips( hitmaps[i][d], [0,31,10,12,27,15,3], [0,31,15,16,21] ) 


    db.fit_hitmaps( hitmaps, point_source_counts, params_guesses,
                    source_names = source_names )
    
