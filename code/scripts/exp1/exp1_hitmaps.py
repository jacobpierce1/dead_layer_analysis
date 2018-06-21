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

# config
# use_variable_thickness = 1
mesh_size = 4


                   

num_basis_functions = mesh_size


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






def point_source_counts( params, hitmap ) :

    
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

    counts = point_source_counts( params, hitmap ) 

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
        # print( ret )
        # return ret 
        
    
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
                
                # if integral[0] == 0 :
                #     print( 'int is 0' )
                #     print( i )
                #     print( j )
                #     print( k ) 
                
                # print( integral[0] ) 
                # print( ret[i,j,k] ) 
                
    # print( 'Done.' )

    compute_basis_function_integrals.integrals_prev = ret 
    
    return ret 
compute_basis_function_integrals.first_displacement_prev = np.zeros( 3 )  
compute_basis_function_integrals.integrals_prev = None


def P1_basis_function( center, mesh_space, x ) :

    if x > center :
        return ( center + mesh_space - x ) / mesh_space

    else :
        return ( x - ( center - mesh_space ) ) / mesh_space 


    
        

# def indicator( x1, x2, x ) :

#     if x >= x1 and x <= x2 :
#         return 1

#     else :
#         return 0 




# evaluate P1 finite element basis of x

# def compute_P1_basis( x1, x2, mesh_size ) :

#     mesh = np.linspace( x1, x2, mesh_size )
    
#     P1_basis = []

#     mesh_space = mesh[1]
#     num_basis_functions = mesh_size + 1
    
#     for i in range( num_basis_functions ) :

#         if i == 0 :
#             P1_basis.append( lambda x : - x / mesh_space if ( x >= x1 and x <= mesh[1] )
#                               else 0 )

#         elif i == N :
#             P1_basis.append( lambda x : ( x - mesh[-1] ) / mesh_space ) * indicator( mesh

#         else :
#             P1_basis.append( lambda x : ( x - mesh[i-1) / mesh_space ) if x < mesh
    
                            
                             
# def P1_finite_element_basis( mesh, i, x ) :

#     if i == 0 :

#         if x > mesh[0] and x < mesh[1] :
#             return ( mesh[1] - x ) / mesh[1]

#         else :
#             return 0

#     elif i == len( mesh ) :

#         if x > mesh[-2] and x < mesh[-1] :
#             return ( x - mesh[-1] ) / mesh[1]

#         else :
#             return 0

#     else :

#         if x > mesh[i-1] and x <= mesh[i] :
#             return ( x - mesh[i-1] ) / mesh[1]

#         elif x > mesh[i] and x < mesh[i+1] :
#             return ( mesh[i+1] - x ) / mesh[1]

#         else :
#             return 0

#     print( 'ERROR: should not be reached.' )
#     raise( IndexError ) 
    




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


bpt_data_path = '../../../bpt-data/extracted_root_tree_data'



hitmaps = spec.compute_hitmaps( (32,32), data_retriever,
                                group_ranges, num_peaks_to_detect, primary_peak_detector,
                                save_path = '../../../storage/hitmaps/centered.dill',
                                plot = 0,
                                reset = 0,
                                filter_data = 1,
                                debug_coords = ( 20, 20 ),
                                rel_plot_bounds = [-100, 120] ) 





# f, axarr = plt.subplots( 3, 3, figsize = ( 10, 8 ) ) 


# params_guesses = [ np.array( [ [ 3.0e7, 87.10, -9.31, 58.35 ],
#                                [ 3.0e7, 32.95, -31.45, 57.88 ],
#                                [ 3.0e7, 85.35, -46.06, 57.74 ] ] ),
#                    np.array( [ [ 87.10, -9.31, 58.35  ] + [-1.0] * num_basis_functions,
#                                [ 32.95, -31.45, 57.88 ] + [-1.0] * num_basis_functions,
#                                [ 85.35, -46.06, 57.74 ] + [-1.0] * num_basis_functions ] ) ]
                   
# eps = [ None, 10.0 ] 


# for i in range( 3 ) :

#     fitfuncs = [ point_source_resid, variable_thickness_resid ]
#     names = [ 'point source calibration', 'variable thickness calibration' ] 

#     filter_data( hitmaps[i] )
#     remove_strips( hitmaps[i], [0,31,10,12,27,15,3], [0,31,15,16,21] ) 
    
#     axarr[0,i].imshow( hitmaps[i] ) 


#     hitmap = hitmaps[i]
#     d_hitmap = np.sqrt( hitmaps[i] )
#     d_hitmap[ d_hitmap == 0 ] = 1

#     args = [ ( hitmap, d_hitmap ),
#              ( hitmap, d_hitmap, 1.0, mesh_size ) ]
    

#     for j in [0] : 
#     # for j in range(2) :

#         print( names[j].upper() )
                                 
#         params_guess = np.copy( params_guesses[j][i] )
        
#         ret = scipy.optimize.leastsq( fitfuncs[j], params_guess,
#                                       full_output = 1, args = args[j] )
        
#         # else :
#         #     ret = scipy.optimize.least_squares( fitfuncs[j], params_guess,
#         #                                         args = args[j], x_scale = [ 50., 50., 50. ]
#         #                                         + [ 15 ] * num_basis_functions )

#         result, cov, info, mesg, success = ret
#         success = success >= 1


#         chisqr = np.sum( info['fvec']**2 )
#         ndata = np.count_nonzero( hitmaps[i][ ( hitmaps[i] != np.nan ) ] ) 
#         nfree =  ndata - len( result ) 
#         redchisqr = chisqr / nfree


#         print( result )
#         print( redchisqr )
#         pvalue = 1 - chi2.cdf( chisqr, nfree )
#         print( 'pvalue: %.2f' % pvalue )
#         print( success )


#         if cov is not None :
#             errors = np.sqrt( np.diag( cov ) * redchisqr )
#             print( errors ) 

#         print( '\n\n' ) 

#         if j == 0 : 
#             resid = ( hitmap - point_source_counts( result, hitmap ) ) / d_hitmap 

#         else:
#             resid = ( hitmap - variable_thickness_source_counts( result, hitmap,
#                                                                  1.0, mesh_size ) ) / d_hitmap 
            
#         im = axarr[ 1+j, i ].imshow( resid )
#         f.colorbar( im, ax = axarr[ 1+j, i ] ) 


cut_strips = 0


savepath = '../../../plots_important/point_source/'

titles = [ '240 Pu (moved)', '238 Pu (centered)', '249 Cf (moved)' ]

if cut_strips :
    title = 'Point Source Fits: With Strip Cuts' 
    savepath += 'with_strip_cuts.eps'
    
else :
    title = 'Point Source Fits: No Strip Cuts'
    savepath += 'no_strip_cuts.eps'



params_guesses = np.array( [ [ 3.0e7, 87.10, -9.31, 58.35 ],
                               [ 3.0e7, 32.95, -31.45, 57.88 ],
                               [ 3.0e7, 85.35, -46.06, 57.74 ] ] )


f, axarr = plt.subplots( 2, 3, figsize = ( 10, 8 ) ) 

f.suptitle( title ) 

f.subplots_adjust( wspace = 0.5 )

for i in range( 3 ) :

    
    filter_data( hitmaps[i] )
        
    if cut_strips : 
        remove_strips( hitmaps[i], [0,31,10,12,27,15,3], [0,31,15,16,21] ) 
    
    # axarr[0,i].imshow( hitmaps[i] )

    masked_array = np.ma.array( hitmaps[i], mask=np.isnan( hitmaps[i] ) )
    cmap = matplotlib.cm.jet
    cmap.set_bad('white',1.)
    im = axarr[ 0, i ].imshow( hitmaps[i], cmap = cmap  )
    divider = make_axes_locatable( axarr[0,i] )
    cax = divider.append_axes("right", size="5%", pad=0.05)
    f.colorbar(im, cax=cax)
    
    # f.colorbar( im, ax = axarr[ 0, i ] ) 

    axarr[0,0].set_ylabel( 'Hit maps' )
    axarr[1,0].set_ylabel( 'Point Source \nFit Residuals' )
    axarr[0,i].set_title( titles[i] ) 

    hitmap = hitmaps[i]
    d_hitmap = np.sqrt( hitmaps[i] )
    d_hitmap[ d_hitmap == 0 ] = 1

    args = ( hitmap, d_hitmap )

    ret = scipy.optimize.leastsq( point_source_resid, params_guesses[i],
                                  full_output = 1, args = args )

    result, cov, info, mesg, success = ret
    success = success >= 1


    chisqr = np.sum( info['fvec']**2 )
    ndata = np.count_nonzero( hitmaps[i][ ~ np.isnan( hitmaps[i] ) ] ) 
    nfree =  ndata - len( result ) 
    redchisqr = chisqr / nfree


    print( 'guess: ' + str( params_guesses[i] ) ) 
    print( 'result: ' + str( result ) )
    print( 'redchisqr: ' + str( redchisqr ) )
    pvalue = 1 - chi2.cdf( chisqr, nfree )
    print( 'pvalue: %.2f' % pvalue )
    print( 'fit converged: ' + str( success ) )


    if cov is not None :
        errors = np.sqrt( np.diag( cov ) * redchisqr )
        print( 'errors: ', errors ) 

    print( '\n\n' ) 

    resid = ( hitmap - point_source_counts( result, hitmap ) ) / d_hitmap 

    im = axarr[ 1, i ].imshow( resid, cmap = cmap )
    divider = make_axes_locatable( axarr[1,i] )
    cax = divider.append_axes("right", size="5%", pad=0.1)
    f.colorbar(im, cax=cax)

    axarr[1,i].set_title( r'$\tilde \chi^2 = %d / %d = %.2f$'
                          % ( np.rint( chisqr ), ndata, redchisqr ) )


plt.savefig( savepath, format='eps', dpi=1000 )
    
plt.show() 
