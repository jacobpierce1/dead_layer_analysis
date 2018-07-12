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
        hitmap[ coords[0], coords[1] ] = np.nan




def point_source_counts( params ) :

    
    nhat = np.array( [0,0,1]  )

    A = np.zeros(2)
    first_displacements = np.zeros((2,3))
    cuts = np.zeros(2)
    
    A[0] = params[0]
    first_displacements[0,0:2] = params[1:3]
    first_displacements[0,2] = params[8]
    cuts[0] = params[3]
    
    A[1] = params[4]
    first_displacements[1,0:2] = params[5:7]
    first_displacements[1,2] = params[8]
    cuts[1] = params[7]
    
    counts = np.zeros( (32,32) ) 
    
    for x in range(32) :
        for y in range( 32 ) :
            for k in range(2) :

                if ( k==0 and x > cuts[k] ) or ( k==1 and x < cuts[k] ) :
                    displacement = 2.0 * np.array( [ y, x, 0 ] ) + first_displacements[k] 
                    costheta = displacement.dot( nhat ) / np.linalg.norm( displacement )
                    counts[x,y] += ( A[k] * np.abs( costheta )
                                     / displacement.dot( displacement ) )

    return counts 





def point_source_resid( params, hitmap, d_hitmap ) :
    counts = point_source_counts( params ) 
    resid = ( counts - hitmap ) / d_hitmap 
    resid[ np.isnan(resid) ] = 0
    return resid.flatten()
    



group_ranges = [ [ -140, 32 ], [-110, 50] ]


num_peaks_to_detect = 2


bpt_data_path = '../../../bpt-data/extracted_root_tree_data'



hitmaps = spec.compute_hitmaps( (32,32), data_retriever,
                                group_ranges, num_peaks_to_detect, primary_peak_detector,
                                save_path = '../../../storage/full_bkgd_tot/tmp/hitmaps.dill',
                                plot = 1,
                                reset = 0,
                                filter_data = 1,
                                debug_coords = None, #( 20, 20 ),
                                rel_plot_bounds = [-200, 120] ) 




savepath = '../../../plots_important/point_source/'

titles = [ 'Gd 148', 'Cm 244' ]

if cut_strips :
    title = 'Point Source Fits: With Strip Cuts' 
    savepath += 'with_strip_cuts.eps'
    
else :
    title = 'Point Source Fits: No Strip Cuts'
    savepath += 'no_strip_cuts.eps'



# params_guesses = np.array( [ [ 3.0e7, -10.0, -10.0, 100.0, 6.0,
#                                3.0e7, 10.0, 10.0, 100.0, 20.0 ],
#                              [3.0e7, -10.0, -10.0, 100.0, 9.0,
#                                3.0e7, 10.0, 10.0, 10.0, 24.0 ] ] )

# extra_strips_to_remove = [ [ [6,18,17,19,20], [7] ],
#                            [ [1,9,10], [7,8] ] ]

params_guesses = np.array( [ [ 3.0e7, -10.0, -10.0, 9.0,
                               3.0e7, 10.0, 10.0, 19.0, 100.0 ],
                             [3.0e7, -10.0, -10.0, 5.0,
                               3.0e7, 10.0, 10.0, 24.0, 100.0 ] ] )

extra_strips_to_remove = [ [ [9,16,17,18,19], [7,8,16] ],
                           [ [1,5,6,23,24,29,30], [7,8,16] ] ]

extra_data_to_remove = [ [ [30,9] ],
                         [] ]
                         


f, axarr = plt.subplots( 2, 2, figsize = ( 10, 8 ) ) 

f.suptitle( title ) 

f.subplots_adjust( wspace = 0.5 )

for i in range( 2 ) :

    
    # filter_data( hitmaps[i] )

    remove_strips( hitmaps[i], [0,31], [0,31] )

    if cut_strips : 
        remove_strips( hitmaps[i], * extra_strips_to_remove[i] )
        remove_data( hitmaps[i], extra_data_to_remove[i] )

        # print( hitmaps[i][30,9] ) 
        
    # axarr[0,i].imshow( hitmaps[i] )

    masked_array = np.ma.array( hitmaps[i], mask=np.isnan( hitmaps[i] ) )
    cmap = colorcet.m_rainbow
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

    else :
        print( 'cov is none' ) 
        

    print( '\n\n' ) 

    resid = ( hitmap - point_source_counts( result ) ) / d_hitmap 

    cmap = colorcet.m_diverging_bkr_55_10_c35
    cmap.set_bad('white',1.)
    im = axarr[ 1, i ].imshow( resid, cmap = cmap )
    divider = make_axes_locatable( axarr[1,i] )
    cax = divider.append_axes("right", size="5%", pad=0.1)
    f.colorbar(im, cax=cax)

    axarr[1,i].set_title( r'$\tilde \chi^2 = %d / %d = %.2f$'
                          % ( np.rint( chisqr ), ndata, redchisqr ) )


plt.savefig( savepath, format='png' )
    
plt.show() 
