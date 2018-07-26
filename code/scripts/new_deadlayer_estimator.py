# this script computes and returns the least squares estimate of the sum of
# calibrated source and detector deadlayers and the calibration coefs for each strip
# in the detector. each strip is handled separately. 

import numpy
import scipy.optimize
import jutils.meas as meas
import numpy as np
import matplotlib.pyplot as plt 
import jspectroscopy as spec 
import os
import sys





def strip_objective( params, channels, secants, actual_energies, d, x ) :
    resid = compute_resid( params, channels, secants, actual_energies, d, x )    
    ret = np.sum( resid ** 2 )
    return ret 




def compute_resid( params, channels, secants, actual_energies, d, x ) :
    a, b =  params[0:2]
    dl = params[2:] 

    idx = 0
    resid = np.zeros( np.sum( channels.num_data_per_group ) * channels.dimy )
    
    for i in range( channels.num_groups ) :
        for j in range( channels.num_data_per_group[i] ) :
            calibrated_energies = calibrate_strip( a, b, dl[i], channels[i][j][d][x], secants[i][d][x] )
            resid[ idx : idx + channels.dimy ] = ( ( actual_energies[i][j] - calibrated_energies.x )
                                                   / calibrated_energies.dx )
            idx += channels.dimy
            
    resid = resid[ ~ np.isnan( resid ) ]
    return resid 



# calibrate one strip or pixel.
def calibrate_strip( a, b, dl, channels, secants ) :
    ret = float( a ) * channels  +b  + dl * secants 
    return ret 
    

def independent_dl_regression(  channels, secants, actual_energies,
                                show = 0, save = 0, analysis_mgr = None ) :

    # dl_estimates = spec.dssd_data_container.empty( channels.num_data_per_group, channels.num_dets,
                                                   # channels.dimx, channels.dimy ) 
    dl_estimates = [ [ [ meas.empty( channels.dimx )
                         for d in range( channels.num_dets ) ] 
                       for j in range( channels.num_data_per_group[i] ) ]
                     for i in range( channels.num_groups ) ] 

    for i in range( channels.num_groups ) :
        for j in range( channels.num_data_per_group[i] ) :
            for d in range( channels.num_dets ) :
                dl_estimates[i][j][d][:] = meas.nan
    
    # handle each det separately.
    for d in range( channels.num_dets ) :
        for x in range( channels.dimx ) :

            print( d,x ) 

            flattened_data = np.array( [ channels[i][j][d][x].x
                                         for i in range( channels.num_groups )
                                         for j in range( channels.num_data_per_group[i] ) ] )

            flattened_secants = np.array( [ secants[i][d][x]
                                            for i in range( channels.num_groups ) ] )
                                            
            
            ndata = np.sum( ~( np.isnan( flattened_data )
                               | np.isnan( flattened_secants ) ) )

            if ndata < 2 :
                print( 'not enough data' )
                continue
            
            params_guess = [ 2.0, -100.0] + [40 for i in range( channels.num_groups ) ] 
            args = ( channels, secants, actual_energies, d, x ) 
            ret = scipy.optimize.basinhopping( strip_objective, params_guess,
                                               minimizer_kwargs = { 'args' : args },
                                               niter = int( 1e2 ) )

            # print( ret ) 
            print( 'ndata: ', ndata ) 
            resid = compute_resid( ret.x, channels, secants, actual_energies, d, x ) 
            # print( resid )
            
            params_result = ret.x
            a = params_result[0]
            chisqr = ret.fun
            nfree = ndata - len( params_result )
            redchisqr = chisqr / nfree
            print( chisqr ) 
            print( redchisqr )
            
            if ret.fun / ndata < 2.1 and a > 1.8 and a < 2.3: 
                cov = ret.lowest_optimization_result.hess_inv
                                            
                print( params_result )
                            
                params_result_errors = np.sqrt( np.diag( cov ) * redchisqr )
                print( params_result_errors )

                for i in range( channels.num_groups ) :
                    for j in range( channels.num_data_per_group[i] ) :
                        dl_estimates[i][j][d][x] = meas.meas( params_result[ 2+i ],
                                                              params_result_errors[ 2+i ] )
                
            else :
                print( 'fit failed' )
                params_result = None
                
            if show or save : 
                plot_strip(  channels, secants, actual_energies, params_result, d, x,
                            show = show, save = save, analysis_mgr = analysis_mgr )  

            print( '\n\n\n' )
    analysis_mgr.save_dill( dl_estimates, 'dl_estimates' )         



    
  
def plot_strip( channels, secants, actual_energies, params,
                d, x, savepath = None, show = 0, save = 0, analysis_mgr = None ) : 

    # channel_fit = edet_to_channel_fit( edet, params ) 
    # num_data_per_group = [ len( data[i] ) for i in range( len( data ) ) ]
    if params is not None : 
        a, b = params[:2]
        dl = params[2:]
    
    num_groups = channels.num_groups 
    num_data_per_group = channels.num_data_per_group
    max_peaks = max( num_data_per_group ) 

    no_data = 1 

    f, axarr = plt.subplots( max_peaks, num_groups,
                                     figsize = ( 12,8 ), squeeze = 0 )

    f.subplots_adjust( wspace = 0.5, hspace = 0.5 )

    for i in range( num_groups ) :
        for j in range( max_peaks ) :

            if j < num_data_per_group[i] :

                tmp1 = channels[i][j][d][x].x
                tmp2 = channels[i][j][d][x].dx
                sec = secants[i][d][x]

                axarr[j,i].errorbar( sec, tmp1, tmp2, ls = 'none' )
                
                if params is not None :
                    mask = ~ ( np.isnan( tmp1 ) | np.isnan( sec ) )

                    if np.sum( mask ) == 0 :
                        return
                    
                    min_sectheta = min( sec[ mask ] ) 
                    max_sectheta = max( sec[ mask ] ) 
                    min_chan = ( ( actual_energies[i][j] - dl[i] * min_sectheta ) - b ) / a
                    max_chan = ( ( actual_energies[i][j] - dl[i] * max_sectheta ) - b ) / a
                    axarr[j,i].plot( [ min_sectheta, max_sectheta ], [ min_chan, max_chan ] ) 

                nan_count = np.count_nonzero(
                    np.isnan( sec )
                    | np.isnan( tmp1 ) ) 

                if nan_count < len( tmp1 ) :
                    no_data = 0

            else :
                axarr[j,i].axis( 'off' )

    if save : 
        outdir = analysis_mgr.storage_path + '/dl_regression/%d/' % analysis_mgr.detidx_to_detnum( d )
        os.makedirs( outdir, exist_ok = 1 )
        if not no_data : 
            plt.savefig( outdir +  '%d.png' % x )

    if show :
        plt.show() 
            
    plt.close( 'all' ) 

