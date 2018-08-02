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
                                show = 0, savename = None, analysis_mgr = None,
                                dets = None, fstrips = None, dpi = 50,
                                source_names = None) :

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
    if dets is None :
        dets = range( channels.num_dets )

    if fstrips is None :
        fstrips = range( channels.dimx )
        
    for d in dets :
        for x in fstrips :

            print( d,x ) 

            # flattened_data_nan_counts = np.array( [ [ np.isnan( channels[i][j][d][x].x )
            #                                           for j in range( channels.num_data_per_group[i] ) ]
            #                                         for i in range( channels.num_groups ) ] )
                                           

            # flattened_secants = np.array( [ secants[i][d][x]
            #                                 for i in range( channels.num_groups ) ] )
            sec_mask = [ np.sum( ~ np.isnan( secants[i][d][x] ) ) < 3 
                         for i in range( channels.num_groups ) ]
            # sec_mask = ~ np.isnan( secants[:,d,x] ) 
                       
            print( sec_mask ) 
            if any( sec_mask ) :
                print( 'not enough data' )
                continue
            
            # ndata = np.sum( ~( np.isnan( flattened_data )
            #                    | np.isnan( flattened_secants ) ) )

            # if ndata < 2 :
            #     print( 'not enough data' )
            #     continue
            
            params_guess = [ 2.0, -100.0] + [40 for i in range( channels.num_groups ) ] 
            args = ( channels, secants, actual_energies, d, x ) 
            ret = scipy.optimize.basinhopping( strip_objective, params_guess,
                                               minimizer_kwargs = { 'args' : args },
                                               niter = int( 1e2 ) )

            # print( ret ) 
            #print( 'ndata: ', ndata ) 
            resid = compute_resid( ret.x, channels, secants, actual_energies, d, x ) 
            # print( resid )
            ndata = len( resid )
            print( 'ndata: ', ndata ) 

            
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
                params_result_errors = None
                
            if show or savename is not None : 
                plot_strip(  channels, secants, actual_energies, params_result,
                             params_result_errors, redchisqr,
                             d, x,
                             show = show, savename = savename, analysis_mgr = analysis_mgr,
                             dpi = dpi, source_names = source_names  )  

            print( '\n\n\n' )
    analysis_mgr.save_dill( dl_estimates, 'dl_estimates_indep' )         
    
    

    
  
def plot_strip( channels, secants, actual_energies, params, params_errors, redchisqr,
                d, x, savename = None, show = 0, analysis_mgr = None,
                dpi = 50, source_names = None ) : 

    # channel_fit = edet_to_channel_fit( edet, params ) 
    # num_data_per_group = [ len( data[i] ) for i in range( len( data ) ) ]
    if params is not None : 
        a, b = params[:2]
        dl = params[2:]
        a_err, b_err = params_errors[:2]
        dl_err = params_errors[2:] 
        
    num_groups = channels.num_groups 
    num_data_per_group = channels.num_data_per_group
    max_peaks = max( num_data_per_group ) 

    no_data = 1 

    f, axarr = plt.subplots( max_peaks, num_groups,
                                     figsize = ( 10, 6 ), squeeze = 0 )

    f.subplots_adjust( wspace = 0.5, hspace = 0.5 )

    for i in range( num_groups ) :
        if source_names is not None :
            axarr[0,i].set_title( source_names[i], fontsize = 15 )

        if params is not None :
            label = r'$m_%d = %.1f \pm %.1f$ keV' % ( i, dl[i], dl_err[i] ) 
            axarr[0,i].text( 0.15, 0.9, label, transform = axarr[0,i].transAxes,
                             fontsize = 15 ) 

            
        for j in range( max_peaks ) :
            axarr[j,i].set_xlabel( r'$\sec( \theta )$', fontsize = 16 )
            axarr[j,i].set_ylabel( r'Source Distribution Average', fontsize = 16 )

            axarr[j,i].tick_params(axis='both', which='major', labelsize = 14)
                     
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

    title = 'Det %d, Strip %d: Average Source Channels' % ( d, x )
    if params is not None : 
        title +=  r', $\tilde \chi ^2 = %.2f$' % redchisqr
    f.suptitle( title, fontsize = 20, y = 1.0 )

                
    if savename is not None : 
        outdir = analysis_mgr.storage_path + '/%s/%d/' % (savename, analysis_mgr.detidx_to_detnum( d ) )
        os.makedirs( outdir, exist_ok = 1 )
        if not no_data :
            tmp_name = outdir +  '%d.eps' % x
            plt.savefig( tmp_name, dpi = dpi, format = 'eps' )

    if show :
        plt.show() 
            
    plt.close( 'all' ) 

