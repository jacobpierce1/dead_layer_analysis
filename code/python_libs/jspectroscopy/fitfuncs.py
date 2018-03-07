# this module contains a class for constructing fitfuncts
# with conversions between their input array params
# and dicts for interpretation / data analysis later on.


# TODO:
# * option to automatically construct a sum of fitfuncs
# * * suboption to hold some params constant when doing this

import libjacob.jmath as jmath
from libjacob.jmath import xcut

import numpy as np
import scipy.special as special
# from lmfit import Model

import libjacob.jmeas as meas

import matplotlib.pyplot as plt

import scipy.optimize
from scipy.stats import chi2



from enum import Enum

class peak_types( Enum ) :
    a = 0
    b = 1 
    g = 2



_fitters_num_params = { 'a' : (2,4),
                         'b' : None,
                         'g' : None,
                         'cb' : None } 

     
# __alpha_num_peak_params = 2
# __alpha_num_det_params = 4


def _resid( f, params, x, y, dy ) :
    ret = ( f( params, x ) - y ) / dy

    # print( dy ) 
    # print( ret )
    
    return ret 






def alpha_fitter( peak_params, det_params, x ) :
    return alpha_fit( *peak_params, *det_params, x ) 


def beta_fitter( peak_params, det_params, x ) :
    raise NotImplementedError()


def gamma_fitter( peak_params, det_params, x ) :
    raise NotImplementedError()


_available_fits = [ 'a', 'b', 'g' ]

_num_available_fits = len( _available_fits ) 


_all_fitters = { 'a' : alpha_fitter,
                'b' : beta_fitter,
                'g' : gamma_fitter } 



class spectrum_fitter( object ) :


    def __init__( self, peak_types, constrain_det_params = None ) :

        if constrain_det_params is None :
            constrain_det_params = { 'a' : 0, 'b' : 0, 'g' : 0 } 

        self.constrain_det_params = constrain_det_params
            
        self.peak_types = peak_types

        self.num_peaks = len( peak_types )
        
        # self.peak_type_indices = {}

        # for t in peak_types  :
        #     peak_type_indices[t] = [ i for i, x in enumerate( peak_types )
        #                              if x == t ]
            
        
        self.fitters = [ _all_fitters[t] for t in peak_types ] 

        self.params_array_peak_indices = np.empty( ( self.num_peaks, 2 ), dtype='int' ) 
        self.params_array_det_indices = np.empty( ( self.num_peaks, 2 ), dtype='int' )

        # set the indices of the relevant params for the flattened array

        peak_idx = 0 
        det_idx = 0

        constrained_det_param_indices = { 'a' : -1, 'b' : -1, 'g' : -1 }
        
        for i in np.arange( self.num_peaks ) :

            peak_type = self.peak_types[i]

            num_peak_params, num_det_params = _fitters_num_params[ peak_type ]
            
            self.params_array_peak_indices[i] = np.array( [ peak_idx,
                                                            peak_idx + num_peak_params ] )
            peak_idx += num_peak_params

            if constrain_det_params[ peak_types[i] ] :

                constrained_det_param_idx = constrained_det_param_indices[ peak_type ]
                
                if constrained_det_param_idx == -1 :

                    self.params_array_det_indices[i] = np.array( [ det_idx,
                                                                   det_idx + num_det_params ] )

                    constrained_det_param_indices[ peak_type ] = det_idx
                    det_idx += num_det_params

                else : 
                    self.params_array_det_indices[i] = np.array( [ constrained_det_param_idx,
                                                                   constrained_det_param_idx
                                                                   + num_det_params ] )

            else : 
                self.params_array_det_indices[i] = ( det_idx, det_idx + num_det_params )
                det_idx += num_det_params
            
    
        self.params_array_det_indices += peak_idx
        
        self.params_array_size = peak_idx + det_idx
                
        self.peak_params = [0] * self.num_peaks
        self.det_params = [0] * self.num_peaks
        
        


        

    def fit( self, x, y, dy, 
             peak_params_guess, det_params_guess, xbounds = None,
             fit_acceptor = None, params_shuffler = None,
             ax = None, plot_bounds = None, logscale = 1 ) :

        if ax is not None:
            ax.plot( x, y, ls = 'steps', zorder = 1, c = 'k', linewidth = 0.1 ) 

            if( logscale ):
                ax.set_yscale('log')
                ax.set_ylim( bottom = 1 ) 
            
            if plot_bounds is not None :
                ax.set_xlim( * plot_bounds ) 
            
        self.fit_acceptor = fit_acceptor
        self.params_shuffler = params_shuffler 
        self.ax = ax
        
        self.fit_attempted = 1
        
        if xbounds is not None :
            y = xcut( x, y, xbounds )
            dy = xcut( x, dy, xbounds )
            x = xcut( x, x, xbounds ) 

        self.x = x
        
        for peak_type, det_params in det_params_guess.items() : 

            # convert everything to an array
            det_params_guess_arr = np.array( det_params )
            det_params_guess[ peak_type ] = det_params_guess_arr
            
            if self.constrain_det_params[ peak_type ] : 
                det_params_guess[ peak_type ] = det_params_guess_arr
            
        for i in np.arange( self.num_peaks ) :

            peak_type = self.peak_types[i]
            
            self.peak_params[i] = np.asarray( peak_params_guess[i] )

            if self.constrain_det_params[ peak_type ] : 
                self.det_params[i] = det_params_guess[ peak_type ]  
            else :
                self.det_params[i] = np.copy( det_params_guess[ peak_type ] ) 
                
        params_array = self.__construct_params_array()

        print( params_array ) 

        # print( 'det_params: ' + str( self.det_params ) )
        # print( 'peak_params: ' + str( self.peak_params ) )
        # print( 'params_array: ' + str( params_array ) ) 

        # # call scipy.curve_fit with appropriate input
        # params_final, cov = scipy.optimize.curve_fit( self.__fit_eval, x, y,
        #                                               p0 = params_array,
        #                                               sigma = dy )

        objective = lambda _params, _x, _y, _dy : _resid( self.__fit_eval, _params,
                                                          _x, _y, _dy ) 

        # print( objective( x, y, dy ) )

        ret = scipy.optimize.leastsq( objective, params_array, args = (x,y,dy),
                                      full_output = 1 )

        params_result, cov, info, msg, status = ret

        if fit_acceptor is not None :
            fit_accept_status = fit_acceptor( x, y, dy, self ) 

        else :
            fit_accept_status = 1 
            
        success = ( status >= 1 and status <= 4
                    and ( cov is not None )
                    and fit_accept_status )

        self.success = success

        if success :

            print( 'successful fit' ) 

            print( params_result )
            # print(cov) 

            self.chisqr = np.sum( info['fvec']**2 )
            self.nfree = len( x ) - len( params_array ) 
            self.redchisqr = self.chisqr / self.nfree

            self.params_result_error = np.sqrt( np.diag( cov ) * self.redchisqr )

            print( self.params_result_error ) 

            print( 'redchisqr: ' + str(self.redchisqr) )

            self.pvalue = 1 - chi2.cdf( self.chisqr, self.nfree )
            print( 'pvalue: ' + str( self.pvalue ) )
            
            # print( self.peak_params )
            # print( self.det_params ) 
            
            if ax is not None:
                ax.plot( x, self.eval( x ), c = 'r', zorder = 2 ) 

        

                
    def eval( self, x ) :

        ret = 0
        for i in np.arange( self.num_peaks ) :
            ret += self.fitters[i]( self.peak_params[i], self.det_params[i], x ) 

        return ret

    

    
    def __fit_eval( self, params_array, x ) :

        self.__set_params_from_params_array( params_array )

        return self.eval(x)

    
            
    def __construct_params_array( self ) :

        p = np.empty( self.params_array_size )

        for i in np.arange( self.num_peaks ) :

            # print( self.params_array_peak_indices[i] )
            # print( self.peak_params[i] ) 

            p[ slice( * self.params_array_peak_indices[i] ) ] = self.peak_params[i]
            p[ slice( * ( self.params_array_det_indices[i] ) ) ] = self.det_params[i]

        return p 


    

    def __set_params_from_params_array( self, params_array ) :

        for i in np.arange( self.num_peaks ) :

            self.peak_params[i][:] = params_array[ slice( * self.params_array_peak_indices[i] ) ]
            self.det_params[i][:] = params_array[ slice( * self.params_array_det_indices[i] ) ]
            


    def plot( self, ax, c = 'r', **kw_args ) :
        ax.plot( self.x, self.eval( self.x ), c=c, **kw_args )

        



        

# reference: equation 10 in Bortels 1987  
# this function is meant to be applied to scalar x, not list

def alpha_fit( A, mu, sigma, eta, tau1, tau2, x ):
            
    # prevent overflow by computing logs and then exponentiating, at the expense of some
    # floating pt error. logtmpz is the log of the 2 analagous terms in the integrand.
    tauz = [tau1, tau2]

    logtmpz = [ (x-mu)/tau + sigma**2.0/(2*tau**2.0)
                + np.log( special.erfc( (x-mu)/sigma + sigma/tau) / np.sqrt(2) )
                for tau in tauz ]

    return ( A / 2.0 ) * ( (1-eta)/tau1 * np.exp(logtmpz[0])
                           + (eta/tau2) * np.exp(logtmpz[1])   ) 











def fit_spectrum( peak_types, x, y, dy, xbounds,
                  peak_params_guess, det_params_guess,
                  constrain_det_params = None,
                  fit_acceptor = None,
                  params_shuffler = None,
                  ax = None, plot_bounds = None, logscale = 1 ) :

    spec_fitter = spectrum_fitter( peak_types, constrain_det_params )

    spec_fitter.fit( x, y, dy, peak_params_guess,
                     det_params_guess, xbounds,
                     fit_acceptor, params_shuffler,
                     ax, plot_bounds, logscale = logscale ) 
    
    return spec_fitter










def auto_fit_spectrum( x, y, dy,
                       group_ranges, peak_locations,
                       num_peaks_to_detect, primary_peak_detector,
                       peak_sizes_guesses, det_params_guesses, peak_mu_offset,
                       peak_position_tolerance = None, db_path = None,
                       fit_acceptor = None,
                       params_shuffler = None,
                       ax = None,
                       rel_plot_bounds = None,
                       logscale = 1 ) :

    num_groups = len( group_ranges )

    if len(peak_locations) != num_groups : 

        print( '''error: inconsistent size of group_ranges, peak_structures, or peak_locations.
        they should all have length equal to group_ranges''' )
        sys.exit(0)

    peaks_per_group = [ len(peak_locations[i]) for i in range( num_groups ) ]
    
    # find main peaks
    our_peaks = jmath.get_n_peak_positions( num_peaks_to_detect, y )

    print( our_peaks ) 
    
    primary_peaks = primary_peak_detector( our_peaks, y )

    if primary_peaks is None :
        return None 
    
    print( primary_peaks )
    
    # print( primary_peaks ) 
    
    num_peaks_found = len( our_peaks )


    # do a check on peak values: energy differenc for the pairs should
    # be constant. no check necessary on the largest peak since it
    # always dominates, the only potential problem is really the
    # smallest peak in the lowest energy pair. this occurs after
    # plotting so that we can return with the plot alread made.
    # basically, we are saying that if we cannot determine where the
    # peak positions are to 0th order, we are not going to bother
    # fitting them since it will definitely not work if we misidentify
    # a peak position.
    
        
        
    # determine which fits to perform, 1 = fit must be attempted. by default if no 
    # conn is supplied we process all the fits.

    fit_attempts = [ -1 ] * num_groups  # first fit to try, which is the one we left off on.

    fits_to_perform = [ 1 ] * num_groups 

    reduc_chisq_all = [] # store all reduced chisq values.
    

    fit_bounds = [ np.array( group_ranges[a] ) + primary_peaks[a]
                   for a in range( num_groups )  ]

    if rel_plot_bounds is not None :
        plot_bounds = [ fit_bounds[0][0] + rel_plot_bounds[0],
                        fit_bounds[-1][1] + rel_plot_bounds[1] ] 

    else :
        plot_bounds = None
        
    # list of detected peaks, to be added to DB (not yet implemented) 
    peak_detect = [ our_peaks[0:2], our_peaks[2:4], our_peaks[4:] ]

    ret = [0] * num_groups 
    
    # loop through the fits that were not in the db and add them if successful.
    for i in range( num_groups ):

        if fits_to_perform[i]:

            # print( fit_bounds[i] )
            print( primary_peaks[i] ) 

            mu_array_guess = ( primary_peaks[i]
                               + np.array( peak_locations[i] )
                               + peak_mu_offset )
            
            peak_params_guess = [ [ peak_sizes_guesses[i][d], mu_array_guess[d] ]
                                  for d in range( peaks_per_group[i] ) ]


            print( peak_params_guess )

                        
            spec_result = fit_spectrum( [ 'a' ] * peaks_per_group[i],
                                        x, y, dy, fit_bounds[i],
                                        peak_params_guess, det_params_guesses[i],
                                        constrain_det_params = { 'a' : 1 },
                                        params_shuffler = params_shuffler,
                                        fit_acceptor = fit_acceptor,
                                        ax = ax,
                                        plot_bounds = plot_bounds,
                                        logscale = logscale )

            ret[i] = spec_result 
        
    return 1






        
        



def auto_fit_many_spectra() :

    pass

    # # option to pass None as the DB 
    # if db is not None:

    #     for i in range( num_groups ):
            
    #         # extract
    #         db_data = db.read_fit_data( x, y, i )

    #         successful_fit = db_data[ 'successful_fit' ]

    #         # signal that this fit will have to be re-attempted.
    #         fits_to_perform[i] = not successful_fit

    #         # last attempt that we left off on 
    #         fit_attempts[i] = db_data[ 'last_attempt' ]
            
    #         if successful_fit:
    #             # fitfunc = spec.sum_n_fitfuncs( spec.fitfunc_n_alpha_peaks, NUM_PEAKS_PER_FEATURE[i] )
                
    #             model = db_data[ 'model' ] 
    #             jplt.add_fit_to_plot( ax, x, db_data[ 'fit_bounds' ],
    #                                   jmath.model_func( model ),
    #                                   logscale = 1 )

    #             # print( model.params ) 
                
    #             reduc_chisq_all.append( model.redchi )
                

    #     # no further processing required if all the fits converged.
    #     if not any( fits_to_perform ):
    #         _add_text( ax, x, y, reduc_chisq_all )
    #         return 1
                
         

    pass 
