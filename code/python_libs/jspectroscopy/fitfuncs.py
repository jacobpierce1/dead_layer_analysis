# this module contains a class for constructing fitfuncts
# with conversions between their input array params
# and dicts for interpretation / data analysis later on.


# TODO:
# * option to automatically construct a sum of fitfuncs
# * * suboption to hold some params constant when doing this


import numpy as np
import scipy.special as special
from lmfit import Model
from libjacob import meas 

# from functools import partial 


# global 'vars' that can be used for rapidly
# applying fitfuncs.

# alpha_fitfunc = fitfunc( 
# beta_fitfunc =

# gamma_fitfunc is synonymous with gaussian.
# gamma_fitfunc = 



    
    
# # class for all the fitfuncs. for maximal efficiency
# # and compabibility with scipy functions, they all
# # take array parameters as input. a function can also
# # be provided to convert the array to a dict
# # and if this is specified, another func to convert
# # such a dict to an array for future use. fitfunc
# # must as its first argument 

# class fitfunc( object ):

#     def __init__( self, fitfunc,
#                   params_array_to_dict_func = None,
#                   params_dict_to_array_func = None ):


#         if not callable( fitfunc ):
#             raise ValueError( 'fitfunc must be a function' )


#         # verify that f is a function that accepts 2 params.
#         if f.__code__.co_argcount != 2:
#             raise ValueError( 'fitfunc must take 2 parameters: ' +
#                               '( params_array, input_array )' )
        
        
#         self._fitfunc = fitfunc


#         # define the parameter conversion functions:
#         if params_array_to_dict_func is not None:

#             if not callable( params_array_to_dict_func ):
#                 raise ValueError( 'params_array_to_dict_func must be a function' )
        
#             self._array_to_dict_func = params_array_to_dict_func

            
#         if params_dict_to_array_func is not None:

#             if not callable( params_dict_to_array_func ):
#                 raise ValueError( 'params_dict_to_array_func must be a function' )

#             self._dict_to_array_func = params_dict_to_array_func


#     def dict_to_array( self, params_array ):
#         return self._dict_to_array_func( params_array )


#     def array_to_dict( self, params_array ):
#         return self._array_to_dict_func( params_array )
    

#     def apply( self, pf, input_array ):
#         return self._fitfunc( pf, input_array ) 




    
    
# # in this case you supply all args as fit params, even the det ones which should be fixed. 
# def fitfunc_alpha_free_det_params( p, x ): 

#     return alpha_fit( p[0],p[1],p[2],p[3],p[4],p[5], x ) 

# #    return np.array( map( lambda z: alpha_fit( *p, z ), x ) )
#     #return np.array( map( lambda z: alpha_fit( p[0],p[1],p[2],p[3],p[4],p[5],z ), x ) )
#     # return alpha_fit( p[0],p[1],p[2],p[3],p[4],p[5] , x )  # expand array as function args 
    



# return an inline of fitfunc_n_alpha_peaks given the number of peaks to fit.
# the first 4 params are the detector params, and the last 2*n are (A1, u1, ... )
# previously fitfunc_n_alpha_peaks_abstract
def sum_n_fitfuncs( fitfunc, n ):
    return lambda x, **params: fitfunc( n, x, **params )
    # return partial( fitfunc, n )




# # fit n alpha peaks given vector p and array x 
# # p format: sigma, eta, tau1, tau2, A1, mu1, ..., A_n, mu_n
# def fitfunc_n_alpha_peaks( n, p, x ):

#     ret = 0

#     for i in range( n ):
#         ret += alpha_fit( p[4+2*i],p[4+2*i+1],p[0],p[1],p[2],p[3], x )

#     return ret





# return an lmfit model constructed from n alpha
# fit functions. keep sigma, eta, tau1, tau2 the same
# in all of them.

def fitfunc_n_alpha_peaks( n, x, **params ):
    
    ret = 0
    for i in range( n ):
        ret += alpha_fit( params[ 'A' + str(i) ], params[ 'mu' + str(i) ],
                          params[ 'sigma' ], params[ 'eta' ], 
                          params['tau1'], params['tau2'], x )

    return ret





# construct dict for input to fitfunc_n_alpha_peaks.

def construct_n_alpha_peaks_params(
        sigma, eta, tau1, tau2, A_array, mu_array ):

    n = len( mu_array )

    param_fit_bounds = { 'sigma' : [ 0, None ],
                         'tau1' : [ 0, None ],
                         'tau2' : [ 0, None ],
                         'eta' : [ 0, 1 ] }
    
    params =  { 'sigma' : sigma, 'eta' : eta,
                'tau1' : tau1, 'tau2' : tau2 }

    for i in range( n ):

        # keys for the 2 dicts 
        mu = 'mu' + str(i)
        A = 'A' + str(i) 
        
        params[ mu ] = mu_array[i]
        params[ A ] = A_array[i]

        param_fit_bounds[ mu ] = [ 0, None ]
        param_fit_bounds[ A ] = [ 0, None ] 

    return ( params, param_fit_bounds ) 





# port the parameters of the model to a meas.meas
def get_alpha_params_dict( model ): 

    params = model.params

    npeaks = ( len( params ) - 4 ) // 2  

    ret = { 'sigma' : 0, 'tau1' : 0, 'tau2' : 0, 'eta' : 0 } 

    for key in ret:
        ret[ key ] = meas( params[key].value, params[key].stderr )

    for key_base in [ 'mu', 'A' ]:

        vals = np.empty( npeaks )
        deltas = np.empty( npeaks ) 
            
        for i in range( npeaks ):

            param = params[ key_base + str(i) ] 
            vals[i] = param.value
            deltas[i] = param.stderr

        ret[ key_base ] = meas( vals, deltas ) 
            

        #    print( ret ) 
    return ret

    
# # fit n alpha peaks given vector p and array x 
# # p format: sigma, tau, A1, mu1, ..., A_n, mu_n
# def fitfunc_n_alpha_peaks_eta1( n, p, x ):

#     ret = 0
    
#     for i in range( n ):
#         ret += alpha_fit_eta1( p[2+2*i],p[2+2*i+1],p[0],p[1], x )

#     return ret


#     # return np.sum( [ map( lambda y: alpha_fit_eta1(p[2+2*i],p[2+2*i+1],p[0],p[1], y), x )  for i in range(n) ], axis=0 )  



    
# # same as alpha_fit but with eta fixed at one. this means we are taking only one of the terms. 
# # experimented with this after noting that the fits seem indep. of eta, which can cause trouble.
# def alpha_fit_eta1( A, mu, sigma, tau, x ):
#     logtmp = (x-mu)/tau + sigma**2.0/(2*tau**2.0) + np.log( special.erfc( (x-mu)/sigma + sigma/tau) / np.sqrt(2) ) 
#     return A/(2.0*tau) * np.exp( logtmp )



# reference: equation 10 in Bortels 1987  
# this function is meant to be applied to scalar x, not list
def alpha_fit( A, mu, sigma, eta, tau1, tau2, x ):        
            
    # prevent overflow by computing logs and then exponentiating, at the expense of some
    # floating pt error. logtmpz is the log of the 2 analagous terms in the integrand.
    tauz = [tau1, tau2]
    logtmpz = [ (x-mu)/tau + sigma**2.0/(2*tau**2.0) + np.log( special.erfc( (x-mu)/sigma + sigma/tau) / np.sqrt(2) ) for tau in tauz ]
    return ( A/2.0 )  * ( (1-eta)/tau1 * np.exp(logtmpz[0])
                            + (eta/tau2) * np.exp(logtmpz[1])   ) 


