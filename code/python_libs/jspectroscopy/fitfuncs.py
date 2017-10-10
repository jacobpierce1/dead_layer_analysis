# this module contains a class for constructing fitfuncts
# with conversions between their input array params
# and dicts for interpretation / data analysis later on.


# TODO:
# * option to automatically construct a sum of fitfuncs
# * * suboption to hold some params constant when doing this


import numpy as np
import scipy.special as special



# global 'vars' that can be used for rapidly
# applying fitfuncs.

# alpha_fitfunc = fitfunc( 
# beta_fitfunc =

# gamma_fitfunc is synonymous with gaussian.
# gamma_fitfunc = 



    
    
# class for all the fitfuncs. for maximal efficiency
# and compabibility with scipy functions, they all
# take array parameters as input. a function can also
# be provided to convert the array to a dict
# and if this is specified, another func to convert
# such a dict to an array for future use. fitfunc
# must as its first argument 

class fitfunc( object ):

    def __init__( self, fitfunc,
                  params_array_to_dict_func = None,
                  params_dict_to_array_func = None ):


        if not callable( fitfunc ):
            raise ValueError( 'fitfunc must be a function' )


        # verify that f is a function that accepts 2 params.
        if f.__code__.co_argcount != 2:
            raise ValueError( 'fitfunc must take 2 parameters: ' +
                              '( params_array, input_array )' )
        
        
        self._fitfunc = fitfunc


        # define the parameter conversion functions:
        if params_array_to_dict_func is not None:

            if not callable( params_array_to_dict_func ):
                raise ValueError( 'params_array_to_dict_func must be a function' )
        
            self._array_to_dict_func = params_array_to_dict_func

            
        if params_dict_to_array_func is not None:

            if not callable( params_dict_to_array_func ):
                raise ValueError( 'params_dict_to_array_func must be a function' )

            self._dict_to_array_func = params_dict_to_array_func


    # convert params dict to an array for input to 
    def dict_to_array( self, params_array ):
        return self._dict_to_array_func( params_array )


    def array_to_dict( self, params_array ):
        return self._array_to_dict_func( params_array )
    

    def apply( self, pf, input_array ):
        return self._fitfunc( pf, input_array ) 




    
    
# in this case you supply all args as fit params, even the det ones which should be fixed. 
def fitfunc_alpha_free_det_params( p, x ): 
    return alpha_fit( p[0],p[1],p[2],p[3],p[4],p[5], x ) 

#    return np.array( map( lambda z: alpha_fit( *p, z ), x ) )
    #return np.array( map( lambda z: alpha_fit( p[0],p[1],p[2],p[3],p[4],p[5],z ), x ) )
    # return alpha_fit( p[0],p[1],p[2],p[3],p[4],p[5] , x )  # expand array as function args 
    



# return an inline of fitfunc_n_alpha_peaks given the number of peaks to fit.
# the first 4 params are the detector params, and the last 2*n are (A1, u1, ... )
# previously fitfunc_n_alpha_peaks_abstract
def n_fitfuncs_abstract( fitfunc, n ):
    return lambda p, x: fitfunc( n, p, x )




# fit n alpha peaks given vector p and array x 
# p format: sigma, eta, tau1, tau2, A1, mu1, ..., A_n, mu_n
def fitfunc_n_alpha_peaks( n, p, x ):

    ret = 0

    for i in range( n ):
        ret += alpha_fit( p[4+2*i],p[4+2*i+1],p[0],p[1],p[2],p[3], x )

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
    #if eta < 0 or eta > 1:
    #    return -1000
        
    # prevent overflow by computing logs and then exponentiating, at the expense of some
    # floating pt error. logtmpz is the log of the 2 analagous terms in the integrand.
    tauz = [tau1, tau2]
    logtmpz = [ (x-mu)/tau + sigma**2.0/(2*tau**2.0) + np.log( special.erfc( (x-mu)/sigma + sigma/tau) / np.sqrt(2) ) for tau in tauz ]
    return ( A/2.0 )  * ( (1-eta)/tau1 * np.exp(logtmpz[0])
                            + (eta/tau2) * np.exp(logtmpz[1])   ) 


