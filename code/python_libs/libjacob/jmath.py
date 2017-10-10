## includes 
# import array
import numpy as np
import scipy.optimize
import scipy.special
from math import log10, floor
import peakdetect.peakdetect as peakdetect


import sys
import warnings
from scipy.optimize import OptimizeWarning
warnings.simplefilter("error", OptimizeWarning)


# curve fitting library
from lmfit import Model, Parameters    



                                     
## gaussian 
#def fitfunc(p,x):    
#    return p[0]*np.exp(-x/p[1]) + p[2]
#

def residual_from_fitfunc( p, x, y, yerr, fitfunc ):
    return ((fitfunc(p, x)-y)/yerr)





def xcut( x, y, newx_bounds, xsorted=0 ):

    # if x is sorted, then find newx_bounds[0] from the left
    # and newx_bounds[1] from the right

    if xsorted:
        
        left = np.searchsorted( x, newx_bounds[0] )
        right = np.searchsorted( x, newx_bounds[1], side = 'right' )
        return y[ left : right ]

    # otherwise find applicable indices by brute force.
    # should be avoided.

    return np.asarray( y[ (x >= newx_bounds[0]) & (x <= newx_bounds[1]) ] )
    





# detect the positions of the 5 highest peaks. peak_positions is an
# array of tuples (peakpos, value)
# https://stackoverflow.com/questions/6910641/how-to-get-indices-of-n-maximum-values-in-a-numpy-array
# even if we found all the fits, for simplicity we still do the peak
# detection since it is the easiest way to deteremine a good location
# for the plot boundaries.

def get_n_peak_positions( n, data ):

    output_peak_positions = [0] * n
    
    # peakdetect returns 2 tuples: positions and counts of peaks
    peak_positions = peakdetect.peakdetect( data, lookahead=10 )[0]


    # print( 'peakpositions type data: ' )
    # print( type(peak_positions) )
    # print( type( peak_positions [0][1] ) ) 

    # it is possible that not all n peaks are found. in that case
    # we will still populate our_peaks with the ones that were
    # found as it may still be useful. 
    num_peaks_found = min( n, len( peak_positions ) )

        
    # now find the 5 largest and sort by x position.
    # indices is the indices of the peaks as found in data
    indices = np.argpartition( [ z[1] for z in peak_positions ],
                               -num_peaks_found )[ -num_peaks_found : ]
    
    output_peak_positions[:] = [ np.asscalar( peak_positions[z][0] )
                                 for z in sorted( indices ) ]

    
    # return number of peaks detected.
    return ( output_peak_positions, num_peaks_found )








# perform fit. return None if there is error.
np.seterr(divide='ignore', invalid='ignore')

def jacob_least_squares( x, y, dy, params_guess, fitfunc, dx = None,
                         reduc_chisq_max = np.inf, fit_bounds = None,
                         params_guess_bounds = None):
    
    
    model = Model( fitfunc, indepdendent_vars = [ 'x' ] ,
                   param_names = params_guess.keys() )

    print( (model.independent_vars, model.param_names ) )
    
    # now set up and apply the lmfit
    params = Parameters()
    for key in params_guess:
        params.add( key, value = params_guess[ key ] )

    if params_guess_bounds is not None:
        for key in params_guess_bounds:
            left, right = params_guess_bounds[key]

            if left > right :
                raise ValueError( 'The specified bounds for param %s are '
                                  + 'not ordered.' % ( key, ) )
            
            if left is not None:
                params[key].set( min=left )

            if right is not None:
                params[key].set( max=right ) 

    # model.independent_vars = [ 'x' ]
    # model.param_names = params.keys() 

    print( (model.independent_vars, model.param_names ) ) 

    
    result = model.fit( y, params, x = x, dy = dy  ) 

    

    
    # model = 

    # # construct args for the fitfunc
    # precut = [ x, y ]

    # if dx is not None:
    #     precut += dx

    #     # residual_function =
    #     raise NotImplementedError

    # else:
    #     pass
    #     # residual_function = (
    #     #     lambda params, x, y, dy:
    #     #     residual_from_fitfunc( params, x, y, dy, fitfunc )
    #     # )

    # precut.append( dy )

    # # cut the args as appropriate
    # if fit_bounds is None:
    #     fit_bounds = [ min(x), max(x) ]  # used for plots 
    #     args = precut

    # else:
    #     args = [ xcut( x, _z, fit_bounds ) for _z in precut ]

    # args = tuple( args )




    
    # # using scipy.leastsq. does not allow you to specify bounds for the fit params. 
    # try:       
    #     pf, cov, info, mesg, success =    \
    #         scipy.optimize.leastsq( residual_function, p0, args=args, full_output=1 )  
    # except ValueError:
    #     status = 0
    #     return None
    
    # dof = len( args[0] )-len(pf)
    # reduc_chisq = sum(info["fvec"]*info["fvec"]) / dof

    # # detect and handle errors 
    # error = success > 4 or len(pf)==0 or cov is None or reduc_chisq > reduc_chisq_max
    # if error:
    #     status = 0 
    #     return None

    # if (reduc_chisq_max is not None) and reduc_chisq > reduc_chisq_max:
    #     return None
    
    # pferr = np.sqrt( np.diag(cov) * reduc_chisq )
    # return ( reduc_chisq, dof, pf.tolist(), pferr.tolist() )

      
#

    # params: sigma, eta, tau1, tau2, A1, u1, ..., An, un )    
    #lower_bounds = [0.01, 0.01, 0.01, 0.01, 0.01, p0[5]-10 ]#, 0, p0[7]-10 ]
    #upper_bounds = [np.inf, 1.0, np.inf, np.inf, np.inf, p0[5]+10 ]# , np.inf, p0[7]+10 ]
    # lower_bounds = [0.01, 0.01, 0.01, p0[3]-10 ]#, 0, p0[7]-10 ]
    # upper_bounds = [np.inf, np.inf, np.inf, p0[3]+10 ]# , np.inf, p0[7]+10 ]
    
#    ###using optimize.least_squares, which allows for bounds. see http://scipy.github.io/devdocs/generated/scipy.optimize.least_squares.html
#    
#    # minimize using least_squares 
#    result = optimize.least_squares( residual_function, p0, args=(newx, newy, newdy), bounds=(lower_bounds, upper_bounds) )
#    
#    # check that fit was successful. return status detailed in above link.
#    if( result.status < 1 ):
#        print "ERROR: fit failed to converge."
#        sys.exit()
#
#    # extract results of the fit.
#    pf = result.x
#    dof = len(x) - len(pf) 
#    chisq = np.dot( result.fun, result.fun )
#    pferr = [ np.sqrt( np.diag(result.jac)) for i in range(len(pf)) ]
      
    ### optimize using curve_fit 
    #try:
    #    pf, cov = optimize.curve_fit( f=fitfunc, xdata=newx, ydata=newy, sigma=newdy, p0=p0, bounds=(lower_bounds, upper_bounds) )
    #except ValueError:
    #    print "ValueError: fit failed to converge."
    #    sys.exit()        
    #except RuntimeError:
    #    print "RuntimeError: fit failed to converge."
    #    sys.exit()        
    #except OptimizeWarning:
    #    print "OptimizeWarning: fit failed to converge."
    #    sys.exit()        
    #    
    ## obtain other info about the fit.
    #chisq = np.sum( np.square( fitfunc(pf, newx) - newy ) )
    #dof = len(x) - len(pf)  
    #pferr= np.sqrt( np.diag( cov ) )     

    # return ( chisq / dof, dof, pf, pferr )
    





    

    
# input: a function that takes array of parameters and a scalar variable x, same as input of
# optimize.least_sq; pf and pferr, obtained from jacob_least_squares; peakpos_guess, estimate of
# the peak positions; number of iterations to perform.
#
# behavior: assume that pferr are standard deviations of a normal distribution of which pf values
# are the means; do a monte carlo simulation in which an array is formed with values from those normal
# distributions. then find the maximum of the function and add it to a list.
#
# return: peakpos (average), peakpos_delta (std of mean), peakval (function at peakpos),
# peakval_delta (estimated using 2nd order taylor expansion of f; first order normally works,
# but in this case we know that f'(x) = 0 at the max so it will give 0.
def estimate_peakpos( f, p, p_delta, peakpos_guess, num_iterations=1000 ):
    peakpos_arr = np.empty( num_iterations, dtype=np.float64 )
    # print( 'p: ' + str(p) ) 
    
    for i in range(num_iterations):

        # keep picking p until the amplitude is positive
        while 1:
            current_p = np.random.normal( p, p_delta )

            # break out of the while if certain entries are
            # not physical. TODO: abstract this.
            if current_p[4] > 0 and current_p[1] < 1:
                break

        # now we construct a function from the random p
        current_inverted_f = lambda x_: 0 - f( current_p, x_ )  
        result = scipy.optimize.fmin( current_inverted_f, peakpos_guess, disp=0 )

        peakpos_arr[i] = result
        
    return peakpos_arr
    



# evaluate the fwhm on the set of points provided. up to used to specify an appropriate
# number of entries to balance speed and precision. function must be scalar outupt for 
# scalar input. assumes that x is a sorted array, otherwise this function would take forever.
def width_at_frac_of_max_from_func( function, bounds, N=2, num_samples=1000, dx=[] ):
    x = np.array(x)
    y = function( x )
    return width_at_frac_of_max( y, bounds, N, num_samples, dx )
    



# this returns the fwhm if N=2, otherwise it is the full width at 1/N fraction of the max.
def width_at_frac_of_max( y, bounds, N=2, num_samples=1000, dx=[] ):
    x = np.linspace( bounds[0], bounds[1], num_samples ) 
    y = np.array(y) 
    
    if( x.size != y.size ):
        print( "ERROR: x and y have different size." )
        print( "x size = " + str(x.size) )
        print( "y size = " + str(y.size) )
        sys.exit(0)
        
    # find positions of x 
    ymax = np.amax(y)
    xmax_arr = x[ np.argwhere( y == ymax ).flatten() ]
        
    # find rightmost and leftmost / top and bottom half max positions
    halfmax = ymax*1.0 / N
    halfmax_positions = np.array( [[0.0,0.0]]*2 )
        
    leftx = xcut( x, x, [x[0], xmax_arr[0] ] )
    rightx = xcut( x, x, [xmax_arr[-1], x[-1] ] )        
    
    lefty = xcut( x, y, [x[0], xmax_arr[0] ] )
    righty = xcut( x, y, [xmax_arr[-1], x[-1] ] )
            
    
    halfmax_positions[0][0] = leftx[ np.argwhere( lefty >= halfmax )[0] ]
    halfmax_positions[0][1] = leftx[ np.argwhere( lefty <= halfmax )[-1] ]
    halfmax_positions[1][0] = rightx[ np.argwhere( righty >= halfmax )[-1] ]
    halfmax_positions[1][1] = rightx[ np.argwhere( righty <= halfmax )[0] ]

    print( halfmax_positions )
    print( halfmax_positions[0][1] )

    average_halfmax_positions = [ np.average(halfmax_positions[i]) for i in range(2) ]
    average_halfmax_uncertainties = [ np.abs( np.ediff1d( halfmax_positions[i] )[0] ) / 2.0 for i in range(2) ]
    fwhm = np.ediff1d( average_halfmax_positions )[0]
    
    return ( fwhm, average_halfmax_positions, average_halfmax_uncertainties )
    
    
