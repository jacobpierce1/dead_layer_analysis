## includes 
# import array
import numpy as np
from scipy import optimize, special
from math import log10, floor

import sys
import warnings
from scipy.optimize import OptimizeWarning
warnings.simplefilter("error", OptimizeWarning)



    



                                     
## gaussian 
#def fitfunc(p,x):    
#    return p[0]*np.exp(-x/p[1]) + p[2]
#

def residual_from_fitfunc( p, x, y, yerr, fitfunc ):
    return ((fitfunc(p, x)-y)/yerr)


def xcut( x, y, newx_bounds ):
    x = np.array(x) 
    y = np.array(y)
    return np.array( y[ np.logical_and( (x >= newx_bounds[0]), (x <= newx_bounds[1]) ) ] )




# perform fit. return None if there is error.
def jacob_least_squares( x, y, dy, xcut_bounds, p0, fitfunc, reduc_chisq_max=None ):
    
    residual_function = lambda p, x, y, dy: residual_from_fitfunc( p, x, y, dy, fitfunc )
    
    newx = xcut( x, x, xcut_bounds) 
    newy = xcut( x, y, xcut_bounds ) 
    newdy = xcut( x, dy, xcut_bounds )
    
    # using optimize.leastsq. does not allow you to specify bounds for the fit params. 
    try:       
        pf, cov, info, mesg, success =    \
            optimize.leastsq( residual_function, p0, args=(newx, newy, newdy), full_output=1 )  
    except ValueError:
        status = 0
        return None
    
    dof = len(newx)-len(pf)
    reduc_chisq = sum(info["fvec"]*info["fvec"]) / dof
    
    if( success > 4 or len(pf)==0 or cov is None or reduc_chisq > 4):
        status = 0 
        return None

    if (reduc_chisq_max is not None) and reduc_chisq > reduc_chisq_max:
        return None
    
    pferr = np.sqrt( np.diag(cov) * reduc_chisq )
    return ( reduc_chisq, dof, pf.tolist(), pferr.tolist() )

      
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
        print "ERROR: x and y have different size."
        print "x size = " + str(x.size)
        print "y size = " + str(y.size)
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

    print halfmax_positions
    print halfmax_positions[0][1]

    average_halfmax_positions = [ np.average(halfmax_positions[i]) for i in range(2) ]
    average_halfmax_uncertainties = [ np.abs( np.ediff1d( halfmax_positions[i] )[0] ) / 2.0 for i in range(2) ]
    fwhm = np.ediff1d( average_halfmax_positions )[0]
    
    return ( fwhm, average_halfmax_positions, average_halfmax_uncertainties )
    
    
