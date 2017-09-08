import numpy as np

from jacob_math import xcut

# evaluate the fwhm on the set of points provided. up to used to specify an appropriate
# number of entries to balance speed and precision. function must be scalar outupt for 
# scalar input. assumes that x is a sorted array, otherwise this function would take forever.
def fwhm_from_func( function, x, dx=[] ):
    x = np.array(x)
    print x
    print function(1) 
    y = function( x )
    return fwhm_from_data( x, y, dx )
    
    

def fwhm_from_data( x, y, dx=[] ):
    x = np.array(x) 
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
    halfmax = ymax / 2.0
    halfmax_positions = np.array( [[0.0,0.0]]*2 )
    lefty = xcut( x, y, [x[0], xmax_arr[0] ] )
    righty = xcut( x, y, [xmax_arr[-1], x[-1] ] )
    
    halfmax_positions[0][0] = x[ np.argwhere( lefty >= halfmax )[0] ]
    halfmax_positions[0][1] = x[ np.argwhere( lefty <= halfmax )[-1] ]
    halfmax_positions[1][0] = x[ np.argwhere( righty >= halfmax )[-1] ]
    halfmax_positions[1][1] = x [np.argwhere( righty <= halfmax )[0] ]

    average_halfmax_positions = [ np.average(halfmax_positions[i]) for i in range(2) ]
    average_halfmax_uncertainties = [ np.ediff1d( halfmax_positions[i] )[0] / 2.0 for i in range(2) ]
    fwhm = 0-np.ediff1d( average_halfmax_positions )[0]
    
    return ( fwhm, average_halfmax_positions, average_halfmax_uncertainties )
    
        
    
function = lambda x: 1- x**2.0
print fwhm_from_func( function, np.linspace(-1,1,100) )
