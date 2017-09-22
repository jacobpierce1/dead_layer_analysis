# this is a module for donig general commonly used uncertainty calculations


import numpy as np
import pandas as pd


values = [ 'value', 'delta' ]
values_index = pd.Index( values )



# class meas( object ):

#     def __init__( self, x, dx ):
#        self.x = x
#        self.dx = dx

       




# x and y must be uncorrelated for this to make sense.
def esum( x, y ):
    return measurement( x['value'] + y['value'],
                        np.sqrt( x['delta']**2 + y['delta']**2 )  )




# in: 2 options, 1: array of measurements with value and delta.
# 2: a pd.Series of value and delta arrays.
def esum_n( measurements, axis=None ):
    x, dx = measvec_to_arrays( measurements )
    return measurement( np.sum(x), np.sqrt( np.sum( dx**2 ) ) ) 




# entries of xlist must be uncorrelated 
def emean( xlist, dx=None, axis=None ):

    if dx is not None:

        # this protects against typo 
        if not np.isscalar( dx ):
            print 'ERROR: dx must be scalar'
            return None

        return measurement(  np.mean(xlist, axis=axis ), dx / np.sqrt(len(xlist)) )

    else:
        return esum_n( xlist, axis ) / xlist['value'].size
                            


    
# entries must be uncorrelated. must both be measurements.
def edivide( num, denom ):
    numx, numdx = measvec_to_arrays( num )
    denomx, denomdx = measvec_to_arrays( denom )
    value = numx / denomx
    return measurement( value,
                        np.abs(value) * np.sqrt(
                            (numdx / numx)**2.0 +
                            (denomdx / denomx )**2.0 ) )
                        
        



# convert either list of measurements or measurement of lists into two tuples.
def measvec_to_arrays( measurements ):
    
    # construct the appropriate lists.
    if isinstance( measurements, pd.Series ):
        x = measurements['value']
        dx = measurements['delta']
    else:
        x = np.array( [ measurements[i]['value'] for i in range(len(measurements)) ] )
        dx = np.array( [ measurements[i]['delta'] for i in range(len(measurements)) ] )

    return ( x, dx )



# "constructor" offers several different modes.
def measurement( x, dx ):

    # if not scalar, attempt to convert to np.array.
    if not np.isscalar(x):
        x = np.asarray(x)

    if not np.isscalar(dx):
        dx = np.asarray(dx)
        
    if type(x) is np.ndarray and np.isscalar(x):
        dx = dx * np.ones_like( x )

    return pd.Series( [x,dx], values_index )





def ndarray_or_scalar( x ):
    return np.isscalar(x) or ( type(x) is np.ndarray )
