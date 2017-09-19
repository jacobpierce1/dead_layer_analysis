# this is a module for donig general commonly used uncertainty calculations


import numpy as np
import pandas as pd


values = [ 'value', 'delta' ]
values_index = pd.Index( values )






# x and y must be uncorrelated for this to make sense.
def esum( x, y ):
    return measurement( x['value'] + y['value'],
                        np.sqrt( x['delta']**2 + y['delta']**2 )  )


# in: array of measurements with value and delta.
def esum_n( xlist ):
    return measurement( np.sum( xlist['value'] ),
                        np.sqrt( np.sum( np.asarray( xlist['delta'] )**2 ) ) )


# entries of xlist must be uncorrelated 
def emean( xlist, dx=None ):

    if dx is not None:

        # this protects against typo 
        if not np.isscalar( dx ):
            print 'ERROR: dx must be scalar'
            return None

        return measurement(  np.mean(xlist), dx / np.sqrt(len(xlist)) )

    else:
        return esum_n( xlist ) / xlist['value'].size
                            



def measurement( x, dx ):
    return pd.Series( [x,dx], values_index )
