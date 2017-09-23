# this is a class for doing all common error analysis operations. the class meas
# is constructed from a numpy ndarray x and ndarray dx, or an ndarray x and a scalar x
# in which case dx is set to an array with all values dx of same shape as x.


import numpy as np



class InputError( Exception ):
    pass


class meas( object ):

    # CONSTRUCTOR
    # https://stackoverflow.com/questions/29318459/python-function-that-handles-scalar-or-arrays
    def __init__( self, x, dx, checksize=1 ):

        x = np.asarray( x )
        dx = np.asarray( dx )

        # will be set if we detect a scalar input for x.
        scalar_input = 0
        
        # if dx.ndims is 0 then dx was input as scalar. in this case, no matter
        # what we assign dx to have same shape as x, even if x is also scalar.
        if dx.ndims == 0:

            # if both x and dx are scalars then suppress a dimension and use them.
            if x.ndims == 0:
                x = squeeze( x )
                dx = squeeze( dx )

            # if not set dx to have same shape as x
            else: 
                dx = dx * np.ones( x.shape )

        # by default do a check on the size of x and dx. i only
        # recommend skipping the checksize in the arithmetic functions
        # where efficiency matters. note that it only has to be called if
        # dx is not a scalar.
        else:
            if checksize:
                if x.shape != y.shape:
                    raise InputError( 'x and dx supplied have different shape.' )

        # set instance vars and return
        self.x = x
        self.dx = dx
        return self                    
        
        # # if not scalar, attempt to convert to np.array.
        # if not np.isscalar(x):

        #     x = np.asarray(x)

        #     # you are allowed to pass a vector x and scalar dx.
        #     if np.isscalar(dx):
        #         dx = dx * np.ones_like( x )

        #     else:
        #         dx = np.asarray( dx )
                
        # # next case: scalar supplied for x
        # else:

        #     # you are not allowed to pass a scalar x and vector dx.
        #     if not np.isscalar(dx):
        #         raise InputError( 'Cannot supply scalar x and vector dx.' )

            
        # set the instance vars and return 



    
    # overload all arithmetic operations.
    # x and y must be uncorrelated for any of these to make sense.
    def __sum__( self, y ):
        return meas( self.x + y.x,
                     np.sqrt( self.dx **2 + y.dx **2 ),
                     checksize=0 )
    

    def __sub__( self, y ):
        return meas( self.x - y.x,
                     np.sqrt( self.dx **2 + y.dx **2 ),
                     checksize=0 )


    def __mul__( self, y ):
        return meas( self.x * y.x,
                     np.sqrt( ( self.x * y.dx ) ** 2 +
                              ( y.x * self.dx ) ** 2 ),
                     checksize=0 )

    def __div__( self, y ):
        val = self.x / y.x
        return meas( val,
                     np.sqrt( ( self.dx / self.x ) **2 +
                              ( y.dx / y.x ) **2 ),
                     checksize=0 )

    def __abs__( self ):
        return meas( np.abs( self.x ),
                     self.dx,
                     checksize=0 )

    def __neg__( self ):
        return meas( 0-self.x, self.dx, checksize=0 )

    def __str__( self ):
        return 'x: %s\ndx: %s' % ( str( self.x ), str( self.dx ) )

    def __eq__( self, y ):
        return ( self.x == y.x ) & ( self.dx == y.dx )

    def __ne__( self, y ):
        return not self.__eq__( y )

    def __repr__( self ):
        return 'Measurement:\n' + str( self )

    # use the power law of calculus
    def __pow__( self, n ):
        return meas( self.x ** n,
                     n * self.x ** (n-1) * self.dx,
                     checksize=0 ) 

        
    
    # EXTENSIONS OF NUMPY / ARRAY OPERATIONS 

    # use chain rule of calculus. f and fprime must be callable.
    def apply( self, f, fprime ):
        return meas( f( self.x ),
                     fprime( self.x ) * self.dx,
                     checksize=0 )


    
    # 2: a pd.Series of value and delta arrays.
    @staticmethod
    def sum( measurements, axis=None ):
        #x, dx = measvec_to_arrays( measurements )
        return meas( np.sum( x.x, axis=axis ),
                     np.sqrt( np.sum( x.dx **2, axis=axis ) ) ) 
    
    
    
    # input: xlist, a list or numpy array of measurements.
    # entries of xlist must be uncorrelated for error to make sense.
    def mean( xlist, axis=None ):
        return meas( np.mean( x.x, axis=axis ),
                     np.sqrt( np.sum( x.dx ** 2, axis=axis ) / x.x.size )
            
            


                
    
    # # convert either list of measurements or measurement of lists into two tuples.
    # def measvec_to_arrays( measurements ):
        
        
    #     if isinstance( measurements, pd.Series ):
    #         x = measurements['value']
    #         dx = measurements['delta']
    #     else:
    #         x = np.array( [ measurements[i]['value'] for i in range(len(measurements)) ] )
    #         dx = np.array( [ measurements[i]['delta'] for i in range(len(measurements)) ] )
            
    #         return ( x, dx )
                
                

# # deprecated
# # "constructor" offers several different modes.
# def measurement( x, dx ):

#     # if not scalar, attempt to convert to np.array.
#     if not np.isscalar(x):
#         x = np.asarray(x)

#     if not np.isscalar(dx):
#         dx = np.asarray(dx)
        
#     if type(x) is np.ndarray and np.isscalar(x):
#         dx = dx * np.ones_like( x )

#     return pd.Series( [x,dx], values_index )





# def ndarray_or_scalar( x ):
#     return np.isscalar(x) or ( type(x) is np.ndarray )
