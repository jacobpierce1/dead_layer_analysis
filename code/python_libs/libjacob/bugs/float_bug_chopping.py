# this is a class for doing all common error analysis operations. the
# class meas is constructed from a numpy ndarray x and ndarray dx, or
# an ndarray x and a scalar x in which case dx is set to an array with
# all values dx of same shape as x.


import numpy as np



class InputError( Exception ):
    pass

class NotImplemented( Exception ):
    pass



class meas( object ):


    # CONSTRUCTOR
    # https://stackoverflow.com/questions/29318459/python-function-that-handles-scalar-or-arrays
    def __init__( self, x, dx, checksize=1, checktype=1, bypass_checks=0 ):

        if bypass_checks :

            self.x = x
            self.dx = dx
            return None
        

        print( 'called __init__' )

        x = np.asarray( x )
        dx = np.asarray( dx )

        # type check on input
        if checktype:

            if not np.issubdtype( x.dtype, np.number ):
                raise InputError( 'x must have numeric dtype.' )

            if not np.issubdtype( dx.dtype, np.number ):
                raise InputError( 'dx must have numeric dtype.' )

                
        # if dx.ndims is 0 then dx was input as scalar. in this case, no matter
        # what we assign dx to have same shape as x, even if x is also scalar.
        if dx.ndim == 0:

            # if both x and dx are scalars then suppress a dimension and use them.
            if x.ndim == 0:
                x = np.squeeze( x )
                dx = np.squeeze( dx )

            # if not set dx to have same shape as x
            else: 
                dx = dx * np.ones( x.shape )

                
        # by default do a check on the size of x and dx. i only
        # recommend skipping the checksize in the arithmetic functions
        # where efficiency matters. note that it only has to be called if
        # dx is not a scalar.
        else:
            if checksize:
                if x.shape != dx.shape:
                    raise InputError( 'x and dx supplied have different shape.' )
                
        # set instance vars and return
        self.x = x
        self.dx = dx
        return None                    



    def __len__( self ) :
        return len( self.x ) 

    

    
            
    # OVERLOAD ALL ARITHMETIC OPERATIONS.
    # x and y must be uncorrelated for any of these to make sense.
    def __add__( self, y ):

        print( type( y ) ) 

        print( 'called __add__' ) 

        # if hasattr( y, 'x' ):
        #     return _meas_no_checks( self.x + y.x,
        #                             np.sqrt( self.dx ** 2 + y.dx ** 2 ) )

        # else:
        #     return _meas_no_checks( self.x + y, self.dx )

        if hasattr( y, 'x' ):
        
            return meas( self.x + y.x,
                         np.sqrt( self.dx ** 2 + y.dx ** 2 ),
                         bypass_checks = 1 )
        
        else:
            return meas( self.x + y, self.dx, bypass_checks = 1 )

                                                 

    def __radd__( self, y ) :
        print( 'called __radd' )
        return self + y

    

    # other right-hand operations: commutative. 
    # __radd__ = __add__
    
    
    def __str__( self ):
        return 'x: %s\tdx: %s\n' % ( str( self.x ), str( self.dx ) )

    def __repr__( self ):
        return str( self )



    # access functions: when pulling out an index of a measurement
    # storing an ndarray, return a measurement with the corresponding
    # indexed x and dx.
    def __getitem__( self, key ):
        return meas( self.x[key], self.dx[key], bypass_checks = 1 )

    
    # value must be a meas
    def __setitem__( self, key, value ):
        self.x[key] = value.x
        self.dx[key] = value.dx
        return self




    
#########################################
#### FAST ALLOC SUBCLASS ################
#########################################

# this class is equivalent to meas except the constrcutor is more
# efficient. note that type errors resulting from using this will
# most likely result in unintelligible errors. designed for absolute
# efficiency

class _meas_no_checks( meas ):

    def __init__( self, x, dx ):

        print( 'called' ) 

        if np.isscalar(x):
            self.x = x
            self.dx = dx
            return None

        else:            
            self.x = np.asarray( x )
            self.dx = np.asarray( dx )
            return None



        





# nan = meas( np.nan, np.nan ) 

x = meas( np.arange(4), np.zeros(4) )

y = np.array( [ 0.0, 1.0 ] )

print( 'performing addition' ) 
print( y[0] + x )
