# this is a class for doing all common error analysis operations. the class meas
# is constructed from a numpy ndarray x and ndarray dx, or an ndarray x and a scalar x
# in which case dx is set to an array with all values dx of same shape as x.


import numpy as np



class InputError( Exception ):
    pass

class NotImplemented( Exception ):
    pass



class meas( object ):

    # CONSTRUCTOR
    # https://stackoverflow.com/questions/29318459/python-function-that-handles-scalar-or-arrays
    def __init__( self, x, dx, checksize=1, checktype=1 ):

        x = np.asarray( x )
        dx = np.asarray( dx )

        # type check on input
        if checktype:

            if not np.issubdtype( x.dtype, np.number ):
                raise InputError( 'x must have numeric dtype.' )

            if not np.issubdtype( dx.dtype, np.number ):
                raise InputError( 'dxx must have numeric dtype.' )

                
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


    # construct a new measurement from a list of other measurements
    @classmethod
    def from_list( cls, measlist ):
        return _meas_no_checks( [ measlist[i].x for i in range(len(measlist)) ],
                                [ measlist[i].dx for i in range(len(measlist)) ] )
    
        
    
    # OVERLOAD ALL ARITHMETIC OPERATIONS.
    # x and y must be uncorrelated for any of these to make sense.
    def __add__( self, y ):

        if np.isscalar( y ):
            return _meas_no_checks( self.x + y,
                                    self.dx )
                                     
        else:
            return _meas_no_checks( self.x + y.x,
                                    np.sqrt( self.dx **2 + y.dx **2 ) )
        
        
    def __sub__( self, y ):
        
        if np.isscalar( y ):
            return _meas_no_checks( self.x - y,
                                    self.dx )
        
        else:
            return _meas_no_checks( self.x - y.x,
                                    np.sqrt( self.dx ** 2 + y.dx ** 2 ) )
        
        
    def __mul__( self, y ):
        
        if np.isscalar( y ):
            return _meas_no_checks( self.x * y,
                                    self.dx * y )
        
        else:
            return _meas_no_checks( self.x * y.x,
                                    np.sqrt( ( self.x * y.dx ) ** 2 +
                                             ( y.x * self.dx ) ** 2 ) )
        
        
    def __div__( self, y ):

        if np.isscalar(y):
            return _meas_no_checks( self.x / y,
                         self.dx / y )

        else:
            val = self.x / y.x
            return _meas_no_checks( val,
                                    np.abs( val ) * np.sqrt( ( self.dx / self.x ) ** 2 +
                                                             ( y.dx / y.x ) ** 2 ) )
        

    # only arithmetic operation that needs to be defined differently if
    # the measurement is on the right side.
    def __rdiv__( self, y ):

        if np.isscalar( y ):
            return y * ( self ** -1 )

        else:
            return __div__( y, self )
        
    # other right-hand operations: commutative. 
    __radd__ = __add__
    __rsub__ = __sub__
    __rmul__ = __mul__
    
    
    # other common operations 
    def __abs__( self ):
        return _meas_no_checks( np.abs( self.x ),
                                self.dx )

    def __neg__( self ):
        return _meas_no_checks( 0-self.x, self.dx )

    def __str__( self ):
        return 'x: %s\ndx: %s' % ( str( self.x ), str( self.dx ) )

    def __eq__( self, y ):
        return ( self.x == y.x ) and ( self.dx == y.dx )

    def __ne__( self, y ):
        return not self.__eq__( y )

    def __repr__( self ):
        return 'Measurement:\n' + str( self )

    # use the power law of calculus
    def __pow__( self, n ):
        return _meas_no_checks( self.x ** n,
                                np.abs( n * self.x ** (n-1) ) * self.dx )

        
    
    # EXTENSIONS OF NUMPY / ARRAY OPERATIONS 

    # use chain rule of calculus. f and fprime must be callable.
    def apply( self, f, fprime ):
        return _meas_no_checks( f( self.x ),
                                np.abs( fprime( self.x ) ) * self.dx )


        
    # sum the entries of the measurement along specified axis.
    # input is one measurement with an x value that is a vector,
    # not a list of measurements. for that see the non-class method
    # sum.
    def sum( measurements, axis=None ):
        #x, dx = measvec_to_arrays( measurements )
        return _meas_no_checks( np.sum( x.x, axis=axis ),
                                np.sqrt( np.sum( x.dx **2, axis=axis ) ) )
    
    
    # analagous to above sum(). a non-class mean is implemented which can
    # be used for a list input.
    def mean( xlist, axis=None ):
        return _meas_no_checks( np.mean( x.x, axis=axis ),
                                np.sqrt( np.sum( x.dx ** 2, axis=axis ) / x.x.size ) )



    # access functions: when pulling out an index of a measurement
    # storing an ndarray, return a measurement with the corresponding
    # indexed x and dx.
    def __getitem__( self, key ):
        return _meas_no_checks( self.x[key], self.dx[key] )

    # value must be a meas
    def __setitem__( self, key, value ):
        self.x[key] = value.x
        self.dx[key] = value.dx
        return self

    def __delitem__( self, key ):
        raise NotImplemented( 'Have not decided on best functionality here.' )


    

#########################################
#### FAST ALLOC SUBCLASS ################
#########################################
# this class is equivalent to meas except the constrcutor is more
# efficient. note that type errors resulting from using this will
# most likely result in unintelligible errors. designed for absolute
# efficiency
class _meas_no_checks( meas ):

    def __init__( self, x, dx ):

        if np.isscalar(x):
            self.x = x
            self.dx = dx
            return None

        else: 
            self.x = np.asarray( x )
            self.dx = np.asarray( dx )
            return None




############################################
# NON CLASS METHODS ########################
############################################
# common functions
def cos( _meas ):
    return _meas.apply( np.cos, np.sin )

def sin( _meas ):
    return _meas.apply( np.sin, np.cos )

def tan( _meas ):
    return _meas.apply( np.tan, lambda x: 1 / ( np.cos(x) ** 2 ) )

def arccos( _meas ):
    return _meas.apply( np.arccos,
                        lambda x: 1 / np.sqrt( 1 - x ** 2 ) )

def arcsin( _meas ):
    return _meas.apply( np.arcsin,
                        lambda x: 1 / np.sqrt( 1 - x ** 2 ) )

def arctan( _meas ):
    return _meas.apply( np.arctan,
                        lambda x: 1 / ( 1 + x ** 2 ) )

def log( _meas ):
    return _meas.apply( np.log,
                        lambda x: 1 / x )

# sum and mean are overloaded with the class instance methods. the difference
# is those methods act on the entries of a single measurement object, whereas these
# act on an input list along the specified axis. both of these return a single
# measurement with the same dimensions as the measurements in the input list.
# note that the 
def sum( measlist, axis=0 ):
    return _meas_no_checks( np.sum( [ measlist[i].x for i in np.arange(len(measlist)) ],
                                    axis=axis ),
                            np.sqrt( np.sum( [ measlist[i].dx
                                               for i in np.arange( len( measlist ) ) ],
                                             axis = axis ) ) )

def mean( measlist, axis=0 ):
    return _meas_no_checks( np.mean( [ measlist[i].x for i in np.arange(len(measlist)) ],
                                     axis=axis ),
                            np.sqrt( np.sum( [ measlist[i].dx
                                               for i in np.arange( len( measlist ) ) ],
                                             axis = axis ) ) )

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
