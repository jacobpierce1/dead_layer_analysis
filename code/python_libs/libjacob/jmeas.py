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


    
    @classmethod
    def empty( cls, tup ):
        return cls( np.empty( tup ), np.empty( tup ) )
    

    
    # construct a new measurement from a list of other measurements.
    @classmethod
    def from_list( cls, measlist ):
        return _meas_no_checks( [ measlist[i].x for i in range(len(measlist)) ],
                                [ measlist[i].dx for i in range(len(measlist)) ] )

    @classmethod
    def from_array( self, meas_array ):

        ret = np.empty( meas_array.shape )

        return self.__from_array_recurs( meas_array, meas.empty( meas_array.shape ) )
        

    @classmethod
    def __from_array_recurs( self, meas_array, ret ):

        if hasattr( meas_array[0], 'dx' ):
            return self.from_list( meas_array )
        
        for i in range( meas_array.shape[0] ):
            ret[i] = self.__from_array_recurs( meas_array[i], ret[i] )

        return ret
    
    # # construct a new measurement from an ndarr along specified axis. if not specified,
    # # then construct by assuming that measurements are stored in the most deeply
    # # nested entries.
    # def from_ndarray( cls, array, axis = -1 ):

    #     if axis == -1:
    #         axis = array.ndim

    #     return _meas_no_checks( 

    
            
    # OVERLOAD ALL ARITHMETIC OPERATIONS.
    # x and y must be uncorrelated for any of these to make sense.
    def __add__( self, y ):

        if hasattr( y, 'x' ):
            return _meas_no_checks( self.x + y.x,
                                    np.sqrt( self.dx ** 2 + y.dx ** 2 ) )

        else:
            return _meas_no_checks( self.x + y,
                                    self.dx )
                                                 
        
        
    def __sub__( self, y ):
        
        if hasattr( y, 'x' ):
            return _meas_no_checks( self.x - y.x,
                                    np.sqrt( self.dx ** 2 + y.dx ** 2 ) )

        else:
            return _meas_no_checks( self.x - y,
                                    self.dx )
        
        
        
    def __mul__( self, y ):
        
        if hasattr( y, 'x' ):
            return _meas_no_checks( self.x * y.x,
                                    np.sqrt( ( self.x * y.dx ) ** 2 +
                                             ( y.x * self.dx ) ** 2 ) )

        else:
            return _meas_no_checks( self.x * y,
                                    self.dx * y )
                
        
        
    def __truediv__( self, y ):

        if hasattr( y, 'x' ):
            val = self.x / y.x
            return _meas_no_checks( val,
                                    np.abs( val ) * np.sqrt( ( self.dx / self.x ) ** 2 +
                                                             ( y.dx / y.x ) ** 2 ) )
        else:
            return _meas_no_checks( self.x / y,
                                    self.dx / y )
        
    # __floordiv__ is not defined since it doesn't make sense for this
    # class.
        

    # only arithmetic operation that needs to be defined differently if
    # the measurement is on the right side.
    def __rtruediv__( self, y ):

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


    
    # very useful function here. handles uncertainty calculation for a map
    # f: R^n to R with continuous first derivatives. the coordinates must
    # be indepdendent for this to make sense, but a typical function can always
    # be written in such a way that this is the case.
    #
    # input:
    #   1: function from R^n to R. No reference to x or dx should occur. 
    #   2: tuple (preferred) or list containing the first partial derivatives of f in order.
    #      So the first coordinate needs to be df/dx_0, and so on.
    # output:
    #   uncertainty using the chain rule of several variables.
    #
    # info: there is no check on the size of fprime_tuple for maximal efficiency. you
    # will probably get an index out of range error if the tuple is too small
    # and no error otherwise.
    
    def apply_nd( self, f, fprime_tuple ):

        val = f( self.x )
        partials_evaluated = fprime_tuple( self.x )

        delta = np.sqrt( np.sum( [ ( partials_evaluated[i] * self.dx[i] ) ** 2
                                   for i in np.arange( len( partials_evaluated ) ) ],
                                 axis = 0 ) )

        return _meas_no_checks( val, delta )
                         


    
    # sum the entries of the measurement along specified axis.
    # input is one measurement with an x value that is a vector,
    # not a list of measurements. for that see the non-class method
    # sum.
    def sum( self, axis=None ):
        #x, dx = measvec_to_arrays( measurements )
        return _meas_no_checks( np.sum( self.x, axis=axis ),
                                np.sqrt( np.sum( self.dx ** 2, axis=axis ) ) )
    
    
    # analagous to above sum(). a non-class mean is implemented which can
    # be used for a list input.
    def mean( self, axis=None ):

        if np.isscalar( self.x ):
            return self
        
        if axis is None:
            num_entries = self.x.size
        else:
            num_entries = self.x.shape[ axis ]

        return _meas_no_checks( np.mean( self.x, axis=axis ),
                                np.sqrt( np.sum( self.dx ** 2, axis=axis ) ) / num_entries )


    
    # same as mean but call np.nanmean
    def nanmean( self, axis = None, option = None ):

        if axis is not None:
            raise NotImplemented
        
        if np.isscalar( self.x ):
            return self
        
        if axis is None:
            num_entries = self.x.size
        else:
            num_entries = self.x.shape[ axis ]


        if option is None:
            mean_x = np.nanmean( self.x, axis=axis )
            mean_dx = ( np.sqrt( np.nansum( self.dx ** 2, axis=axis ) )
                        / len( self.x != np.nan ) )

        elif option == 'weighted' :
            one_over_dx_sq = self.dx ** 2  # save some arithmetic
            mean_x = np.nanmean( self.x, axis = axis, weights = one_over_dx_sq )
            mean_dx = ( np.sqrt( np.nanmean( self.x ** 2,
                                             axis = axis,
                                             weights = one_over_dx_sq ** 2 ) ) 
                        / np.nansum( one_over_dx_sq, axis = axis ) )
            
        return _meas_no_checks( mean_x, mean_dx ) 
                                

                        
    
    # take std of the elements about specified axis.
    # assume that the self.x is a np.ndarray.
    # todo: someone verify that this is correct.
    def std( self, axis=None ):

        if axis is not None:
            raise NotImplemented

        std_x = np.std( self.x, axis = None )

        if axis is not None:
            N = self.x.shape[axis]
        else:
            N = self.x.size

        # formula derived from error analysis.
        std_dx = 1 / ( N * std_x ) * np.sqrt(
            np.sum( ( self.x - np.mean( self.x ) * self.dx ) ** 2,
                    axis = axis ) )
        
        return _meas_no_checks( std_x, std_dx )
                                

    # same as std but call nanstd
    def nanstd( self, axis=None ):

        if axis is not None:
            raise NotImplemented

        std_x = np.nanstd( self.x, axis = None )

        if axis is not None:
            N = self.x.shape[axis]
        else:
            N = ( self.x != np.nan ).size

        # formula derived from error analysis.
        std_dx = 1 / ( N * std_x ) * np.sqrt(
            np.nansum( ( ( self.x - np.nanmean( self.x ) ) * self.dx ) ** 2,
                       axis = axis ) )
        
        return _meas_no_checks( std_x, std_dx )

    

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

    def shape( self ):
        return self.x.shape

    def size( self ):
        return self.x.size

    

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


# sum and mean are overloaded with the class instance methods. the
# difference is those methods act on the entries of a single
# measurement object, whereas these act on an input list along the
# specified axis. both of these return a single measurement with the
# same dimensions as the measurements in the input list.  note that
# the
def sum( measlist, axis=0 ):
    return _meas_no_checks( np.sum( [ measlist[i].x for i in np.arange(len(measlist)) ],
                                    axis=axis ),
                            np.sqrt( np.sum( [ measlist[i].dx ** 2 
                                               for i in np.arange( len( measlist ) ) ],
                                             axis = axis ) ) )

# take a list of measurements and report the mean about
# specified axis.
def mean( measlist, axis=0 ):
    return _meas_no_checks( np.mean( [ measlist[i].x for i in np.arange(len(measlist)) ],
                                     axis=axis ),
                            np.sqrt( np.sum( [ measlist[i].dx
                                               for i in np.arange( len( measlist ) ) ],
                                             axis = axis ) ) )

# only defined for 1D and 2D matrices.
def dot( x, y ):

    ndim = x.x.ndim

    # general case of 2x2 matrices.
    if ndim == 2:
        xshape = x.x.shape
        yshape = y.x.shape
        
        # value is easy to get.
        val = np.dot( x.x, y.x )
        
        # compute the uncertainty
        delta = np.empty( (xshape[0], yshape[1] ) )

        for i in range( xshape[0] ):
            for j in range( yshape[1] ):
                delta[i,j] = np.sqrt( np.sum( ( x.x[ i,: ] * y.dx[ :,j ] ) ** 2 )  + 
                                      np.sum( ( x.dx[ i,: ] * y.x[ :,j] ) ** 2 )  )

        return meas( val, delta )


    # same thing but easier here.
    elif ndim == 1:

        val = np.dot( x.x, y.x )
        size = len( x.x )
        delta = np.emppty( size ) 
        
        for i in range( size ):
            delta[i] = np.sqrt( np.sum( ( x.x * y.dx ) ** 2 )  + 
                                np.sum( ( x.dx * y.x ) ** 2 )  )
            

    # scalar case
    elif ndim == 0:
        return x * y

    else:
        raise NotImplemented( 'meas.dot not implemented for 2 or greater dimensions.' )

    
        

    
# concatenate the current measurement to another
# measurement. since all measurements are constructed with the
# right shape, we don't have to do a check on the shapes of these
# things.  todo: determine the most efficient way to implement this.

def append( x, y ):
    
    retx = np.append( x.x, y.x )
    retdx = np.append( x.dx, y.dx )

    return _meas_no_checks( retx, retdx )



# overload abs. x.dx is already positive.
def abs( x ):
    return _meas_no_checks( np.abs( x.x ), x.dx )


def isnan( x ):
    return np.isnan( x.x )
                        

nan = meas( np.nan, np.nan ) 
