import numpy as np

class double_array( object ) :

    # x and y assumed to be numpy arrays of
    # same shape
    
    def __init__( self, x, y ):
        self.x = x
        self.y = y


    def __add__( self, z ) :
        return double_array( self.x + z, self.y )


    # handle RHS addition the same way
    def __radd__( self, z ) :
        return self + z
    

    def __repr__( self ) :
        return '%s %s' % ( str( self.x ), str( self.y ) )

    
    def __getitem__( self, key ):
        return double_array( self.x[key], self.y[key] )

    
    # value must be a meas
    def __setitem__( self, key, value ):
        self.x[key] = value.x
        self.y[key] = value.y
        return self

    def __len__( self ) :
        return len( self.x ) 

    
darray = double_array( np.arange(4) ,
                       np.zeros(4) )

print( darray )

# print( darray + 6 )

# print( 6 + darray )

x = np.array( [ 0.0, 2.0] )

print( darray + x[0] )
print( x[0] + darray )
print( np.float64(0) + darray ) 
