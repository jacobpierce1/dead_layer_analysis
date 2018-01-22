# status: unresolved

import libjacob.jmeas as meas
import numpy as np

x = meas.meas( np.arange(4), np.zeros(4) )

y = np.array( [ 0.0, 1.0 ] )

# print( 1.0 - x )

# print( 1 - x )

# print( -x )


print( 'performing addition' ) 
print( x + y[0] )


print( 'performing addition' ) 
print( np.asscalar( y[0] ) + x )


print( 'performing addition' ) 
print( y[0] + x )


print( 'performing addition' ) 
print( 0.0 + x ) 


print( type( y[0] ) )
print( type( x ) ) 
