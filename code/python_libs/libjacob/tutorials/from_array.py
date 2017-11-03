import libjacob.jmeas as meas

import numpy as np

x = np.array( [ [ [ meas.meas(2,3)] * 3] * 2] * 3 )

print( 'before from_array: ' )
print( x )
print( '\n' )

y = meas.meas.from_array( x  )

print( 'after from_array: ' )
print( y ) 
