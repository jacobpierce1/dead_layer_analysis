import numpy as np

x = np.zeros( (32,32 ), dtype = 'float64' )

print( x )
print( x.shape )
print( x.dtype ) 
print( x.dtype.byteorder ) 

np.savetxt( 'test.csv', x, delimiter = ',' )
