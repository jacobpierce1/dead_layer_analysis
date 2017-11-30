
import numpy as np 
import deadlayer_helpers.geometry as geom 


z = np.array( [1,1] )
print( geom.rotate_2d( z, np.pi / 4 ) )


x = np.array( [0, 0, 1 ] )

print( geom.rotate_3d( 2, 0, x ) )


y = np.array( [ 30, 0, 0] )
print( geom.rotate_3d( 2, np.pi / 4, y ) )
