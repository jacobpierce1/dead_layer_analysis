import numpy as np
import jutils 



det_center_x = 1.08 * 25.4
# det_center_x = -1.08 * 25.4

det_center_y = 0.2 * 25.4 # plus minus of this

det_center_z = 4.155 * 25.4


# first are the strips up to which source1 illuminates (not inclusively)
# and second are strips above which source2 illuminates. e.g. source 1
# for det 2 illuminates fstrips 0, ..., 8.

illuminated_fstrips = [ [ 9, 9, 6, 8 ],
                        [ 24, 25, 20, 20 ] ]

# displacements = np.array( [ [ -4.08990552e+00 -2.26874715e+01  9.79663739e+01 ],
#                             [ -1.09875774e+01 -4.06935820e+01, 9.21738382e+01 ] ] )

def get_secant( detnum, x, y ) :

    fstrip = x + 1
    bstrip = y + 1

    theta1 = np.arctan2( np.sqrt( (det_center_x + ( bstrip - 16.5) * 2 ) ** 2 
                                     + ( det_center_y + ( fstrip - 16.5) * 2 ) ** 2 ),
                                     det_center_z )

    theta2 = np.arctan2( np.sqrt( ( det_center_x + (bstrip - 16.5 ) * 2) ** 2
                                     + ( -det_center_y + ( fstrip - 16.5 ) * 2) ** 2),
                                     det_center_z )

    # theta1 = 
    
    # if x < illuminated_fstrips[ 0 ][ detnum ] :
    #     return 1 / np.cos( theta2 )
        
    # elif x > illuminated_fstrips[ 1 ][ detnum ] :
    #     return 1 / np.cos( theta1 ) 

    else :
        return np.nan

    

# # test secant ( remove )
# theta1 = np.arctan2( np.sqrt( (det_center_x + ( np.arange(32) - 16.5) * 2 ) ** 2 
#                               + ( det_center_y + ( 0 - 16.5) * 2 ) ** 2 ),
#                      det_center_z )

# theta2 = np.arctan2( np.sqrt( (det_center_x + ( np.arange(32) - 16.5) * 2 ) ** 2 
#                               + ( - det_center_y + ( 0 - 16.5) * 2 ) ** 2 ),
#                      det_center_z )



# print( 'theta1 :', theta1 ) 

# print( 'theta2 :', theta2 ) 

# sys.exit(0);


def get_secant_matrices() : 
    
    det_sectheta = np.zeros( (2, 4, 32,32))

    for i in range( 2 ) : 
        for d in range(4) :
            for x in range(32) :
                for y in range(32) :
                    det_sectheta[ i, d, x, y ] = get_secant(
                        d, x, y ) 

    return det_sectheta
