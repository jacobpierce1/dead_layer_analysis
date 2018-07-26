import numpy as np
import jutils 



det_center_x = 1.08 * 25.4
# det_center_x = -1.08 * 25.4

det_center_y = 0.2 * 25.4 # plus minus of this

det_center_z = 4.155 * 25.4


# first are the strips up to which source1 illuminates (not inclusively)
# and second are strips above which source2 illuminates. e.g. source 1
# for det 2 illuminates fstrips 0, ..., 8.

# displacements = np.array( [ [ -4.08990552e+00 -2.26874715e+01  9.79663739e+01 ],
#                             [ -1.09875774e+01 -4.06935820e+01, 9.21738382e+01 ] ] )






# key: group, det, sourcenum, coords
displacements =  [
    # trapped ion cloud 
    [[ 32.0, -32.0, 60.0 ]] * 4,
    # group 1 
    [
        # det 1 
        [
            # left source
            [ -5.76844528e+00, -2.57996212e+01, 9.57466449e+01 ],
            [  3.68404794e-01, -3.71099292e+01,  9.57466449e+01 ]  
        ],
        [
            [ -6.27236723e+00, -2.54154184e+01, 9.31458402e+01 ],
            [ -3.39243777e+00, -3.08870082e+01,  9.31458402e+01 ],
        ],
        [
            [ -2.14729892e+00, -2.80700556e+01,  8.70656498e+01 ], # bad data 
            [ -1.41022536e+01, -4.02107244e+01,  8.70656498e+01 ]
        ],
        [
            [ -5.47707985e+00, -2.29581393e+01, 9.57494496e+01 ],
            [ -8.88687108e+00, -4.22352995e+01,  9.57494496e+01 ]
        ]
    ],
    # group 2 
    [
        # det 1 
        [
            [ -6.04991268e+00, -2.88246923e+01, 9.80023281e+01 ],
            [  -1.00809834e+01, -2.85989816e+01,  9.80023281e+01 ]
        ],
        [
            [  -2.14729892e+00, -2.80700556e+01, 9.05313381e+01 ],
            [ -1.00684580e+01, -3.54901722e+01,  9.05313381e+01 ]
        ],
        [
            [  -7.10437249e+00, -2.39445899e+01,  9.38582831e+01 ], # bad data 
            [  -6.76254885e+00, -3.09655228e+01,  9.38582831e+01 ]
        ],
        [
            [ -3.10025524e+00, -2.81885239e+01, 9.06208294e+01 ],
            [ -8.61585443e+00, -3.18512124e+01, 9.06208294e+01 ]
        ]
    ]
] 



illuminated_fstrips = [
    [],
    # group
    [
        # det / source
        [10,23],
        [11,25],
        [6,20],
        [10,18]
    ],
    [
        [5,16],
        [10,21],
        [9,24],
        [6,23]
    ]
] 

    


# i = group (cm or gd )
# d = detnum
# x = fstrip
# y = bstrip

def get_secant( i, j, d, x, y ) :

    # be 8 continuum
    if i == 0 :
        ret = 1 
        displacement = displacements[i][d] + 2.0 * np.array( [ -y, x, 0 ] )

    else : 
        # if x < illuminated_fstrips[i][d][0] :
        #     k = 1

        # elif x > illuminated_fstrips[i][d][1] :
        #     k = 0

        # else :
        #     return np.nan
        
        
        displacement = displacements[i][d][j] + 2.0 * np.array( [ y, x, 0 ] )

        ret = 0 
        
        if j == 0 and x < illuminated_fstrips[i][d][0] :
            ret = 1
        elif j == 1 and x > illuminated_fstrips[i][d][1] :
            ret = 1
            
    if ret :
        return  np.linalg.norm( displacement ) / np.abs( displacement[2] ) 
    else :
        return np.nan
    



def get_secant_matrices(  ) : 

    # if not separate_sources : 
    
    #     det_sectheta = np.zeros( (3, 4, 32,32))

    #     for i in range( 3 ) : 
    #         for d in range(4) :
    #             for x in range(32) :
    #                 for y in range(32) :
    #                     det_sectheta[ i, d, x, y ] = get_secant( i, d, x, y ) 

    det_sectheta = np.zeros( ( 3,2,4,32,32 ) )
    for i in range(3):
        for j in range(2):
            for d in range(4):
                for x in range(32):
                    for y in range(32) :
                        det_sectheta[ i,j,d,x,y ] = get_secant( i,j,d,x,y )
                        
    return det_sectheta
