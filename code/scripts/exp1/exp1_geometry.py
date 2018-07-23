import numpy as np


def get_secant_matrices() :


    displacements = {
        'pu_240' : np.array( [ 8.59875408e+01, -9.86065654e+00,  6.00662818e+01 ] ),
        'cf_249' : np.array( [ 8.48117806e+01, -4.96340923e+01,  5.76559514e+01 ] ),
        'pu_238_centered' : np.array( [ 3.295e+01, -3.145e+01,  5.788e+01 ] ),
        'pu_238_flat' : np.array( [ 3.46898841e+01, -3.16480392e+01,  5.89136921e+01 ] ),
        'pu_238_moved' : np.array( [ -5.25014590e+01, -2.82241896e+01,  6.31551817e+01 ] ),
        'pu_238_angled' : np.array( [  3.60180932e+01, -3.06134161e+01,  5.89913745e+01 ] )
    }

    secant_matrices = {
        'centered' : np.zeros((3,1,32,32)),
        'flat' : np.zeros((3,1,32,32)),
        'moved' : np.zeros((3,1,32,32)),
        'angled' : np.zeros((3,1,32,32)),
        'det3_moved' : np.zeros((3,1,32,32)),
        'det3_flat' : np.zeros((3,1,32,32))
    }

    pu240 = np.zeros((32,32))
    cf249 = np.zeros((32,32))

    for x in range(32) :
        for y in range(32) :
            pu240_displacement = displacements[ 'pu_240' ] + 2.0 * np.array( [ -y, x, 0 ] ) 
            pu240[x,y] = np.linalg.norm( pu240_displacement ) / pu240_displacement[2]
            cf249_displacement = displacements[ 'cf_249' ] + 2.0 * np.array( [ -y, x, 0 ] ) 
            cf249[x,y] = np.linalg.norm( cf249_displacement ) / cf249_displacement[2]
            
    for name in secant_matrices.keys() :
        secant_matrices[ name ][0,0] = pu240
        secant_matrices[ name ][2,0] = cf249
        tmp = secant_matrices[ name ][1,0]

        if name in [ 'centered', 'angled', 'flat', 'moved' ] : 
            pu238_id = 'pu_238_%s' % name
        else :
            tmp_name = name.split('_')[1]
            pu238_id = 'pu_238_%s' % tmp_name
            
        for x in range(32) :
            for y in range(32) :
                pu238_displacement = ( displacements[ pu238_id ]
                                       + 2.0 * np.array( [ -y, x, 0 ] ) )
                tmp[x,y] = np.linalg.norm( pu238_displacement ) / pu238_displacement[2]

    secant_matrices[ 'det1_centered' ] = secant_matrices[ 'centered' ]
    secant_matrices[ 'det1_moved' ] = secant_matrices[ 'moved' ]
                
    return secant_matrices 
