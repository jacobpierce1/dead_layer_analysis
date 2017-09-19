# my includes 
import sql_db_manager
# import jacoblib.jacob_file_parsers as fparsers

## includes 
from mpl_toolkits.mplot3d import axes3d
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import sqlite3
import json



CM_PER_INCH = 2.54 


sources_index = [ 'pu_240', 'cf_249', 'pu_238_centered', 'pu_238_moved', 'pu_238_flat', 'pu_238_tilted' ] 
all_objects_index = sources_index.append( 'detector' )
values_index = [ 'value', 'delta' ]


# these are the distances from the cylinder containing the source to the edge of the 
# container. each element of the list is a single measurement, so the actual measurement
# is the result of averaging the numbers. i am putting all measurements here instead of 
# average value out of principle.
source_data = pd.DataFrame
( 
    { 
        'bottom'   :  [ [0.1570,0.1613], [1.2335,1.2300], [0.6830,0.6920], [0.5725,0.5775],
                        [], []   ],
        
        'top'      :  [ [1.8220,1.8285], [0.0580,0.0615], [0.5840,0.5885], [0.6883,0.6853],
                        [], []   ],
        
        'left'     :  [ [2.0470,2.0355], [1.7343,1.7383], [3.8443,3.8480], [6.9030,6.8973],
                        [], []   ],
        
        'right'    :  [ [6.4445,6.4420], [6.0085,6.0040], [3.8980,3.8965], [0.8380,0.8380],
                        [], []   ],
        
        'diameter' :  [ [1.0040,1.0045,1.0040], [1.7520,1.7510,1.7517], [1.7525,1.7515,1.7530],
                        [], [], [] ],
        
        'height'   :  [ [3.0020,3.0030,3.0025], [3.0055,3.0070,3.0055], [3.0015,3.0070, 3.0055],
                        [], [], [] ], 
        
        'wafer'    :  [ [0.0320,0.0326,0.0300], [0.052,0.0515,0.0530], [0.051,0.050,0.049],
                        [], [], [] ]
    },
    index = sources_index
)

source_data_uncertainty = 0.005




# all in cm
detector_data = pd.Series(
    {
        'x_offset' : [5.0, 5.0 ],
        'y_offset' : [ 0.1 ],
        'height' : [0.1]
    }
)



# all in cm
enclosure_data = pd.Series(
    {
        'total_z' : [ 100.0 ],
        'top_offset' : [ 5.0 ],
        'bottom_offset' : [ 5.0 ],
        'total_x' : [ 50.0 ],
        'total_y' : [ 50.0 ]
    }
)




# fill in the redundant pu_238 entries 
def fill_redundant_source_data( source_data ):
    reference_source = 'pu_238_centered'

    redundant_sources = ['pu_238_moved', 'pu_238_flat' ] 
    redundant_cols = [ 'diameter', 'height', 'wafer' ]
            
    for col in redundant_cols:
        source_data[col].loc[redundant_sources] = source_data[col].loc[reference_source]

    #  more to be dnoe here ...



    


# return angle between source and detector pixel (i,j) given all relevant parameters. 
# this angle and the uncertainty in it are returned as a tuple (angle, delta_anle ).
# x, y, z must all be a series or dict with entries for 'value' and 'delta'
def get_costheta( x, y, z ):
    value= z / np.sqrt( x['value']**2.0 + y['value']**2.0 + z['value']**2.0 )
    delta = get_costheta_delta( x['value'], x['delta'], y['value'], y['delta'], z['value'], z['delta'] )
    return pd.Series( { 'value' : value, 'delta' : delta } )
    
def get_costheta_delta( x,dx, y,dy, z,dz ) :
    return ( ( z**2 * ( (x*dx)**2 + (y*dy)**2 )  +
                    ( x**2 + y**2 )**2 * dz**2 )  /
                        (x**2 + y**2 + z**2 )**3 )




# # get the coords detector referenced to 0
# det_coords = pd.Series(
#     [
#         ( enclosure_data['detector_xy'][0],
#           enclosure_data['detector_xy'][1],
#           enclosure_data['total_height'],
#           - enclosure_data['top_offset'],
#           - enclosure_data['detector_height']
#         ),
        
#         ( 0,0,0 ),
#     ],
#     index = values_index
# )


all_coords = pd.DataFrame(
    columns = values_index,
    index = all_objects_index 
)




# use the measurements in source_data to obtain the coordinates of each source and a
# edge of the detector. 
def populate_source_coords( source_coords ):

    # xyz is 3 tuple storing coordinates, xyz_delta stores uncertainties.
    xyz = np.empty(3)
    xyz_delta = np.empty(3)

    # first handle the detector separately
    ## TODO URGENT
    all_coords[values_index[0]].loc['detector'] = ( 1.0, 2.0, 3.0 ) # xyz
    all_coords[values_index[1]].loc['detector'] = ( 0.1, 0.2, 0.3 )  # xyz_delta* CM_PER_IN
    

    # now populate all the source indices
    for source in sources_index:
        
        # first get x and y from the diameters / distance to edges.
        configurations = [ ['left','right','diameter'], ['top','bottom','diameter'] ]
        for i in range(len(configurations)):

            # these coords are the sum of the diameter plus distance to walls over 2 
            xyz[i] = np.sum( [ np.mean( source_data[col].loc[source] )
                               for col in configurations[i] ] ) / 2
            
            # for uncertainty: square root of sums squared. the uncertainty of each mean
            # is source_data_uncertainty / sqrt(num_measurements).
            xyz_delta[i] = source_data_uncertainty * np.sqrt(
                np.sum( [ 1 / len( source_data[col].loc[source] ) )
                    for col in configurations[i] ] ) ) / 2

        # add the coordinates and uncertainties to the source_data. convert from in to cm.
        all_coords[values_index[0]].loc[source] = xyz * CM_PER_IN
        all_coords[values_index[1]].loc[source] = xyz_delta* CM_PER_IN

        

    
# in: 3-tuple of x,y,z coordinates
# out: rotated vector about z axis by theta
def rotate_z( theta, x ):
    return np.dot( np.array( [ np.cos(theta), 0-np.sin(theta), 0 ],
                       [ np.sin(theta),   np.cos(theta), 0 ],
                             [ 0, 0, 1 ]  ),  x  )

def rotate_x( phi, x ):
    return np.dot( np.array( [1,0,0],
                             [ 0, np.cos(phi), 0-np.sin(phi) ],
                             [ 0, np.sin(phi), np.cos(phi) ] ),  x  )
    
# unused
def rotate_y( theta, x ):
    pass 


# all uncertainty calculations with the rotation matrices use this formula for different
# permutations of the 3-tuple vector x. not meant to be used more generally outside these
# uncertainty calculations. x assumed to have 3 'value' and 'delta' entries.
def rotated_delta_entry( theta, x ):
    return ( ( x['value'][0] * np.cos(theta['value']) * theta['delta'] )**2.0  +
             ( np.cos(theta['value'] ) 
             x['value'][1] * np.sin



# in: 3-tuple coords of detector, coords of a particular source, and theta / phi 2-tuple of the
# source if it is rotated (only used for source deadlayer calculation )
# out: pd.Series containing 32x32 matrix of penetration angle for each
# detector pixel and another for the angle of penetration of the source layer.
# theta is measured from 0 to pi with 0 on the z axis, phi from 0 to 2pi with 0 on the x axis,
# same as the direction of 'right'

def get_costheta_grid( det_coords, source_coords, theta_phi=(0,0) ):

    grid = pd.Series(
        {
            'det_costheta_values' : np.empty(32,32),
            'det_costheta_deltas' : np.empty(32,32),
            'source_costheta_values' : np.empty(32,32),
            'source_costheta_deltas' : np.empty(32,32)
        }
    )

    # keep shifting by 0.1cm to get next coord 
    for i in range(32):
        for j in range(32):

            displacement_value = det_coords['value'] + ( i*0.1, j*0.1, 0 ) - source_coords['value']
            displacement_delta = np.sqrt( det_coords['delta']**2 + source_coords['delta']**2 ) 

            # get angle rel to detector pixel, put in the arrays
            costheta = get_costheta( pd.Series( { 'value' : displacement_value,
                                                  'delta' : displacemnet_delta } ) )
            grid['det_costheta_values'][i][j] = costheta['value']
            grid['det_costheta_deltas'][i][j] = costheta['delta']

            # now find penetration angle rel to source deadlayer.
            # see function description for definition of theta and phi

            # avoid the computation if we can. if it is (0,0) then the planes of the detector and
            # source deadlayer are parallel. in this case we have that the angle through the detector
            # is the same as through the source deadlayer. this is the case most of the time.
            if theta_phi is (0,0):
                grid['source_costheta_values'][i][j] = costheta['value']
                grid['source_costheta_deltas'][i][j] = costheta['delta']

            # if not, then we rotate the displacement vector by 0-theta about the z axis, then
            # 0-phi about the x axis, then take costheta. this is because rotating by the negative
            # angles is equivalent to rotating the source by positive angles, which is the
            # definition of the theta/phi angles.
            rotated_displacement_value = rotate_x( theta_phi[1],
                                                   rotate_z( theta_phi[0], displacement_value ) ) 
