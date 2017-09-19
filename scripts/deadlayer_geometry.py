# the purpose of this module is to construct several global variables that are extremely
# useful, such as grids of cosine theta values for all the detector pixels.
# the goal is that geometry issues and calculations will never need to be considered outside the
# module.


# my includes 
import sql_db_manager
import error 


## includes 
from mpl_toolkits.mplot3d import axes3d
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import sqlite3
import json


CM_PER_INCH = 2.54 


sources = [ 'pu_240', 'cf_249', 'pu_238_centered', 'pu_238_moved', 'pu_238_flat', 'pu_238_tilted' ]
all_objects = sources + [ 'detector' ]

# indices to be used for constructing DataFrames
sources_index = pd.Index( sources )
all_objects_index = pd.Index( all_objects )



# uncertainty for all measurements. note that this is 

source_data_delta = 0.005



# these are the distances from the cylinder containing the source to the edge of the 
# container. each element of the list is a single measurement, so the actual measurement
# is the result of averaging the numbers. i am putting all measurements here instead of 
# average value out of principle.

source_data = pd.DataFrame.from_dict( 
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
    } 
)


# measurements used to construct tilt angle of pu_238_tilted
pu_238_tilted_data = pd.Series(
    {
        'lower_height' : [1,1,1],
        'upper_height' : [2,2,2]
    }
)




# measurements relevant to the detector
# all in cm
detector_data = pd.Series(
    {
        'x_offset' : [0.266,0.273,0.263],  # left, or right. todo: figure out.
        'y_offset' : [0.1555, 0.1575],
        'height' : [0.1]
    }
)



# all in cm
enclosure_data = pd.Series(
    {
        'top_offset' : [ 1.0 ],
        'bottom_offset' : [ 1.0 ],
        'total_x' : [ 151.696 ],
        'total_y' : [ 50.0 ],
        'total_z' : [ 100.0 ]
    }
)




# fill in the redundant pu_238 entries 
def fill_redundant_source_data( source_data ):

    # these sources share the same diam, height, and wafer because they are the same, just moved.
    reference_source = 'pu_238_centered'
    redundant_sources = ['pu_238_moved', 'pu_238_flat', 'pu_238_tilted' ] 
    redundant_cols = [ 'diameter', 'height', 'wafer' ]
    for col in redundant_cols:
        for redundant_source in redundant_sources:
            source_data.set_value( redundant_source, col, source_data.loc[ reference_source, col ] )
            
    # in addition these two share the same position with the centered source:
    reference_source = 'pu_238_centered'
    redundant_sources = [ 'pu_238_flat', 'pu_238_tilted' ]
    redundant_cols = [ 'top', 'bottom', 'left', 'right' ]
    for col in redundant_cols:
        for redundant_source in redundant_sources:
            source_data.set_value( redundant_source, col, source_data.loc[ reference_source, col ] )
        
    


    


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





# use the measurements in source_data to obtain the coordinates of each source and a
# edge of the detector. 
def populate_all_coords( source_coords, all_coords ):

    # # xyz is 3 tuple storing coordinates, xyz_delta stores uncertainties.
    # xyz = np.empty(3)
    # xyz_delta = np.empty(3)

    # # first handle the detector separately
    # ## TODO URGENT
    # todo_det_coords =  ( 1.0, 2.0, 3.0 )
    # all_coords.set_value( 'detector', (  # xyz
    # # all_coords[values_index[1]].loc['detector'] = ( 0.1, 0.2, 0.3 )  # xyz_delta* CM_PER_IN
    


    
    # now populate all the source indices
    for source in sources_index:
        
        # first get x and y from the diameters / distance to edges.
        #  configurations = [ ['left','right','diameter'], ['top','bottom','diameter'] ]
        # for i in range(len(configurations)):
        
        # start off only using the left and bottom  measurements. TODO: use value of total_x
        x, y = [ error.esum( error.emean( source_data.loc[ source, col ], source_data_delta ),
                             error.emean( np.asarray( source_data.loc[ source, 'diameter' ] )  / 2,
                                          source_data_delta ) )
                 for col in [ 'left', 'bottom' ] ]


        # for the top measurement, reference to the bottom of the enclosure.
        z = error.esum_n( error.emean( enclosure_data['bottom_offset'], source_data_delta ),
                          error.emean( source_data.loc[ source, 'height' ],
                                       source_data_delta )
                          error.emean( soure_data.loc[ source, 'wafer' ], source_data_delta )

        xyz = np.array( [ var['value'] for var in [x,y,z] ] ) 
        xyz_delta =  np.array( [ var['delta'] for var in [x,y,z] ] )
        m = error.measurement( xyz, xyz_delta )
        all_coords.loc[source] = m * CM_PER_INCH
                                       
        # # these coords are the sum of the diameter plus distance to walls over 2 
        # xyz[i] = np.sum( [ np.mean( source_data[col].loc[source] )
        #                    for col in configurations[i] ] ) / 2
        
        # # for uncertainty: square root of sums squared. the uncertainty of each mean
        # # is source_data_uncertainty / sqrt(num_measurements).
        # xyz_delta[i] = source_data_uncertainty * np.sqrt(
        #     np.sum( [ 1 / len( source_data[col].loc[source] ) 
        #         for col in configurations[i] ] ) ) / 2

        # add the coordinates and uncertainties to the source_data. convert from in to cm.
        # all_coords[values_index[0]].loc[source] = xyz * CM_PER_INCH
        # all_coords[values_index[1]].loc[source] = xyz_delta* CM_PER_INCH

        



    
# in: 3-tuple of x,y,z coordinates
# out: rotated vector about z axis by theta
def rotate_z( theta, x ):
    value = np.dot( np.array( [ np.cos(theta['value']), 0-np.sin(theta['value']), 0 ],
                       [ np.sin(theta['value']),   np.cos(theta['value']), 0 ],
                             [ 0, 0, 1 ]  ),  x  )
    delta = np.empty(3)
    delta[0:1] = rotation_matrix_delta( theta, x )
    delta[2] = x['delta'][2]
    return pd.Series( [ value, delta ], values_index )
                                             
                        
    
def rotate_x( phi, x ):
    value = np.dot( np.array( [1,0,0],
                             [ 0, np.cos(phi['value']), 0-np.sin(phi['value']) ],
                             [ 0, np.sin(phi['value']), np.cos(phi['value']) ] ),  x  )
    delta = np.empty(3)
    delta[0] = x['delta'][0]
    delta[1:2] = rotation_matrix_delta( theta, x[1:2] )
    return pd.Series( [ value, delta ], values_index )

    
# unused
def rotate_y( theta, x ):
    pass 


# all uncertainty calculations with the rotation matrices use this formula for different
# permutations of the 3-tuple vector x. not meant to be used more generally outside these
# uncertainty calculations. x assumed to have 3 'value' and 'delta' entries.
def rotated_delta_entry( theta, x ):
    return np.sqrt(
        ( ( x['value'][0] * np.cos(theta['value']) + x['value'][1] * np.sin(theta['value']) ) * theta['delta'] ) **2 + \
        ( x['value'][0] * np.cos(theta['value']))**2 + \
        ( x['value'][1] * np.sin(theta['value']))**2 
)



# return the uncertainty for a standard 2x2 rotation matrix. since the 2x2 rotation matrix
# is embedded in the higher order versions, this can be used to construct their uncertainties.
def rotation_matrix_delta( theta, x ):
    return np.array( [  rotated_delta_entry( theta, x ),
                        rotated_delta_entry( theta, x[1,0,2] * np.array( -1,1,1) )
                        ] )


                        
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
            det_costheta = get_costheta( pd.Series( { 'value' : displacement_value,
                                                  'delta' : displacemnet_delta } ) )
            grid['det_costheta_values'][i][j] = det_costheta['value']
            grid['det_costheta_deltas'][i][j] = det_costheta['delta']

            # now find penetration angle rel to source deadlayer.
            # see function description for definition of theta and phi

            if theta_phi is (0,0):

                # avoid the computation if we can. if it is (0,0) then the planes of the detector and
                # source deadlayer are parallel. in this case we have that the angle through the detector
                # is the same as through the source deadlayer. this is the case most of the time.

                grid['source_costheta_values'][i][j] = det_costheta['value']
                grid['source_costheta_deltas'][i][j] = det_costheta['delta']

                
            else:
                
                # if not, then we rotate the displacement vector by 0-theta about the z axis, then
                # 0-phi about the x axis, then take costheta. this is because rotating by the negative
                # angles is equivalent to rotating the source by positive angles, which is the
                # definition of the theta/phi angles. note that this is independent of the other
                # det_costheta value.

                rotated_displacement = rotate_x( 0-theta_phi[1],
                                                 rotate_z( 0-theta_phi[0], displacement_value ) ) 
                source_costheta = get_costheta( rotated_displacement )
                grid['source_costheta_values'][i][j] = source_costheta['value']
                grid['source_costheta_deltas'][i][j] = source_costheta['delta']

    return grid








########################################################################
# construct relevant variables that will be accessed from this module.


# inefficient but it works.
source_data = source_data.set_index( sources_index ) 
fill_redundant_source_data( source_data )
# print source_data

# declare dataframe to store all coordinates, which are pd.Series of two 3-tuples, one for
# value and one for delta.
all_coords = pd.DataFrame( columns = error.values_index, index = all_objects_index )
# all_coords = pd.DataFrame( columns=error.values, index = all_objects_index )
print all_coords

populate_all_coords( source_data, all_coords ) 
print all_coords 
