# the purpose of this module is to construct several global variables that are extremely
# useful, such as grids of cosine theta values for all the detector pixels.
# the goal is that geometry issues and calculations will never need to be considered outside the
# module.


# my includes 
import sql_db_manager
import libjacob.meas as meas


## includes 
from mpl_toolkits.mplot3d import axes3d
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import sqlite3
import json


_CM_PER_INCH = 2.54
_MM_PER_INCH = 25.4


sources = [ 'pu_240', 'cf_249', 'pu_238_centered', 'pu_238_moved', 'pu_238_flat', 'pu_238_tilted' ]
all_objects = sources + [ 'detector' ]

source_deadlayers = [ 'si' ] * 6

# indices to be used for constructing DataFrames
sources_index = pd.Index( sources )
all_objects_index = pd.Index( all_objects )



# costheta_labels = [ 'det_costheta_grid', 'source_costheta_grid' ]
# costheta_labels_index = pd.Index( costheta_labels )



# uncertainty for all measurements in inches.
_source_data_delta = 0.005



# these are the distances from the cylinder containing the source to
# the edge of the container. each element of the list is a single
# measurement, so the actual measurement is the result of averaging
# the numbers. i am putting all measurements here instead of average
# value out of principle. empty entries are filled in by
# _fill_redundant_source_data.

def _get_source_data():

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
                            [], [], [] ],
        } 
    )

    return source_data




# measurements used to construct tilt angle of pu_238_tilted
pu_238_tilted_data = pd.Series(
    {
        'lower_height' : [1,1,1],
        'upper_height' : [2,2,2],
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



# pretty sure total_x and total_y are in mm.
enclosure_data = pd.Series(
    {
        'top_offset' : [ 1.0 ],
        'bottom_offset' : [ 1.0 ],
        'total_x' : [ 151.696 ],
        'total_y' : [ 50.0 ],
        'total_z' : [ 5.764 ]  # inches, line 211 of pixelanalysissrc.c
    }
)


###########################################################################################################
#################################### FUNCTIONS FOR POPULATING THE DATAFRAMES ###############################
###########################################################################################################




# fill in the redundant pu_238 entries. the dataframe / series method set_value allows you to add lists into
# the dataframe / series.
def _fill_redundant_source_data( source_data ):

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

            
    



# use the measurements in source_data to obtain the coordinates of each source and a
# edge of the detector. 
def _populate_all_coords( all_coords, source_data ):

    
    # first handle the detector separately
    ## TODO URGENT
    all_coords.loc['detector'] = error.measurement( np.array( [10.0,10.0,10.0]),
                                                    np.array( [0.01,0.01,0.01] ) )

    
    # now populate all the source indices
    for source in sources_index:
        
        # first get x and y from the diameters / distance to edges.
        #  configurations = [ ['left','right','diameter'], ['top','bottom','diameter'] ]
        # for i in range(len(configurations)):
        
        # start off only using the left and bottom  measurements. TODO: use value of total_x
        x, y = [ meas.meas( source_data.loc[ source, col ],
                            _source_data_delta ).mean()
                 + meas.meas( source_data.loc[ source, 'diameter' ] / 2,
                         source_data_delta ) for col in [ 'left', 'bottom' ]  ]

        # for the top measurement, reference to the bottom of the enclosure.
        z = meas.sum( [ meas.meas( enclosure_data['bottom_offset'], _source_data_delta ).mean(),
                        meas.meas( source_data.loc[ source, 'height' ], _source_data_delta ).mean(),                                                                  meas.meas( source_data.loc[ source, 'wafer' ], _source_data_delta ).mean()
        ] )

        xyz = meas.meas.fromlist( [ x, y, z ] )
        all_coords.loc[source] = xyz * _MM_PER_INCH
                                    
        
    


# return angle between source and detector pixel (i,j) given all relevant parameters. 
# this angle and the uncertainty in it are returned as a tuple (angle, delta_anle ).
# x, y, z must all be a series or dict with entries for 'value' and 'delta'
def get_costheta( coords ):

    x, y, z = coords.x
    dx, dy, dz = coords['delta']
    
    value= z / np.sqrt( x**2.0 + y**2.0 + z**2.0 )
    delta = get_costheta_delta( x, dx, y, dy, z, dz )





# def get_costheta_delta( x,dx, y,dy, z,dz ) :
#     return ( ( z**2 * ( (x*dx)**2 + (y*dy)**2 )  +
#                     ( x**2 + y**2 )**2 * dz**2 )  /
#                         (x**2 + y**2 + z**2 )**3 )






    
# in: 3-tuple of x,y,z coordinates
# out: rotated vector about z axis by theta
def rotate_z( theta, x ):

    R = np.array( [ [ np.cos(theta['value']), 0-np.sin(theta['value']), 0 ],
                    [ np.sin(theta['value']),   np.cos(theta['value']), 0 ],
                    [ 0, 0, 1 ] ] )
       
    value = np.dot( R, x['value'] )
    delta = np.empty(3)
    delta[0:2] = rotation_matrix_delta( theta, x )
    delta[2] = x['delta'][2]
    return pd.Series( [ value, delta ], error.values_index )
                                             
                        
    
def rotate_x( phi, x ):
    value = np.dot( np.array( [ [ 1, 0, 0],
                                [ 0, np.cos(phi['value']), 0-np.sin(phi['value']) ],
                                [ 0, np.sin(phi['value']), np.cos(phi['value']) ] ] ),
                    x['value']  )

    delta = np.empty(3)
    delta[0] = x['delta'][0]
    delta[1:3] = rotation_matrix_delta( phi, x )
    return pd.Series( [ value, delta ], error.values_index )

    
# unused
def rotate_y( theta, x ):
    pass 


# all uncertainty calculations with the rotation matrices use this formula for different
# permutations of the 3-tuple vector x. not meant to be used more generally outside these
# uncertainty calculations. x assumed to have 3 'value' and 'delta' entries.
def rotated_delta_entry( theta, x ):
    
    return np.sqrt( ( ( x['value'][0] * np.sin(theta['value'])
                        + x['value'][1] * np.cos(theta['value']) )
                      * theta['delta'] )**2 +
                    + ( x['value'][0] * np.cos(theta['value']))**2 
                    + ( x['value'][1] * np.sin(theta['value']))**2 )



# return the uncertainty for a standard 2x2 rotation matrix. since the 2x2 rotation matrix
# is embedded in the higher order versions, this can be used to construct their uncertainties.
def rotation_matrix_delta( theta, x ):

    # first entry
    y1 = rotated_delta_entry( theta, x )

    # switch two columns in this copy. 
    xcopy = x.copy()
    xcopy_matrix = xcopy.values
    xcopy_matrix[0][0], xcopy_matrix[0][1] = 0-xcopy_matrix[0][1], xcopy_matrix[0][0]
    xcopy_matrix[1][0], xcopy_matrix[1][1] = xcopy_matrix[1][1], xcopy_matrix[1][0]

    # second entry
    y2 = rotated_delta_entry( theta, xcopy )  
    return np.array( [  y1, y2 ] )
                        






def _populate_source_theta_phi( source_theta_phi, source_data, pu_238_tilted_data ):

    # assume all measured (or implicitly measured, e.g. assuming something is flat ) are known to
    # 1 degree.
    angle_data_delta = np.deg2rad( 1.0 / 360 ) 

    
    # get the theta of pu_238_tilted. assuming phi = 0, which mary told me.
    radius = error.emean( error.measurement( source_data.loc[ 'pu_238_tilted', 'diameter' ],
                                             _source_data_delta ) ) / 2.0
    
    height_diff = error.esum( error.emean( pu_238_tilted_data[ 'upper_height' ],
                                           _source_data_delta ),
                              error.emean( pu_238_tilted_data[ 'lower_height' ],
                                           _source_data_delta ) )

    tantheta = error.edivide( height_diff, radius )
    theta_value = np.arctan( tantheta['value'] )

    theta_delta = angle_data_delta  # todo
    
    source_theta_phi.loc[ 'pu_238_tilted' ] = error.measurement( [ theta_value, 0.0 ],
                                                                 [ theta_delta, angle_data_delta ] )
    
    # add in theta_phi values for everything except pu_238_tilted.
    upright_sources = list(sources)
    upright_sources.remove( 'pu_238_tilted' )
    for source in upright_sources:
        source_theta_phi.loc[ source ] = error.measurement( [0,0], 2*[angle_data_delta] )

                                    
                                    




                          
# in: 3-tuple coords of detector, coords of a particular source, and theta / phi 2-tuple of the
# source if it is rotated (only used for source deadlayer calculation )
# out: pd.Series containing 32x32 matrix of penetration angle for each
# detector pixel and another for the angle of penetration of the source layer.
# theta is measured from 0 to pi with 0 on the z axis, phi from 0 to 2pi with 0 on the x axis,
# same as the direction of 'right'

def _populate_costheta_grid( cosine_matrices, all_coords, source_theta_phi ):
       
    det_coords = all_coords.loc['detector']
    
    # loop through all sources
    for sourcenum in range(len(sources)):

        print str(sourcenum) + ' / ' + str(len(sources))

        source_coords = all_coords.loc[ sources[sourcenum] ]
        theta_phi = source_theta_phi.loc[ sources[sourcenum] ]
        
        theta, phi = [ error.measurement( theta_phi['value'][i], theta_phi['delta'][i] )
                       for i in range(2) ]

        # pre-index the current values in order to enhance readability 
        det_costheta_grid_value = cosine_matrices[sourcenum, 0, 0]
        det_costheta_grid_delta = cosine_matrices[sourcenum, 0, 1]
        source_costheta_grid_value = cosine_matrices[sourcenum, 1, 0]
        source_costheta_grid_delta = cosine_matrices[sourcenum, 1, 1]
        
        
        # keep shifting by 1mm to get next coord 
        for i in range(32):
            for j in range(32):
                
                    
                displacement_value = det_coords['value'] + np.array( [i*1.0, j*1.0, 0 ] ) - source_coords['value']
                displacement_delta = np.sqrt( det_coords['delta']**2 + source_coords['delta']**2 ) 
                displacement = error.measurement( displacement_value,
                                                  displacement_delta )
                
                # get angle rel to detector pixel, put in the arrays
                det_costheta = get_costheta( displacement ) 
                det_costheta_grid_value[i][j] = det_costheta['value']
                det_costheta_grid_delta[i][j] = det_costheta['delta']
                
                # now find penetration angle rel to source deadlayer.
                # see function description for definition of theta and phi
                
                # if theta_phi == [0,0]:
                    
                #     # avoid the computation if we can. if it is (0,0) then the planes of the detector and
                #     # source deadlayer are parallel. in this case we have that the angle through the detector
                #     # is the same as through the source deadlayer. this is the case most of the time.
                    
                #     source_costheta_grid_value[i][j] = det_costheta['value']
                #     source_costheta_grid_delta[i][j] = det_costheta['delta']
                    
                    
                if 1:
                    
                    # if not, then we rotate the displacement vector by 0-theta about the z axis, then
                    # 0-phi about the x axis, then take costheta. this is because rotating by the negative
                    # angles is equivalent to rotating the source by positive angles, which is the
                    # definition of the theta/phi angles. note that this is independent of the other
                    # det_costheta value.
                    
                    rotated_displacement = rotate_x( 0-phi,
                                                     rotate_z( 0-theta, displacement ) ) 
                    source_costheta = get_costheta( rotated_displacement )
                    source_costheta_grid_value[i][j] = source_costheta['value']
                    source_costheta_grid_delta[i][j] = source_costheta['delta']
                    
                    
                    
                    
                        
    
                    






######################################################################## MAIN ##########################

# call get_costheta_grid in order to populate the cosine matrices with
# uncertainties.  construct all unnecessary (outside this module)
# intermediaries, such as source coordinates.

# data structure to hold the cosine matrices: ( source ) -> ( detector
# angle, source angle ) -> ( value, delta )

#  i thought for a while about what pandas data structure to use
# before realizing that they all kind of suck. i am just going to
# declare an array here that has the labels for each dimension. if you
# really need to look up by a particular label, just figure out the
# index and do it. good riddance. in retrospect i probably will never
# use pandas again and always take the route of using an array with
# labels stored elsewhere.

# use debug = 1 to return an array of the same type but all entries are 0.5

def get_cosine_matrices( debug=0 ):


    if debug:
        return np.ones( ( len(sources), 2, 2, 32, 32 ) ) / 2.0
    
    print 'INFO: constructing costheta grid...'
    cosine_matrices_labels = [ sources, ['detector', 'source'], ['value', 'delta'] ]
    cosine_matrices = np.empty( ( len(sources), 2, 2, 32, 32 ) )

    # inefficient but it works.
    source_data = _get_source_data() 
    source_data = source_data.set_index( sources_index )
    _fill_redundant_source_data( source_data )

    
    # this shall be populated with the theta_phi angles for all detectors, which are 0
    # for all but the pu_238_tilted
    source_theta_phi = pd.DataFrame( columns = error.values_index, index = sources_index )
    _populate_source_theta_phi( source_theta_phi, source_data, pu_238_tilted_data )
    
        
    # declare dataframe to store all coordinates, which are pd.Series of two 3-tuples, one for
    # value and one for delta.
    all_coords = pd.DataFrame( columns = error.values_index, index = all_objects_index )
    _populate_all_coords( all_coords, source_data ) 


    # get each array of values / uncertainties and add to the grid.
    _populate_costheta_grid( cosine_matrices, all_coords, source_theta_phi )
    
    return cosine_matrices
    



