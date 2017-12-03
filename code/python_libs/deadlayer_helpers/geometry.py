# the purpose of this module is to construct several global variables
# that are extremely useful, such as grids of cosine theta values for
# all the detector pixels.  the goal is that geometry issues and
# calculations will never need to be considered outside the module.


# my includes 
import deadlayer_helpers.sql_db_manager as db
import libjacob.jmeas as meas


## includes 
from mpl_toolkits.mplot3d import axes3d
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt


_CM_PER_INCH = 2.54
_MM_PER_INCH = 25.4


sources = [ 'pu_240', 'cf_249', 'pu_238_centered',
            'pu_238_moved', 'pu_238_flat', 'pu_238_angled' ]

all_objects = sources + [ 'detector' ]



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






# measurements used to construct tilt angle of pu_238_angled
pu_238_angled_data = pd.Series(
    {
        'lower_height' : 2.14,  # inches 
        'upper_height' : 3.88,
    }
)




# measurements relevant to the detector
# all in cm
detector_data = pd.Series(
    {
        'x_offset' : [ 2.6630, 2.6685 ], #2.6843, 
        #'y_offset' : [0.2663, 0.2668 ], 
        'y_offset' : [ 0.1555, 0.1575 ],
        'height' : [1.], # mm
        'total_width' : 64 # mm, this is 32 strips * 2 mm in both dimensions 
    }
)



# pretty sure total_x and total_y are in mm.
# total_z and top / bottom offset are obtained from pixelanalysisnew.c: lines 210 - 212,
# e.g.         Hcf249=5.764-0.1260-0.267-3.006-0.0522; //inches
enclosure_data = pd.Series(
    {
        'top_offset' : [ 0.1260 ],
        'bottom_offset' : [ 0.267 ],
        #        'total_x' : [ 151.696 ],  # pretty sure this is wrong. unused in analysis
        #        'total_y' : [ 50.0 ], # pretty sure this is wrong. unused in analysis.
        'total_z' : [ 5.764 ]  # inches, line 211 of pixelanalysissrc.c
    }
)



# ambiguity: which is x or y? TODO
ceramic_data = pd.Series(
    {
        'total_x' : 106.0,    # mm
        'total_y' : 80.0 # mm
    }
)




# extra_source_offset = pd.Series(
#     {
#         'x' : 1.2,
#         'y' : 7.475
#     }
# )




##############################################################################################
#################################### FUNCTIONS FOR POPULATING THE DATAFRAMES #################
##############################################################################################


# fill in the redundant pu_238 entries. the dataframe / series method
# set_value allows you to add lists into the dataframe / series.

def _fill_redundant_source_data( source_data ):

    # these sources share the same diam, height, and wafer because they are the same, just moved.
    reference_source = 'pu_238_centered'
    redundant_sources = ['pu_238_moved', 'pu_238_flat' ] 
    redundant_cols = [ 'diameter', 'height', 'wafer' ]
    for col in redundant_cols:
        for redundant_source in redundant_sources:
            source_data.set_value( redundant_source, col,
                                   source_data.loc[ reference_source, col ] )
            
    # in addition these two share the same position with the centered source:
    reference_source = 'pu_238_centered'
    redundant_sources = [ 'pu_238_flat', 'pu_238_angled' ]
    redundant_cols = [ 'top', 'bottom', 'left', 'right' ]
    for col in redundant_cols:
        for redundant_source in redundant_sources:
            source_data.set_value( redundant_source, col,
                                   source_data.loc[ reference_source, col ] )

    # final redundant parameter
    source_data.set_value( 'pu_238_angled' , 'diameter',
                           source_data.loc[ 'pu_238_centered', 'diameter' ] )
    
    

    

# use the measurements in source_data to obtain the coordinates of
# each source and a edge of the detector.

_DEBUG_COORDS = 0

def _populate_all_coords( all_coords, source_data ):

    # first handle the detector separately

    # option 1 : "ceramic L/R" to plate edge in the logbook means
    # measured from the left
    # results:

    det_coords = meas.meas.from_list( [ meas.meas( detector_data['x_offset'],
                                                   _source_data_delta ).mean(),
                                        meas.meas( detector_data['y_offset'],
                                                   _source_data_delta ).mean(),
                                        meas.meas( enclosure_data['total_z'],
                                                   _source_data_delta ).mean()
                                        - meas.meas( enclosure_data['top_offset'],
                                                   _source_data_delta ).mean()
                                        - meas.meas( enclosure_data['bottom_offset'],
                                                   _source_data_delta ).mean() ] )


    x_measurement_inverted = 1
    y_measurement_inverted = 0

    if x_measurement_inverted :
        total_x = ( np.mean( source_data.loc[ 'pu_240', 'left' ] )
                    + np.mean( source_data.loc[ 'pu_240', 'diameter' ] )
                    + np.mean( source_data.loc[ 'pu_240', 'right' ] ) )

        det_coords[0] = total_x - det_coords[0]

        
    if y_measurement_inverted :
        total_y = ( np.mean( source_data.loc[ 'pu_240', 'bottom' ] )
                    + np.mean( source_data.loc[ 'pu_240', 'diameter' ] )
                    + np.mean( source_data.loc[ 'pu_240', 'top' ] ) )
        
        det_coords[1] = total_y - det_coords[1]


        
    det_coords *= _MM_PER_INCH


    if x_measurement_inverted :
        shift = - ( ceramic_data[ 'total_x' ]
                    - detector_data[ 'total_width' ] ) / 2 
        det_coords += np.array( [ shift, 0, 0 ] )

    else:
        shift = ( ceramic_data[ 'total_x' ]
                  + detector_data[ 'total_width' ] ) / 2 
        det_coords += np.array( [ shift, 0, 0 ] )
        

    if y_measurement_inverted :
        shift = - detector_data[ 'total_width' ]
        det_coords += np.array( [0, shift, 0 ] )

    else:
        shift = ( ceramic_data[ 'total_y' ]
                  - detector_data[ 'total_width' ] )
        det_coords += np.array( [0, shift, 0 ] )


    # the 1 mm additions center the pixel. 
    det_coords += np.array( [ -1, 1, 0 ] ) 
        
            
    all_coords.loc['detector'] = det_coords

    

        
    # now populate all the source indices
    for source in sources_index:
        
        # first get x and y from the diameters / distance to edges.
        #  configurations = [ ['left','right','diameter'], ['top','bottom','diameter'] ]
        # for i in range(len(configurations)):
        
        # start off only using the left and bottom  measurements. TODO: use value of total_x
        x, y = [ meas.meas( source_data.loc[ source, col ],
                            _source_data_delta ).mean()
                 + meas.meas( source_data.loc[ source, 'diameter' ],
                              _source_data_delta ).mean() / 2 
                 for col in [ 'left', 'bottom' ]  ]

        # for the top measurement, reference to the bottom of the enclosure.

        if source != 'pu_238_angled' : 
            z = meas.sum( [ meas.meas( source_data.loc[ source, 'height' ],
                                       _source_data_delta ).mean(),
                            meas.meas( source_data.loc[ source, 'wafer' ],
                                       _source_data_delta ).mean() ] )

        else:
            z = ( meas.meas( pu_238_angled_data.loc[ 'upper_height' ],
                             _source_data_delta ) +
                  meas.meas( pu_238_angled_data.loc[ 'lower_height' ],
                             _source_data_delta ) ) / 2
        
            
        xyz = meas.meas.from_list( [ x, y, z ] )
        xyz *= _MM_PER_INCH
        # xyz += np.array( [ 2, -2, 0 ] )
        all_coords.loc[source] = xyz

    # for source in all_objects :
    #     print( source + ': ' )
    #     print( all_coords.loc[ source ] )
    #     print( '' ) 
    #     # print( 'all_coords: ' + str(all_coords ) ) 


# these functions to be used as input for apply_nd
def costheta_from_3d_f( coords ):
    return coords[2] / np.sqrt( np.sum( coords ** 2 ) )  


def costheta_from_3d_fprime_tuple( coords ):
    r_sq = np.sum( coords ** 2 )
    return np.asarray ( [ -coords[0]*coords[2] / ( r_sq ** (3/2) ),
                          -coords[1]*coords[2] / ( r_sq ** (3/2) ),
                          ( r_sq - coords[2]**2 ) / ( r_sq ** (3/2) ) ] )

# input: measurement coords.
# output: cosine of angle between coords and z axis,
# along with uncertainty estimate
def costheta_from_3d( coords ):
    return coords.apply( costheta_from_3d_f,
                         costheta_from_3d_fprime_tuple )


    

# rotate a vector x about an axis ( 0 = x, 1 = y, 2 = z ) by theta
# applies to a np.ndarray
def rotate_3d( axis, theta, x, deg = 0 ):

    if deg :
        theta *= np.pi / 180
        
    # select the correct components to be rotated, using the fact that the
    # 3D rotation is identity for one component and 2d rotation of the other
    # two components 

    rot_indices = np.delete( np.arange(3), axis )

    # set the component that doesn't change.

    ret = np.empty(3)
    ret[axis] = x[axis]
    
    # do the 2D rotation and add it to ret
    
    ret[ rot_indices ] = rotate_2d( x[ rot_indices ], theta )

    return ret




# rotate meas x about an axis ( 0 = x, 1 = y, 2 = z ) by theta
# accounts for uncertainty calculation, hence the _meas suffix.
# we do not use the apply_nd function which is well-suited for
# maps from R^n to R. hopefully eventually there will be a similar
# function that can handle maps R^n to R^m
def rotate_3d_meas( axis, theta, x, deg=0 ):

    if deg :
        theta *= np.pi / 180 
    
    combined_meas = meas.append( theta, x )

    f = lambda _combined_meas : _rotate_3d_meas_f( axis, _combined_meas )
    fprime_tuple = lambda _combined_meas : _rotate_3d_meas_fprime_tuple( axis, _combined_meas )

    return combined_meas.apply_nd( f, fprime_tuple )




    
# wrapper for the rotate_3d function for an input of
# (theta, x ) combined into same measurement.
def _rotate_3d_meas_f( axis, combined_meas ):

    theta = combined_meas[0]
    x = combined_meas[1:]

    return rotate_3d( axis, theta, x ) 

        

# get 2d rotation matrix, can be used for different things.
def get_rotate_2d_matrix( theta ):
    return np.array( [ [ np.cos(theta), -np.sin(theta) ],
                       [ np.sin(theta),  np.cos(theta) ] ] )





def rotate_2d( x, theta ):
    return np.dot( get_rotate_2d_matrix( theta), x )



# these are the partial derivatives of the 2x2 rotation
# matrix
def _rotate_2d_fprime_tuple( theta, x ):

    a = x[0]
    b = x[1]
    
    return np.array( [ [ -a*np.sin(theta) - b*np.cos(theta),
                         a*np.cos(theta) - b*np.sin(theta) ],
                       [ np.cos(theta), np.sin(theta) ],
                       [ -np.sin(theta), np.cos(theta) ] ] )

    

def _rotate_3d_meas_fprime_tuple( axis, combined_meas ):

    theta = combined_meas[0]
    x = combined_meas[1:]

    # ret shall store 4 entries of 3-vectors. 0th entry is partial of R_{theta}(
    # x ) with respect to theta, rest are with respect to spatial
    # coords.
    ret = np.empty((4,3))

    rot_indices = np.array( list( range( 0, axis ) )
                            + list( range( axis+1, 3 ) ) )
    
    fprime_tuple_2d = _rotate_2d_fprime_tuple( theta, x[ rot_indices ] )

    # this is the partial derivative of the 3-vec resultant wrt theta.
    ret[0][rot_indices] = fprime_tuple_2d[0] 
    ret[0][axis] = 0

    # now handle the spatial coords.
    ret[ np.ix_( rot_indices + 1, rot_indices ) ] = fprime_tuple_2d[ 1: ]
    ret[ rot_indices + 1, axis ] = 0

    # evaluated at axis + 1 since 0th entry is reserved for theta derivative.
    ret[axis + 1] = unit_vec( axis, 3 )

    return ret





# generate a unit vec in R^n with entry 1 at k 
def unit_vec( k, n = 3 ):
    if k > n:
        raise ValueError( 'k must be less than n' )
    ret = np.zeros( n )
    ret[k] = 1
    return ret





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
                        








# def _populate_source_theta_phi( source_theta, source_phi, source_data, pu_238_angled_data ):

#     # assume all measured (or implicitly measured, e.g. assuming something is flat ) are known to
#     # 1 degree.
#     angle_data_delta = np.deg2rad( 1.0 / 360 ) 

#     radius = meas.meas( source_data.loc[ 'pu_238_angled', 'diameter' ],
#                         _source_data_delta ).mean()
    
#     height_diff = ( meas.meas( pu_238_angled_data[ 'upper_height' ],
#                                _source_data_delta ) +
#                     meas.meas( pu_238_angled_data[ 'lower_height' ],
#                                _source_data_delta ).mean() )

#     # take inverse tan of opposite over adjacent
#     theta = meas.arctan( height_diff / radius ) 
     
#     source_theta[ 'pu_238_angled' ] = theta
#     source_phi[ 'pu_238_angled' ] = meas.meas( 0, angle_data_delta )
    
    
#     # add in theta_phi values for everything except pu_238_angled.
#     upright_sources = list(sources)
#     upright_sources.remove( 'pu_238_angled' )

#     # loop through and add the same angles for the other sources.
#     for source in upright_sources:
#         source_theta[ source ] = meas.meas( 0, angle_data_delta )
#         source_phi[ source ] = meas.meas( 0, angle_data_delta )
                                            



        


                          
# in: 3-tuple coords of detector, coords of a particular source, and theta / phi 2-tuple of the
# source if it is rotated (only used for source deadlayer calculation )
# out: pd.Series containing 32x32 matrix of penetration angle for each
# detector pixel and another for the angle of penetration of the source layer.
# theta is measured from 0 to pi with 0 on the z axis, phi from 0 to 2pi with 0 on the x axis,
# same as the direction of 'right'

def _populate_costheta_grid( cosine_matrices, all_coords, source_data,
                             compute_source_costheta = False, source = None ):
       
    det_coords = all_coords.loc['detector']

    pu_238_angled_normal = meas.meas.from_list(
        np.array( [ - ( meas.meas( pu_238_angled_data[ 'upper_height' ],
                                   _source_data_delta )
                        - meas.meas( pu_238_angled_data[ 'lower_height' ],
                                     _source_data_delta ) ),
                    meas.meas( 0, _source_data_delta ),
                    meas.meas(
                        source_data.loc[ 'pu_238_angled', 'diameter' ],
                        _source_data_delta ).mean() ] ) )

    pu_238_angled_normal *= _MM_PER_INCH
    
    sourcenum = -1

    for source in sources:

        sourcenum += 1
        
        # print( str(sourcenum) + ' / ' + str(len(sources)) + '...' )

        # extract coords and angels 
        source_coords = all_coords.loc[ source ]
        

        # rename matrices in order to enhance readability 
        det_costheta_grid = cosine_matrices[ source ][ 0 ]
        source_costheta_grid = cosine_matrices[ source ][ 1 ]

        first_pixel_coords = det_coords - source_coords

        # print( source + ': ' + str( first_pixel_coords ) )
        
        # keep shifting by 2 mm to get next coord 
        for i in range(32):
            for j in range(32):

                # this works since all the pixels are separated by 2 mm.
                #pixel_displacement = 2.0 * rotate_3d( 1, 6, np.array([ -j, i, 0 ] ), deg=1 )
                pixel_displacement = 2.0 * np.array( [ -j, i, 0 ] ) 
                displacement = first_pixel_coords + pixel_displacement

                                                
                # print 'displacement: ' + str( displacement )
                
                
                # # get angle rel to detector pixel
                # print 'fprime tuple: ' + str( costheta_from_3d_fprime_tuple( displacement.x ) )
                # print 'displacement dx: '  + str( displacement.dx ) 
                costheta = displacement.apply_nd( costheta_from_3d_f,
                                                  costheta_from_3d_fprime_tuple )
                
                det_costheta_grid[i][j] = costheta

                if compute_source_costheta:
                    if source != 'pu_238_angled' :
                        source_costheta_grid[i,j] = costheta

                    else:
                        # rotated_displacement = rotate_3d_meas( 1, -theta, displacement ) ) 
                        # source_costheta_grid[i][j] = costheta_from_3d( rotated_displacement )
                        
                        # todo: proper error anaysis.
                        tmp = meas.dot( pu_238_angled_normal, displacement )
                        tmp /=  ( np.linalg.norm( pu_238_angled_normal.x ) *
                                  np.linalg.norm( displacement.x ) )
                        
                        source_costheta_grid[i][j] = tmp        
                     
                                    
                    
                    
                    
                        




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

def get_cosine_matrices( compute_source_costheta = 0 ):


    # if debug:
    #     return np.ones{ ( len(sources), 2, 32, 32 ) ) / 2.0
    
    print( 'INFO: constructing costheta grid...' )

    cosine_matrices_labels = [ sources, ['detector', 'source'] ]
    cosine_matrices = dict( zip( sources,
                                 meas.meas.empty( ( len(sources), 2, 32, 32 ) ) ) )

    
    # measurements 
    source_data = _get_source_data() 
    source_data = source_data.set_index( sources_index )
    _fill_redundant_source_data( source_data )

    
    # # this shall be populated with the theta_phi angles for all detectors, which are 0
    # # for all but the pu_238_angled
    # source_theta = pd.Series( sources_index )
    # source_phi = pd.Series( sources_index )
    # _populate_source_theta_phi( source_theta, source_phi, source_data, pu_238_angled_data )
    
        
    # declare dataframe to store all coordinates, which are pd.Series of two 3-tuples, one for
    # value and one for delta.
    all_coords = pd.Series( index = all_objects_index )
    _populate_all_coords( all_coords, source_data ) 


    # get each array of values / uncertainties and add to the grid.
    _populate_costheta_grid( cosine_matrices, all_coords, source_data,
                             compute_source_costheta )

    
    return cosine_matrices
    





def get_secant_matrices( compute_source_costheta = 0 ) :

    cosine_matrices = get_cosine_matrices( compute_source_costheta )

    secant_matrices = {}

    for key in cosine_matrices.keys() :
        secant_matrices[ key ] = abs( 1 / cosine_matrices[ key ]  ) 

    return secant_matrices 
