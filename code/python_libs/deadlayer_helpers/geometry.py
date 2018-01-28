# the purpose of this module is to construct several global variables
# that are extremely useful, such as grids of cosine theta values for
# all the detector pixels.  the goal is that geometry issues and
# calculations will never need to be considered outside the module.


# my includes 
# import deadlayer_helpers.sql_db_manager as db
import libjacob.jmeas as meas



## includes 
from mpl_toolkits.mplot3d import axes3d
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import scipy.integrate 

import os




_current_abs_path = os.path.dirname( __file__ ) + '/'





_CM_PER_INCH = 2.54
_MM_PER_INCH = 25.4


USE_MARY_DISPLACEMENTS = 0

# these are the first pixel displacements as computed by Mary
mary_first_pixel_displacements = { 'pu_240' : meas.meas( [87.10, -9.31, 58.35], np.zeros(3) ),
                                   'cf_249' : meas.meas( [85.35, -46.06, 57.74], np.zeros(3) ),
                                   'pu_238_centered': meas.meas( [32.95, -31.45, 57.88], np.zeros(3) ),
                                   'pu_238_flat': meas.meas( [32.95, -31.45, 57.88], np.zeros(3) ),
                                   'pu_238_angled': meas.meas( [32.95, -31.45, 58.72], np.zeros(3) ),
                                   'pu_238_moved' : meas.meas( [-44.59, -28.58, 57.88], np.zeros(3) ) } 


# mary_first_pixel_displacements = { 'pu_240' : meas.meas( [87.10, -9.31, 58.35], np.zeros(3) ),
#                                    'cf_249' : meas.meas( [85.35, -46.06, 57.74], np.zeros(3) ),
#                                    'pu_238_centered': meas.meas( [32.95, -31.45, 57.88], np.zeros(3) ),
#                                    'pu_238_flat': meas.meas( [32.95, -31.45, 57.88], np.zeros(3) ),
#                                    'pu_238_angled': meas.meas( [32.95, -31.45, 58.72], np.zeros(3) ),
#                                    'pu_238_moved' : meas.meas( [-44.59, -28.58, 57.88], np.zeros(3) ) } 


                                   
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
            
            'source_diameter' : [ 2, 8, 2, 2, 2, 2 ]
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


    x_measurement_inverted = 0
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

        
    # if y_measurement_inverted :
    #     shift = - detector_data[ 'total_width' ]
    #     det_coords += np.array( [0, shift, 0 ] )

    # else:
    #     shift = ( ceramic_data[ 'total_y' ]
    #               - detector_data[ 'total_width' ] )
    #     det_coords += np.array( [0, shift, 0 ] )


    
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
                        





def average_sectheta_integrand( R, nhat, r, theta ) :

    rprime = np.array( [ r * np.cos(theta),
                         r * np.sin(theta),
                         0 ] )
    
    x = R - rprime
    
    return np.abs( np.dot( nhat, x ) ) * r / np.dot( x, x ) 




def source_sectheta_distribution( R, nhat, r, theta ) :

    rprime = np.array( [ r * np.cos(theta),
                         r * np.sin(theta),
                         0 ] )
    
    x = R - rprime

    return np.abs( np.dot( nhat, x ) ) * r / np.dot( x, x ) ** 1.5
    





def compute_average_sectheta_over_source( R, nhat, radius ) :
                        
    # todo: do rotation here if necessary. requires normal vector of the
    # source.

    # compute normalizing factor for the distribution.

    A_integrand = lambda r, theta : source_sectheta_distribution( R, nhat, r, theta ) 
    
    A = scipy.integrate.dblquad( A_integrand,
                                 0, 2 * np.pi,
                                 lambda x : 0,
                                 lambda x : radius )

    ave_sectheta_integrand = lambda r, theta : average_sectheta_integrand( R, nhat, r, theta )

    ave_sectheta = scipy.integrate.dblquad( ave_sectheta_integrand,
                                            0, 2 * np.pi,
                                            lambda x : 0,
                                            lambda x : radius )
    
    return ave_sectheta[0] / ( A[0] * R[2] ) 



        



                          
# in: 3-tuple coords of detector, coords of a particular source, and theta / phi 2-tuple of the
# source if it is rotated (only used for source deadlayer calculation )
# out: pd.Series containing 32x32 matrix of penetration angle for each
# detector pixel and another for the angle of penetration of the source layer.
# theta is measured from 0 to pi with 0 on the z axis, phi from 0 to 2pi with 0 on the x axis,
# same as the direction of 'right'

def _populate_sectheta_grid( secant_matrices, all_coords, source_data,
                             compute_source_costheta = False,
                             average_over_source = False ):
       
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

    # print( pu_238_angled_normal )
    
    pu_238_angled_normal *= _MM_PER_INCH

    # print( pu_238_angled_normal ) 
    
    sourcenum = -1

    for source in sources :

        sourcenum += 1
        
        # print( str(sourcenum) + ' / ' + str(len(sources)) + '...' )

        # extract coords and angels 
        source_coords = all_coords.loc[ source ]
        

        # rename matrices in order to enhance readability 
        det_sectheta_grid = secant_matrices[ source ][ 0 ]
        source_sectheta_grid = secant_matrices[ source ][ 1 ]

        first_pixel_coords = det_coords - source_coords

        # debug: replace the first pixel coords with mary's
        if USE_MARY_DISPLACEMENTS : 
            first_pixel_coords = mary_first_pixel_displacements[ source ] 

        print( source + ': ' + str( first_pixel_coords ) )
        print( 'mary: ' + str( mary_first_pixel_displacements[ source ]  ) )
        
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
                
                sectheta = 1 / costheta

                det_sectheta_grid[i][j] = sectheta

                if compute_source_costheta:

                    if not average_over_source : 
                        
                        if source != 'pu_238_angled' :
                            source_sectheta_grid[i,j] = sectheta
                        
                        else:
                            # rotated_displacement = rotate_3d_meas( 1, -theta, displacement ) ) 
                            # source_costheta_grid[i][j] = costheta_from_3d( rotated_displacement )
                        
                            # todo: proper error anaysis.
                            tmp = ( np.linalg.norm( pu_238_angled_normal.x ) *
                                    np.linalg.norm( displacement.x ) )
                            
                            tmp /= meas.dot( pu_238_angled_normal, displacement )
                            
                            source_sectheta_grid[i][j] = tmp        

                    else:

                        # todo: haven't figured out how to do this calculation yet.
                        if source == 'pu_238_angled' :
                            tmp = ( np.linalg.norm( pu_238_angled_normal.x ) *
                                    np.linalg.norm( displacement.x ) )
                            
                            tmp /= meas.dot( pu_238_angled_normal, displacement )
                            
                            source_sectheta_grid[i][j] = tmp
                            
                        else:
                            source_radius = source_data.loc[ source, 'source_diameter' ] / 2 
                            ave_sectheta = (
                                compute_average_sectheta_over_source( displacement.x,
                                                                      np.array( [0.0, 0.0, 1.0] ),
                                                                      source_radius ) )

                            source_sectheta_grid[i][j] = meas.meas( ave_sectheta, 0 ) 
                    
                    

                            

                            
                            

# construct matrices containing secant of penetration
# angle through all sources as well as the detector.

def get_secant_matrices( compute_source_sectheta = 0,
                         average_over_source = 0,
                         reset = 0 ):

    data_path =  _current_abs_path + '../../../storage/secant_matrices/'

    if average_over_source :
        data_path += 'average_over_source/'
    else:
        data_path += 'regular/'
        
    if not os.path.exists( data_path ) :
        os.makedirs( data_path )


    # what to do if not rewriting all the files: read them from
    # existing location 
        
    if not reset:
        print( 'INFO: attempting to read secant matrices from disk...' )

        secant_matrices = {}

        if all( [ os.path.exists( data_path + key + z + '.bin' )
                  for z in [ '_x', '_dx' ]
                  for key in sources ] ) :

            for key in sources :
                
                xpath, dxpath = [ data_path + key + z + '.bin'
                                  for z in [ '_x', '_dx' ] ]
                if os.path.exists( xpath ) and os.path.exists( dxpath ) :

                    secant_matrices[ key ] = meas.meas( np.fromfile( xpath ).reshape( 2, 32, 32 ),
                                                        np.fromfile( dxpath ).reshape( 2, 32, 32 ) )

            print( 'INFO: success.' )
            return secant_matrices

        else:
            print( 'INFO: not all matrices present, reconstructing...' )
    
    else:
        print( 'INFO: reconstructing sectheta grid...' )


        
    secant_matrices_labels = [ sources, ['detector', 'source'] ]
    secant_matrices = dict( zip( sources,
                                 meas.meas.empty( ( len(sources), 2, 32, 32 ) ) ) )

    
    # measurements 
    source_data = _get_source_data() 
    source_data = source_data.set_index( sources_index )
    _fill_redundant_source_data( source_data )

            
    # declare dataframe to store all coordinates, which are pd.Series of two 3-tuples, one for
    # value and one for delta.
    all_coords = pd.Series( index = all_objects_index )
    _populate_all_coords( all_coords, source_data ) 


    # get each array of values / uncertainties and add to the grid.
    _populate_sectheta_grid( secant_matrices, all_coords, source_data,
                             compute_source_sectheta,
                             average_over_source )

    for key, val in secant_matrices.items() :
        secant_matrices[key] = abs( val )


    # write the arrays to files, both bin and csv.
    
    for key in sources :

        xpath, dxpath = [ [ data_path + key + z + suffix
                            for suffix in [ '.bin', '.csv' ] ]
                          for z in [ '_x', '_dx' ] ]

        secant_matrices[ key ].x.tofile( xpath[0] ) 
        secant_matrices[ key ].dx.tofile( dxpath[0] )
        
        
        np.savetxt( xpath[1], secant_matrices[ key ].x[0], delimiter = ',', fmt = '%4f' )

        if key == 'pu_238_angled' :
            np.savetxt( xpath[1].replace( key, key + '_source_sectheta' ),
                        secant_matrices[ key ].x[1],
                        delimiter = ',', fmt = '%4f' )

        
        
    return secant_matrices
    






# def get_secant_matrices( compute_source_costheta = 0,
#                          average_over_source = 0 ) :

#     cosine_matrices = get_cosine_matrices( compute_source_costheta,
#                                            average_over_source )

#     secant_matrices = {}

#     for key in cosine_matrices.keys() :
#         secant_matrices[ key ] = abs( 1 / cosine_matrices[ key ]  ) 

#     return secant_matrices 
