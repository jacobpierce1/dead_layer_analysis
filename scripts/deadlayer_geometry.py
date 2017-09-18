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
def get_costheta( x, y, z ):
    return z / np.sqrt( x**2.0 + y**2.0 + z**2.0 )
    
    
def get_costheta_delta( x,dx, y,dy, z,dz ) :
    return ( ( z**2 * ( (x*dx)**2 + (y*dy)**2 )  +
                    ( x**2 + y**2 )**2 * dz**2 )  /
                        (x**2 + y**2 + z**2 )**3 )




# get the coords detector referenced to 0
det_coords = pd.Series(
    [
        ( enclosure_data['detector_xy'][0],
          enclosure_data['detector_xy'][1],
          enclosure_data['total_height'],
          - enclosure_data['top_offset'],
          - enclosure_data['detector_height']
        ),
        
        ( 0,0,0 ),
    ],
    index = values_index
)


all_coords = pd.DataFrame(
    columns = values_index,
    index = all_objects_index 
)



def populate_source_coords( source_coords ):

    for source in sources_index:

        # xyz is 3 tuple storing coordinates, xyz_delta stores uncertainties.
        xyz = (0,0,0)
        xyz_delta = (0,0,0)
        
        # first get x and y from the diameters / distance to edges.
        configurations = [ ['left','right'], ['top','bottom'] ]
        for i in range(len(configurations)):

            # these coords are the sum of the diameter plus distance to walls over 2 
            xyz[i] = np.sum( [ source_data[col].loc[source] for col in configurations[i] ] )
            xyz[i] += source_data['diameter'].loc[source]
            xyz[i] /= 2

            # for uncertainty
            xyz_delta[i] = 

        # add the coordinates and uncertainties to the source_data
        all_coords[values_index[0]].loc[source] = xyz
        all_coords[values_index[1]].loc[source] = xyz_delta
