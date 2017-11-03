# the purpose of this python module is to provide functions that
# access directory-structure dependent file paths. if the directory is
# rearranged, the hope is that the only changes would need to occur in
# this file instead of modifying many others. note that this was not
# the case when this module was initially made, which is why it is not
# used in some of the original scripts.


import os 
import array
import numpy as np
from scipy import special

# import deadlayer_helpers.sql_db_manager as db



_current_abs_path = os.path.dirname( __file__ ) + '/'



# construct the histogram of a particular
# pixel in the detector. id is 'centered', 'moved',
# 'flat', or 'angled'. these are defined in
# deadlayer_helpers.sql_db_manager.

def get_pixel_histo( db, x, y, histo_in ):
    
    datafile = ( _current_abs_path + '../../../data/extracted_ttree_data/' +
                 db.name + '/' + db.name + '_' + str( x ) + '_' + str( y ) + '.bin' )

    construct_histo_array( datafile, histo_in )

    


    
# extract next row from the file. file assumed to already be opened.
def read_doubles_from_bin( f, buf ):
    bytes = f.read(8*2 )  # since we write 2 doubles for each event: efront and eback.
    if( not bytes ): return 0;
    buf[:] = array.array('d', bytes) 
    return 1




# open binary infile and populate the array efront_histo as a histogram 
def construct_histo_array( infile, efront_histo ):

    try:
        with open( infile, "rb" ) as f:
            buf = array.array( 'd', [0.0, 0.0] )  # data gets read into here.
            while( read_doubles_from_bin( f, buf ) ):
                efront_histo[ int( buf[0] ) ] += 1
                
    except IOError:
        return 0

    return 1
