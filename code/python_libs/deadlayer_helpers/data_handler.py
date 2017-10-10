# the purpose of this python module is to provide functions that
# access directory-structure dependent file paths. if the directory is
# rearranged, the hope is that the only changes would need to occur in
# this file instead of modifying many others. note that this was not
# the case when this module was initially made, which is why it is not
# used in some of the original scripts.



import array
import numpy as np
from scipy import special



from sql_db_manager import centered_db, moved_db, flat_db, angled_db
import deadlayer_helpers.functions as dlf


_current_abs_path = os.path.dirname( __file__ ) + '/'


# construct the histogram of a particular
# pixel in the detector. id is 'centered', 'moved',
# 'flat', or 'angled'.
def get_pixel_histo( id, x, y, histo_in ):

    datafile = ( _current_abs_path + '../../../data/extracted_ttree_data/' +
id + '/' + id + '_' + x + '_' + y + '.bin' )

    dlf.construct_histo_array( datafile, histo_in )

    
