import matplotlib
matplotlib.use( 'agg' ) 


import sys
import numpy as np

import os 

# from .. import deadlayer_estimator as dl_estimator 
sys.path.append('../')
import new_deadlayer_estimator as dl_estimator

# import exp1_geometry as geom 

import jspectroscopy as spec
import jutils 
import jutils.meas as meas 

import scipy.interpolate
import exp1_geometry


filter_above_sectheta = 0

data_name = 'det3_moved' 

cut_data = 1

storage_path = '../../../storage/'


# strip_coords = 'all'

strip_coords = None

num_dets = 1
num_sources = 2



def remove_strips( dssd, fstrips, bstrips ) :
    for i in range( dssd.num_groups ) :
        for j in range( dssd.num_data_per_group[i] ) : 
            for f in fstrips :
                dssd[ i,j,0,f, : ] = meas.nan
            for b in bstrips :
                dssd[i,j,0,:, b ] = meas.nan



        
def remove_data( dssd, data_to_remove ) :
    print( 'removing data: ', data_to_remove ) 
    for i in range( dssd.num_groups ) :
        for j in range( dssd.num_data_per_group[i] ) : 
            for coords in data_to_remove :
        # print( coords ) 
                dssd[ i,j,0, coords ] = meas.nan




analysis_mgr = spec.dssd_analysis_manager( data_name, storage_path, (1,1), [1,1,1] )

    


cf249_energies = np.array( [  5704, 5759.5, 5783.3, 5813.3, 5849.3 ] )
cf249_probs = np.array( [ 0.0003, 0.0469, 0.0026, 0.822, 0.01430  ] )
cf249_probs /= np.sum( cf249_probs )

pu240_energies = np.array( [ 5123.68, 5168.17 ] ) 
pu240_probs = np.array( [ 0.2710, 0.7280 ] )
pu240_probs /= np.sum( pu240_probs ) 

actual_energies = [  [ np.dot( pu240_energies, pu240_probs ) ],
                     [ np.dot( cf249_energies, cf249_probs ) ] ]
                               


channels = analysis_mgr.load_dill( 'means' )
channels = spec.dssd_data_container( [ channels[0], channels[2] ] )

peak_indices = [ [0], [0] ]

print( 'actual_energies', actual_energies ) 


# channels = analysis_mgr.load_dill( 'peaks' )
# del channels[1]
# actual_energies = [ np.array( [ 3182.690 ] ),
#                     np.array( [ 5813.3 ] ) ]

# peak_indices = [ [0], [1] ]

# approx_slope = ( (actual_energies[0][0] - actual_energies[1][0] )
#                  / ( channels[ 0, 0, 0, 16, 16 ] - channels[ 1, 0, 0, 16, 16 ] ) )

# print( 'approx_slope', approx_slope  )



num_peaks_per_source = [ len( actual_energies[i] ) for i in range( num_sources ) ]


det_sectheta = exp1_geometry.get_secant_matrices()[ data_name ] 
det_sectheta = np.delete( det_sectheta, 1, axis = 0 )

# print( 'after cut' ) 
# print( det_sectheta )
# sys.exit(0)


if filter_above_sectheta > 0 :
    
    for det in range( num_dets ) :
        for i in range( num_sources ) :
            for j in range( num_peaks_per_source[i] ) :
                mask = ( det_sectheta[i][ det ] > filter_above_sectheta )
                channels[i][j][ det ][ mask ] = meas.nan

                    
if peak_indices is not None :
    for i in range( num_sources ) :
        channels[i] = [ channels[i][j] for j in peak_indices[i] ]
                    
                

num_peaks_per_source = [ len( peak_indices[i] ) for i in range( num_sources ) ]    
print( 'num_peaks_per_source', num_peaks_per_source )

for d in range( num_dets ) : 
    for i in range( num_sources ) :
        
        mask = np.isnan( det_sectheta[i][d] )

        for j in range( num_peaks_per_source[i] ) :
            channels[i][j][ d ][ mask ] = meas.nan
    

strip_cuts = [ [0,31], range(21) ] # list( range(16) ) + [ 31]  ] 
# data_cuts = [ [18,26] ] 

# print( channels.values )
print( det_sectheta ) 

if cut_data :
    remove_strips( channels, strip_cuts[0], strip_cuts[1] ) 
    # remove_data( channels, data_cuts ) 

    # remove_strips( channels, [0] + list( range( 2,31) ), [] )
    
# if filter_above_channel_delta > 0 :
    
#     for det in range( num_dets ) :
#         for i in range( num_sources ) :
#             for j in range( len( num_peaks_per_source ) ) :
#                 mask = ( channels[i][j][ det ].dx > filter_above_channel_delta )
#                 channels[i][j][ det ][ mask ] = meas.nan







savepath = analysis_mgr.storage_path + '/dl_regression/'

dl_estimator.independent_dl_regression( channels, det_sectheta, actual_energies,
                                        show = 1, save = 1, analysis_mgr = analysis_mgr  ) 
