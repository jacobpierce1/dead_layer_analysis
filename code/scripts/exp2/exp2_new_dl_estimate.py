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
import exp2_geometry


filter_above_sectheta = 0

data_name = 'full_bkgd_tot' 

cut_data = 1

storage_path = '../../../storage/'


# strip_coords = 'all'

strip_coords = None

num_dets = 4
num_sources = 2



# def remove_strips( dssd, fstrips, bstrips ) :
#     for i in range( dssd.num_groups ) :
#         for j in range( dssd.num_data_per_group[i] ) : 
#             for d in range( dssd.num_groups ) : 
#             for f in fstrips :
#                 dssd[ i,j,0,f, : ] = meas.nan
#             for b in bstrips :
#                 dssd[i,j,0,:, b ] = meas.nan



        
# def remove_data( dssd, data_to_remove ) :
#     print( 'removing data: ', data_to_remove ) 
#     for i in range( dssd.num_groups ) :
#         for j in range( dssd.num_data_per_group[i] ) : 
#             for coords in data_to_remove :
#         # print( coords ) 
#                 dssd[ i,j,0, coords ] = meas.nan




analysis_mgr = spec.dssd_analysis_manager( data_name, storage_path, (1,1), [1,1,1] )


actual_energies = [  [ 3182.690 ] ,
                     [ 0.231 * 5762.64  + 0.769 * 5804.77 ] ]

num_peaks_per_source = [ len( actual_energies[i] ) for i in range( num_sources ) ]

                               


channels = analysis_mgr.load_dill( 'means' )
channels = spec.dssd_data_container( [ channels[1], channels[2] ] )
# del channels[0]

peak_indices = [ [0], [0] ]

print( 'actual_energies', actual_energies ) 



num_peaks_per_source = [ len( actual_energies[i] ) for i in range( num_sources ) ]


secant_matrices = exp2_geometry.get_secant_matrices() 
secant_matrices = np.delete( secant_matrices, 0, axis = 0 )


# if filter_above_sectheta > 0 :
    
#     for det in range( num_dets ) :
#         for i in range( num_sources ) :
#             for j in range( num_peaks_per_source[i] ) :
#                 mask = ( secant_matrices[i][ det ] > filter_above_sectheta )
#                 channels[i][j][ det ][ mask ] = meas.nan



if peak_indices is not None :
    for i in range( num_sources ) :
        channels[i] = [ channels[i][j] for j in peak_indices[i] ]
                    


for i in range(2) :
    for j in range(2) :
        mask = ~ np.isnan( secant_matrices[i,j] )
        secant_matrices[i,0][ mask ] = secant_matrices[i,j][mask] 

secant_matrices = [ secant_matrices[i,0] for i in range(2) ]

print( type( secant_matrices ) )
print( type( secant_matrices[0] ) ) 


num_peaks_per_source = [ len( peak_indices[i] ) for i in range( num_sources ) ]    
print( 'num_peaks_per_source', num_peaks_per_source )

for i in range( num_sources ) :
    mask = np.isnan( secant_matrices[i] )
    print( mask.shape ) 
    for j in range( num_peaks_per_source[i] ) :
        channels[i][j][ mask ] = meas.nan

# print( secant_matrices[0][0] )
# print( channels[0] ) 

print( len( secant_matrices ) )
print( len( secant_matrices[0] ) )
print( len( secant_matrices[0][0] ) )



strip_cuts = [ [0,31], range(16), [31] ] 


# print( secant_matrices ) 

if cut_data :
    channels.remove_strips( strip_cuts[0], strip_cuts[1] ) 
    # remove_data( channels, data_cuts ) 

    # remove_strips( channels, [0] + list( range( 2,31) ), [] )
    
# if filter_above_channel_delta > 0 :
    
#     for det in range( num_dets ) :
#         for i in range( num_sources ) :
#             for j in range( len( num_peaks_per_source ) ) :
#                 mask = ( channels[i][j][ det ].dx > filter_above_channel_delta )
#                 channels[i][j][ det ][ mask ] = meas.nan







savename = 'dl_regression_indep'

source_names = [ '$^{148}$Gd', '$^{264}$Cm' ]

dl_estimator.independent_dl_regression( channels, secant_matrices, actual_energies,
                                        show = 1, savename = savename,
                                        analysis_mgr = analysis_mgr,
                                        source_names = source_names,
                                        # dets = [0], fstrips = [25],
                                        # title = 'Det 1, Strip 25: Source Distribution Means',
                                        dpi = 500 ) 
