import matplotlib
matplotlib.use('agg') 

import sys
import numpy as np

import os 

# from .. import deadlayer_estimator as dl_estimator 
sys.path.append('../')
import deadlayer_estimator as dl_estimator

# import exp1_geometry as geom 

import jspectroscopy as spec
import jutils 
import jutils.meas as meas 

import scipy.interpolate
import exp2_geometry

# CONFIG 

cut_data = 0


filter_above_sectheta = None

data_names = [ 'full_bkgd_tot' ]

storage_path = '../../../storage/'


# strip_coords = 'all'

strip_coords = None



num_dets = 4 
num_sources = 2



    
# actual_energies = [ np.array( [ 3182.690 ] ),
#                     np.array( [ 5762.64, 5804.77 ] ) ]

    
actual_energies = [  [ 3182.690 ] ,
                     [ 0.231 * 5762.64  + 0.769 * 5804.77 ] ]

num_peaks_per_source = [ len( actual_energies[i] ) for i in range( num_sources ) ]


det_stopping_power_interp = None



model_params = dl_estimator.deadlayer_model_params( disable_sources = 0,
                                                    vary_det_deadlayer = 1,
                                                    interp_det_stopping_power = 1,
                                                    interp_source_stopping_powers = 0,
                                                    fstrips_requested = np.arange( 1, 31 ),
                                                    bstrips = np.arange( 15, 28 ),
                                                    fix_source_deadlayers = None,
                                                    one_source_constant = 0,
                                                    det_thickness = 0,
                                                    vary_source_thickness = 0,
                                                    constant_energy_offset = 0 )


analysis_mgr = spec.dssd_analysis_manager( 'full_bkgd_tot', storage_path, (32,32), [1,1,1] )

# channels = analysis_mgr.load_dill( 'peaks' )
# peak_indices = [ [1], [0,1] ]
# num_peaks_per_source = [ len( x ) for x in peak_indices ]


channels = analysis_mgr.load_dill( 'means' )[1:]


# peak_indices = [ [0], [0] ]
peak_indices = None

det_sectheta = exp2_geometry.get_secant_matrices()[1:]

source_sectheta = det_sectheta



                    
if peak_indices is not None :
    for i in range( num_sources ) :
        channels[i] = [ channels[i][j]
                        for j in peak_indices[i] ]

    num_peaks_per_source = [ len( peak_indices[i] ) for i in range( num_sources ) ]    
    print( 'num_peaks_per_source', num_peaks_per_source )


    

if filter_above_sectheta is not None :
    
    for det in range( num_dets ) :
        for i in range( num_sources ) :
            for j in range( num_peaks_per_source[i] ) :
                mask = ( det_sectheta[i][ det ] > filter_above_sectheta )
                channels[i][j][ det ][ mask ] = meas.nan

                    
        


def remove_strips_from_dets( dssd, fstrips, bstrips ) :
    for i in range( len(dssd) ) :
        for j in range( len(dssd[i] ) ) : 
            for d in range( len(dssd[i][j] ) ) : 
                for f in fstrips[d] :
                    dssd[ i][j][d,f, : ] = meas.nan
                for b in bstrips :
                    dssd[i][j][d,:, b ] = meas.nan



        
def remove_data_from_dets( dssd, data_to_remove ) :
    for i in range( len(dssd) ) :
        for j in range( len(dssd[i] ) ) : 
            for d in range( len(dssd[i][j] ) ) : 
                for coords in data_to_remove[d] :
                    dssd[ i][j][d, coords[0], coords[1] ] = meas.nan

                
fstrip_cuts = [ [ 28, 29 ],
                [],
                [],
                [] ]

data_cuts =  [ [], [], [], [] ]


if cut_data :
    remove_strips_from_dets( channels, fstrip_cuts, [] )
    remove_data_from_dets( channels, data_cuts )


# INITIAL FIT PARAMETERS 

# source_deadlayer_guesses = [ [ 6., 6.], [3.,3.], [15.,15.,15.,15.] ] 
source_stopping_power_interps = [ None ] * 3 
det_deadlayer_guess = 100.0
calibration_coefs_guess = [ 2.0, 0.0 ]
source_deadlayer_guesses = [ [25.0], [25.0] ] 



ret = meas.empty( ( 2, 2, 4 ) ) 


for d in range(4) :
    for k in range(2) : 

        # det_sectheta_tmp = [ det_sectheta[ d ] ]
        # source_sectheta_tmp = [ source_sectheta[ d ] ]
        # channels_tmp = [ channels[ 0 ][ d ] ]

        det_sectheta_tmp = [ [ det_sectheta[i][k][ d ]
                               for i in range( num_sources ) ] ]
        
        source_sectheta_tmp = det_sectheta_tmp 

        channels_tmp = [ [ [ channels[i][j][ d ]
                             for j in range( num_peaks_per_source[i] ) ] 
                           for i in range( num_sources ) ] ]

        savepath = ( analysis_mgr.storage_path + '/dl_regression/%d/%d/'
                     % ( analysis_mgr.detidx_to_detnum(d), k )   )

        
        dl_result = dl_estimator.estimate_deadlayers( model_params,
                                                      channels_tmp,
                                                      actual_energies,
                                                      det_sectheta_tmp,
                                                      source_sectheta_tmp,
                                                      source_stopping_power_interps,
                                                      source_deadlayer_guesses,
                                                      det_stopping_power_interp, det_deadlayer_guess,
                                                      calibration_coefs_guess,
                                                      names = data_names,
                                                      strip_coords = strip_coords,
                                                      savepath = savepath  )
        tmp1 = dl_result.params[ 'source_constant_0_0' ]
        tmp2 = dl_result.params[ 'source_constant_1_0' ]
        
        ret[ :,k,d ] = meas.meas( [ tmp1.value, tmp2.value ], [ tmp1.stderr, tmp2.stderr ] )


analysis_mgr.save_dill( ret, 'agg_dl_sums' )
