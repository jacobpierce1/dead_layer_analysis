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
import exp1_geometry


filter_above_sectheta = 0


data_name = 'det1_moved' 

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





# np.set_printoptions(threshold=np.nan)
# print( det_sectheta[0][0] )
# sys.exit(0)




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






# actual_energies = [ 5168.17 ]

num_peaks_per_source = [ len( actual_energies[i] ) for i in range( num_sources ) ]


density_si = 2.328 # g / cm^3




# interpolate stopping power of alphas in Si using data
# in the range emin, emax. data is from astar program of NIST.
# uses scipy.interpolate.interp1d. if from_0_keV is 1, then
# interpolate over a lot of data to obtain large interpolation
# and set the 

def construct_si_stopping_power_interpolation( plot = 0 ) :

    # data = np.loadtxt( '../../data/stopping_power_data/alpha_stopping_power_si.txt',
    #                   skiprows = 10, unpack = 1 )

    energy, stopping_power = np.loadtxt(
        '../../../data/stopping_power_data/alpha_stopping_power_si.txt',
        usecols = [0,3], unpack = 1 )

    # tmp = np.loadtxt(
    #     '../../../data/stopping_power_data/alpha_stopping_power_si.txt',
    #     usecols = [0,3], unpack = 1 )

    
    # energy = data[0] * 1000
    
    # energy = energy[ energy <= emax ] 
    
    # stopping_power = data[3][ 0 : len( energy ) ]
    energy *= 1000 
    stopping_power *= density_si * 1000 * 100 / 1e9

    print( 'stopping power interp data:' )
    print( 'energy: ', energy )
    print( 'stopping: ', stopping_power )

    # add particular data points of interest to the interpolation
    
    interp = scipy.interpolate.interp1d( energy, stopping_power, kind = 'cubic' )

    
    if plot :

        ax = plt.axes()

        interp_axis = np.linspace( min( energy ), max( energy ), 100 )
        
        ax.scatter( energy, stopping_power, color='r' )

        ax.plot( interp_axis,
                 interp( interp_axis ),
                 color = 'b' )
        
        plt.show()

        return 1
        

    return interp
    


det_stopping_power_interp = construct_si_stopping_power_interpolation()

# print( det_stopping_power_interp( 5.5 ) )
# print( det_stopping_power_interp( 3.2 ) )
# sys.exit(0) 






model_params = dl_estimator.deadlayer_model_params( disable_sources = 0,
                                                    vary_det_deadlayer = 1,
                                                    interp_det_stopping_power = 1,
                                                    interp_source_stopping_powers = 0,
                                                    fstrips_requested = np.arange(1,30),
                                                    bstrips = np.arange( 25, 30 ),
                                                    fix_source_deadlayers = None,
                                                    one_source_constant = 0,
                                                    det_thickness = 0,
                                                    vary_source_thickness = 0,
                                                    constant_energy_offset = 0 )

# sys.exit( 0 ) 







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

                    
# print( 'det_sectheta dimensions:' )
# print( len( det_sectheta ) )
# print( len( det_sectheta[0] ) )
# print( len( det_sectheta[0][0] ) )
# # print( len( channels[0][0][0] ) )


                    
if peak_indices is not None :
    for i in range( num_sources ) :
        channels[i] = [ channels[i][j]
                        for j in peak_indices[i] ]
                    
                

num_peaks_per_source = [ len( peak_indices[i] ) for i in range( num_sources ) ]    
print( 'num_peaks_per_source', num_peaks_per_source )

for d in range( num_dets ) : 
    for i in range( num_sources ) :
        
        mask = np.isnan( det_sectheta[i][d] )

        for j in range( num_peaks_per_source[i] ) :
            channels[i][j][ d ][ mask ] = meas.nan
    

strip_cuts = [ [0,31], [] ] # list( range(16) ) + [ 31]  ] 
data_cuts = [ [18,26] ] 

# print( channels.values )
print( det_sectheta ) 

if cut_data :
    remove_strips( channels, strip_cuts[0], strip_cuts[1] ) 
    remove_data( channels, data_cuts ) 

    # remove_strips( channels, [0] + list( range( 2,31) ), [] )
    
# if filter_above_channel_delta > 0 :
    
#     for det in range( num_dets ) :
#         for i in range( num_sources ) :
#             for j in range( len( num_peaks_per_source ) ) :
#                 mask = ( channels[i][j][ det ].dx > filter_above_channel_delta )
#                 channels[i][j][ det ][ mask ] = meas.nan





# INITIAL FIT PARAMETERS 

# source_deadlayer_guesses = [ [ 6., 6.], [3.,3.], [15.,15.,15.,15.] ] 
source_stopping_power_interps = [ None ] * 3 
det_deadlayer_guess = 100.0
calibration_coefs_guess = [ 2.0, 0.0 ]
source_deadlayer_guesses = [ [25.0], [25.0] ] 

# print( channels[0,0,0,1] )

for detnum in  range(1) :


    det_sectheta_tmp = [ [ det_sectheta[i][ detnum ]
                           for i in range( num_sources ) ] ]

    source_sectheta_tmp = det_sectheta_tmp 

    channels_tmp = [ [ [ channels[i][j][ detnum ]
                         for j in range( num_peaks_per_source[i] ) ] 
                       for i in range( num_sources ) ] ]
    
    dl_estimator.estimate_deadlayers( model_params,
                                      channels_tmp,
                                      actual_energies,
                                      det_sectheta_tmp,
                                      source_sectheta_tmp,
                                      source_stopping_power_interps,
                                      source_deadlayer_guesses,
                                      det_stopping_power_interp, det_deadlayer_guess,
                                      calibration_coefs_guess,
                                      names = data_name,
                                      strip_coords = strip_coords ) 
                                      # figpath = '../../../storage/current_peaks_vs_sectheta/exp2_aggregate/det_%d/'
                                      # % detnum ) 














