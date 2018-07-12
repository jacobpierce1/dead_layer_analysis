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


# CONFIG 

filter_above_channel_delta = 3.0

# filter_above_sectheta = 1.3 

peak_indices = [ [1], [0,1] ]
num_peaks_per_source = [ len( x ) for x in peak_indices ] 


db_names = [ 'full_bkgd_tot' ]

storage_path = '../../../storage/'


# strip_coords = 'all'

strip_coords = None





num_dbs = len( db_names )
num_sources = len( peak_indices )
num_peaks_per_source = [ len(x) for x in peak_indices ]
num_dets = 4





# np.set_printoptions(threshold=np.nan)
# print( det_sectheta[0][0] )
# sys.exit(0)




# save secants for debugging

# print( source_sectheta[0][0][0] ) 

# for sourcenum in range(2) :
#     for detnum in range(4) :
#         data = np.copy( source_sectheta[ sourcenum ][ detnum ] ) 
#         data[ np.isnan( data ) ] = 0

#         # print( os.path.exists( '../../../debugging_with_mary/exp2/' ) )
        
#         savepath = ( '../../../debugging_with_mary/exp2_secant_matrices/det_%d_source_%d.csv'
#                      % ( detnum + 1, sourcenum ) )

#         # print( savepath ) 
#         # print( data.dtype ) 
        
#         np.savetxt( savepath, data, delimiter=",", fmt = '%.5e' )

# sys.exit(0) 

    
actual_energies = [ np.array( [ 3182.690 ] ),
                    np.array( [ 5762.64, 5804.77 ] ) ]




density_si = 2.328 # g / cm^3



# density_si_dioxide = 2.65

# si_interp_stopping_powers *= density_si * 1000 * 100 / 1e9 
# si_interp_energies *= 1000   # keV / MeV





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






model_params = dl_estimator.deadlayer_model_params( disable_sources = 1,
                                                    vary_det_deadlayer = 0,
                                                    interp_det_stopping_power = 1,
                                                    interp_source_stopping_powers = 0,
                                                    fstrips_requested = np.arange(1,30),
                                                    bstrips = np.arange( 1, 30 ),
                                                    fix_source_deadlayers = None,
                                                    one_source_constant = 0,
                                                    det_thickness = 0,
                                                    vary_source_thickness = 0 )


dbs = [ spec.spectrum_db( db_names[i], storage_path ) for i in range( len( db_names ) ) ]

channels = [ dbs[i].read_peak_values() for i in range( len( db_names ) ) ]

print( 'channels dimensions:')
print( len( channels ) )
print( len( channels[0] ) )
print( len( channels[0][0] ) )
print( len( channels[0][0][0] ) )

print( 'det_sectheta dimensions:' )
print( len( det_sectheta ) )
print( len( det_sectheta[0] ) )
print( len( det_sectheta[0][0] ) )
# print( len( channels[0][0][0] ) )


for d in range( num_dbs ) :
    for detnum in range( num_dets ) : 
        for i in range( num_sources ) :

            mask = np.isnan( det_sectheta[ detnum ][i] )

            for j in range( num_peaks_per_source[i] ) :
                channels[ d ][ detnum ][i][j][ mask ] = meas.nan


        

if filter_above_channel_delta > 0 :
    
    for d in range( num_dbs ) :
        for det in range( num_dets ) :
            for i in range( num_sources ) :
                for j in range( len( num_peaks_per_source ) ) :
                    mask = ( channels[d][ det ][i][j].dx > filter_above_channel_delta )
                    channels[d][ det ][i][j][ mask ] = meas.nan



                    
if peak_indices is not None :
    
    for d in range( num_dbs ) :
        for det in range( num_dets ) :
            for i in range( num_sources ) :
                channels[ d ][ det ][i] = [ channels[d][ det ][i][j]
                                            for x in peak_indices[i] ]
                    
                

num_peaks_per_source = [ len( peak_indices[i] ) for i in range( num_sources ) ]

                    


                    
        




# INITIAL FIT PARAMETERS 

source_deadlayer_guesses = [ [ 6., 6.], [3.,3.], [15.,15.,15.,15.] ] 
source_stopping_power_interps = [ None ] * 3 
det_deadlayer_guess = 100.0
calibration_coefs_guess = [ 2.0, 0.0 ]


# print( len( channels ) )
# print( len( channels[0] ) ) 
# tmp = channels[ 0 ][1]
# print( len( tmp ) )
# print( len( tmp[0] ) )
# print( len( tmp[0][0] ) ) 


# strip_coords = 'all' 
# strip_coords = [0,2]

for detnum in  range(4) :

    det_sectheta_tmp = [ det_sectheta[ detnum ] ]

    source_sectheta_tmp = [ source_sectheta[ detnum ] ]

    channels_tmp = [ channels[ 0 ][ detnum ] ] 

    dl_estimator.estimate_deadlayers( model_params,
                                      channels_tmp,
                                      actual_energies,
                                      det_sectheta_tmp,
                                      source_sectheta_tmp,
                                      source_stopping_power_interps,
                                      source_deadlayer_guesses,
                                      det_stopping_power_interp, det_deadlayer_guess,
                                      calibration_coefs_guess,
                                      names = db_names,
                                      strip_coords = strip_coords,
                                      figpath = '../../../storage/current_peaks_vs_sectheta/exp2_aggregate/det_%d/'
                                      % detnum ) 





# dl_estimator.estimate_deadlayers( dbs, source_indices,
#                                   model_params,
#                                   cut_high_sectheta = 0,
#                                   annotate = 0,
#                                   # view_pixel = [ dbmgr.moved, 5 ],
#                                   subtitle = 'Det 1: Absolute Calibration of Entire Detector\n',
#                                   reset_angles = None,
#                                   residual_scatter_plot = 0,
#                                   plot_3d = 1,
#                                   savefig_dir = '../../../deadlayer_analysis_paper/images/det1_calibration.eps' )













