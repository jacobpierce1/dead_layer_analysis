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

# filter_channels = 1 

# filter_above_sectheta = 1.3 

peak_indices = [ [0], [0,1] ]
num_peaks_per_source = [ len( x ) for x in peak_indices ] 


db_names = [ 'alpharun11-19' ]



# strip_coords = 'all'
strip_coords = None





num_dbs = len( db_names )



det_center_x = 1.08 * 25.4
# det_center_x = -1.08 * 25.4

det_center_y = 0.2 * 25.4 # plus minus of this

det_center_z = 4.155 * 25.4


# first are the strips up to which source1 illuminates (not inclusively)
# and second are strips above which source2 illuminates. e.g. source 1
# for det 2 illuminates fstrips 0, ..., 8.

illuminated_fstrips = [ [ 9, 9, 6, 9 ],
                        [ 24, 25, 20, 18 ] ]


def get_secant( sourcenum, detnum, fstrip, bstrip ) :

    fstrip += 1
    bstrip += 1 

    theta1 = np.arctan2( np.sqrt( (det_center_x + ( bstrip - 16.5) * 2 ) ** 2 
                                     + ( det_center_y + ( fstrip - 16.5) * 2 ) ** 2 ),
                                     det_center_z )

    theta2 = np.arctan2( np.sqrt( ( det_center_x + (bstrip - 16.5 ) * 2) ** 2
                                     + ( -det_center_y + ( fstrip - 16.5 ) * 2) ** 2),
                                     det_center_z )

    if sourcenum == 0 :
        if fstrip < illuminated_fstrips[ sourcenum ][ detnum ] :
            return 1 / np.cos( theta1 )
        else :
            return np.nan

    else :
        if fstrip > illuminated_fstrips[ sourcenum ][ detnum ] :
            return 1 / np.cos( theta2 )
        else :
            return np.nan
    

# test secant ( remove )
theta1 = np.arctan2( np.sqrt( (det_center_x + ( np.arange(32) - 16.5) * 2 ) ** 2 
                              + ( det_center_y + ( 0 - 16.5) * 2 ) ** 2 ),
                     det_center_z )

theta2 = np.arctan2( np.sqrt( (det_center_x + ( np.arange(32) - 16.5) * 2 ) ** 2 
                              + ( - det_center_y + ( 0 - 16.5) * 2 ) ** 2 ),
                     det_center_z )



# print( 'theta1 :', theta1 ) 

# print( 'theta2 :', theta2 ) 

# sys.exit(0);

    
det_sectheta = np.zeros( (2,4, 32,32))

for detnum in range(4) : 
    for sourcenum in range(2) :
        for x in range(32) :
            for y in range(32) :
                det_sectheta[ sourcenum, detnum, x, y ] = get_secant(
                    sourcenum, detnum, x, y ) 

print( 'sec theta1: ', det_sectheta[0,0,0] ) 
# sys.exit(0)
                
source_sectheta = det_sectheta



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

    
actual_energies = [ np.array( [ 5123.68, 5168.17 ] ),
                    np.array( [ 5456.3, 5499.03 ] ),
                    np.array( [ 5759.5, 5813.3, 5903.2, 5946.0 ] ) ]




density_si = 2.328 # g / cm^2



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

    energy, stopping_power = np.loadtxt( '../../../data/stopping_power_data/alpha_stopping_power_si.txt',
                                        unpack = 1 )

    # energy = data[0] * 1000
    
    # energy = energy[ energy <= emax ] 
    
    # stopping_power = data[3][ 0 : len( energy ) ]
    energy *= 1000 
    stopping_power *= density_si * 1000 * 100 / 1e9


    # add particular data points of interest to the interpolation
    
    interp = scipy.interpolate.interp1d( energy, stopping_power, kind = 'cubic' )

    
    if plot :

        ax = plt.axes()

        interp_axis = np.linspace( min( energy ),
                                   max( energy ),
                                   100 )
        
        ax.scatter( energy, stopping_power, color='r' )

        ax.plot( interp_axis,
                 interp( interp_axis ),
                 color = 'b' )
        
        plt.show()

        return 1
        

    return interp
    

            




# construct the mean energy loss divide by height at each 









model_params = dl_estimator.deadlayer_model_params( vary_det_deadlayer = 1,
                                                    interp_det_stopping_power = 1,
                                                    interp_source_stopping_powers = 0,
                                                    fstrips_requested = np.arange(2,30),
                                                    bstrips = np.arange( 2, 30 ),
                                                    fix_source_deadlayers = None,
                                                    one_source_constant = 0,
                                                    det_thickness = 0,
                                                    vary_source_thickness = 0 )



db_path = '../../../storage/databases/'

dbs = [ spec.spectrum_db( db_path + db_names[i] ) for i in range( len( db_names ) ) ]

channels = [ dbs[i].read_values( '../../../storage/peak_values/'
                                    + db_names[i]
                                    + '_peak_values.bin' )
             for i in range( len( db_names ) ) ]



num_sources = len( channels[0] ) 



# print( channels )



# if peak_indices is not None :
#     for i in range( num_dbs ) : 
#         for j in range( num_sources ) :

#             channels[i][j] = [ channels[i][j][x]  for x in peak_indices[j] ]





# if filter_channels :
    
#     for d in range( num_dbs ) :
#         for i in range( len( source_geometries[d] ) ):
#             for j in range( len( channels[d][i] ) ) :
            
#                 mask = ( channels[d][i][j].dx > 1 )
#                 channels[d][i][j][ mask] = meas.nan
                    


                    
        
# CONSTRUCT STOPPING POWER INTERPOLATION 

det_stopping_power_interp = construct_si_stopping_power_interpolation()




# INITIAL FIT PARAMETERS 

source_deadlayer_guesses = [ [ 6., 6.], [3.,3.], [15.,15.,15.,15.] ] 
source_stopping_power_interps = [ None ] * 3 
det_deadlayer_guess = 100.0
calibration_coefs_guess = [ 2.0, -260.0 ]


# print( len( channels ) )
# print( len( channels[0] ) ) 
# tmp = channels[ 0 ][1]
# print( len( tmp ) )
# print( len( tmp[0] ) )
# print( len( tmp[0][0] ) ) 



# strip_coords = 'all' 
# strip_coords = [0,2]

for sourcenum in [0,1] :
    for detnum in  range(4) :

        det_sectheta_tmp = [ [ det_sectheta[ sourcenum, detnum ]
                               for i in range( num_sources ) ] ]

        
        source_sectheta_tmp = [ [ det_sectheta[ sourcenum, detnum ]
                                  for i in range( num_sources ) ] ]

        channels_tmp = [ channels[ 0 ][ detnum ] ] 
        
        dl_estimator.estimate_deadlayers( model_params, channels_tmp,
                                          actual_energies,
                                          det_sectheta_tmp,
                                          source_sectheta_tmp,
                                          source_stopping_power_interps,
                                          source_deadlayer_guesses,
                                          det_stopping_power_interp, det_deadlayer_guess,
                                          calibration_coefs_guess,
                                          names = db_names,
                                          strip_coords = strip_coords,
                                          figpath = '../../../storage/current_peaks_vs_sectheta/exp2/det_%d/source_%d'
                                          % ( detnum, sourcenum ) ) 
        




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













