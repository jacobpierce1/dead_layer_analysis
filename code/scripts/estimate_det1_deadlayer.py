import numpy as np

import deadlayer_helpers.deadlayer_estimator as dl_estimator
import deadlayer_helpers.geometry as geom

import jspectroscopy as spec
import libjacob.jutils as jutils
import libjacob.jmeas as meas 

import scipy.interpolate


filter_channels = 1 


# si_interp_energies = np.array( [ 5.050E+00, 5.090E+00, 5.130E+00, 5.170E+00,
#                                  5.210E+00, 5.250E+00, 5.290E+00, 5.330E+00,
#                                  5.370E+00, 5.410E+00, 5.450E+00, 5.490E+00,
#                                  5.530E+00, 5.570E+00, 5.610E+00, 5.650E+00,
#                                  5.690E+00, 5.730E+00, 5.770E+00, 5.810E+00,
#                                  5.81310,] )

# si_interp_stopping_powers = np.array( [ 6.138E+02, 6.107E+02, 6.075E+02, 6.044E+02,
#                                         6.014E+02, 5.983E+02, 5.953E+02, 5.924E+02,
#                                         5.894E+02, 5.866E+02, 5.837E+02, 5.809E+02,
#                                         5.781E+02, 5.753E+02, 5.726E+02, 5.699E+02,
#                                         5.672E+02, 5.645E+02, 5.619E+02, 5.593E+02,
#                                         5.591E+02 ] )


actual_energies = [ np.array( [ 5123.68, 5168.17 ] ),
                    np.array( [ 5456.3, 5499.03 ] ),
                    np.array( [ 5759.5, 5813.3, 5903.2, 5946.0 ] ) ]


# 5.12368
# 5.16817
# 5.4563
# 5.49903
# 5.7595
# 5.8133
# 5.9032
# 5.9460
# 5.0
# 5.1
# 5.2
# 5.3
# 5.4
# 5.5
# 5.6
# 5.7
# 5.8
# 5.9


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

    energy, stopping_power = np.loadtxt( '../../data/stopping_power_data/alpha_stopping_power_si.txt',
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
    



    


# dbs = dbmgr.all_dbs 
# # dbs = [ dbmgr.angled ]

# source_indices = [ [0,1], [0,1], [1] ]
# # source_indices = [ [0,1], [0,1], [0,1] ]


# subtitle = ''

# fstrips = np.arange( 2, 30 )




#  # next, read the data that will go into the regression. 
# if model_params.mu : 
#     mu_matrices = { db.name : db.get_all_mu_grids( 1 )
#                      for db in dbs }
     
#     peak_matrices = mu_matrices
     
# else:
#     peak_matrices = { db.name : db.get_all_peak_grids( 1 )
#                           for db in dbs }
     
# # filter out bad data 
# for db_name in [ db.name for db in dbs ] :
#     for i in range(3) :
#         for j in range( 2 ) :
#             mask = ( ( peak_matrices[db_name][i][j].dx > 2 ) | 
#                      ( peak_matrices[db_name][i][j].dx < 0.01 ) )
#             peak_matrices[db_name][i][j][ mask ] = meas.nan
            
# if reset_angles is None :
#     reset_angles = not model_params.average_over_source
                
# secant_matrices = geom.get_secant_matrices( compute_source_sectheta = 1,
#                                             average_over_source = model_params.average_over_source,
#                                             reset = reset_angles )

    
#     # cut out large sec theta 
#     if cut_high_sectheta : 
#         for key, val in secant_matrices.items() :
#             mask = ( val[0].x >= 1.45 )
#             secant_matrices[key][0][mask] = meas.nan


            



model_params = dl_estimator.deadlayer_model_params( vary_det_deadlayer = 0,
                                                    interp_det_stopping_power = 1,
                                                    interp_source_stopping_powers = 0,
                                                    fstrips_requested = np.arange(2,30),
                                                    bstrips = np.arange( 2, 30 ),
                                                    fix_source_deadlayers = None,
                                                    one_source_constant = 0,
                                                    det_thickness = 1 )

db_names = [ 'moved', 'angled' ] 
# db_names = [ 'angled' ] # , 'flat' ]
# db_names = [ 'angled', 'moved' ]

num_dbs = len( db_names ) 

db_path = '../../storage/databases/'
dbs = [ spec.spectrum_db( db_path + db_names[i] ) for i in range( len( db_names ) ) ]

channels = [ dbs[i].read_mu_values( '../../storage/mu_values/'
                                    + db_names[i]
                                    + '_mu_values.bin' )
             for i in range( len( db_names ) ) ] 


num_sources = len( channels[0] ) 

peak_indices = [ [1,2], [1,2], [1,3,4] ]
# peak_indices = [ [1,2], [1,2], [2] ]

if peak_indices is not None :
    for i in range( num_dbs ) : 
        for j in range( num_sources ) :

            channels[i][j] = [ channels[i][j][x]  for x in peak_indices[j] ]

            if filter_channels : 
                for x in range( len( peak_indices[j] ) ) :
                    mask = channels[i][j][x].dx > 1
                    channels[i][j][x][mask] = meas.nan


secant_matrices = geom.get_secant_matrices( compute_source_sectheta = 1, reset = 1 )

pu_240_sources = [ dl_estimator.source_geom_data( secant_matrices['pu_240'][0],
                                                  secant_matrices['pu_240'][1], 0 ) ] * len( peak_indices[0] )

pu_238_sources = [ dl_estimator.source_geom_data( secant_matrices['pu_238_' + name ][0],
                                                  secant_matrices['pu_238_' + name][1], 1 if name == 'angled' else 0 )
                   for name in db_names ]

cf_249_sources = [ dl_estimator.source_geom_data( secant_matrices['cf_249'][0],
                                                  secant_matrices['cf_249'][1], 0 ) ] * len( peak_indices[2] )

source_geometries = [ [0] * num_sources  for i in range( num_dbs ) ] 

for i in range( num_dbs ) :

    sources = [ 'pu_240', 'pu_238_' + db_names[i], 'cf_249' ]
    is_angled = [ 0, db_names[i] == 'angled', 0 ] 
    
    for j in range( num_sources ) : 
   
        source_geometries[i][j] = dl_estimator.source_geom_data( secant_matrices[ sources[j] ][0],
                                                                 secant_matrices[ sources[j] ][1],
                                                                 is_angled[j] )


source_stopping_power_interps = [ None ] * 3 

source_deadlayer_guesses = [ [ 6., 6.], [3.,3.], [15.,15.,15.,15.] ] 

det_stopping_power_interp = construct_si_stopping_power_interpolation()

det_deadlayer_guess = 100.0

calibration_coefs_guess = [ 2.0, -260.0 ]


strip_coords = 'all' 
# strip_coords = None
# strip_coords = [0,2]

dl_estimator.estimate_deadlayers( model_params, channels, actual_energies,
                                  source_geometries, source_stopping_power_interps,
                                  source_deadlayer_guesses,
                                  det_stopping_power_interp, det_deadlayer_guess,
                                  calibration_coefs_guess,
                                  names = db_names,
                                  strip_coords = strip_coords ) 




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




