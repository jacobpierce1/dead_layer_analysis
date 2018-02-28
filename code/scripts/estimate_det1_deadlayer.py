import numpy as np
import deadlayer_delpes.deadlayer_estimator as dl_estimator




    

# dbs = [ dbmgr.angled ]
# dbs = [ dbmgr.centered, dbmgr.angled, dbmgr.moved, dbmgr.flat ] 
# dbs = [ dbmgr.flat ]
# dbs = dbmgr.all_dbs 
# dbs = [ dbmgr.angled ]

# dbs = dbmgr.det3_dbs
dbs = [ dbmgr.det3_moved ] 

# source_indices = [ [0,1], [0,1], [1] ]
source_indices = [ [0,1], [0,1], [1] ]



subtitle = ''

fstrips = np.arange( 2, 30 )


model_params = dl_estimator.dead_layer_model_params( vary_det_deadlayer = 0,
                                                     quadratic_source = 0,
                                                     quadratic_det = 0,
                                                     calibrate_each_pixel = 0,
                                                     interp_stopping_power = 1,
                                                     mu = 0,
                                                     average_over_source = 0,
                                                     pulse_height_defect = 0,
                                                     fstrips_requested = fstrips,
                                                     bstrips = np.arange( 2, 30 ),
                                                     different_fstrip_deadlayers = 0,
                                                     different_pixel_deadlayers = 0,
                                                     fix_source_deadlayers = None,
                                                     one_source_constant = 1 )





dl_estimator.linear_calibration_on_each_x_strip( dbs, source_indices,
                                                 model_params,
                                                 cut_high_sectheta = 0,
                                                 annotate = 0,
                                                 # view_pixel = [ dbmgr.moved, 5 ],
                                                 subtitle = 'Det 3: Source Params Extrapolated from Det 1\n',
                                                 reset_angles = None,
                                                 residual_scatter_plot = 0,
                                                 plot_3d = 1 ) #  ,
# savefig_dir = '../../../deadlayer_analysis_paper/images/' )




