import numpy as np
import matplotlib.pyplot as plt
import jspectroscopy as spec 

analysis_mgr = spec.dssd_analysis_manager( 'full_bkgd_tot', '../../../storage', (32,32), [1,1,1] )

indep_dl_sums = analysis_mgr.load_dill( 'dl_estimates_indep' )


f, axarr = plt.subplots( 4, 2, figsize = ( 10, 10 ) )

print( indep_dl_sums )

for d in range(4) :
    for i in range(2) :
        tmp = indep_dl_sums[i][0][d]
        axarr[d,i].errorbar( range(32), tmp.x, tmp.dx, ls = 'none' ) 

plt.show() 
