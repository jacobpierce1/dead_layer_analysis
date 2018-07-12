import jspectroscopy as spec


db_name = 'full_bkgd_tot' 


db = spec.spectrum_db( '../../storage/databases/' + db_name )

success_rate = db.compute_fit_success_rate()
print( success_rate )

db.plot_success_pixels()
    
db.disconnect() 
