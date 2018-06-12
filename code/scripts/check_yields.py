import jspectroscopy as spec


db_names = [ 'moved', 'angled', 'centered', 'flat', 'det3_cent', 'det3_moved',
             'alpharun11-19', 'alpharun20-30'  ]

for name in db_names :

    db = spec.spectrum_db( '../../storage/databases/' + name )


    success_rate = db.compute_fit_success_rate()

    print( name ) 
    print( '%.2f\n' % success_rate ) 
    
    db.disconnect() 
