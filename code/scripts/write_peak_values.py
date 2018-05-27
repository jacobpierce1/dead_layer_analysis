import jspectroscopy as spec
# import libjacob.jutils 
from libjacob.jutils import isint, time_estimator



# db_names = [ 'moved' ]
# db_names = [ 'angled', 'moved', 'centered', 'flat' ] # , 'det3_cent', 'det3_moved' ]

db_names = [ 'det3_cent', 'det3_moved' ] 


time_estimator = time_estimator( len(db_names), 20 )

for name in db_names :

    time_estimator.update()
    
    db = spec.spectrum_db( '../../storage/databases/' + name ) 

    path = '../../storage/peak_values/' + name + '_peak_values.bin'

    db.write_peak_values( path ) 

    # mu_values = test_db.read_mu_values( path )

    # print( mu_values ) 
