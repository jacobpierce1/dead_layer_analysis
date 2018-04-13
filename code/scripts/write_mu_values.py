import jspectroscopy as spec



db_names = [ 'moved' ]
# db_names = [ 'angled', 'moved', 'centered', 'flat', 'det3_cent', 'det3_moved' ]

for name in db_names :

    db = spec.spectrum_db( '../../storage/databases/' + name ) 

    path = '../../storage/mu_values/' + name + '_mu_values.bin'

    db.write_mu_values( path ) 

    # mu_values = test_db.read_mu_values( path )

    # print( mu_values ) 
