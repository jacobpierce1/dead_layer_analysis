import jspectroscopy as spec




# db_names = [ 'moved' ]
# db_names = [ 'angled', 'moved', 'centered', 'flat' ] # , 'det3_cent', 'det3_moved' ]

# db_names = [ 'det3_cent', 'det3_moved' ] 

db_names = [ 'alpharun11-19', 'alpharun20-30' ] 


for name in db_names :
    
    db = spec.spectrum_db( '../../storage/databases/' + name ) 

    path = '../../storage/peak_values/' + name + '_peak_values.bin'

    db.write_peak_values( path, 1 )

    db.disconnect() 

    # mu_values = test_db.read_mu_values( path )

    # print( mu_values ) 
