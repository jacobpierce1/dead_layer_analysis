import dill

# this is a function for the general case, e.g. if there were multiple datasets
# storage path = directory where databases are located (this will make life easier
# if i gave you more data, then that argument is the same )
# name = name of the dataset (full_bkgd_tot)
# param_str = name of parameter, e.g. secant_matrices

def load_dill( storage_path, name, param_str ) :
    path = '%s/%s/%s.%s.dill' % ( storage_path, name, name, param_str ) 
    with open( path, 'rb' ) as f :
        return dill.load( f ) 


def my_load_dill( param_str ) :
    return load_dill( '../../storage/', 'full_bkgd_tot', param_str ) 
    


# secant_matrices
# key: [ source, detidx, fstrip, bstrip ]
# source: 0 = gd, 1 = cm
# detidx: 0 to 3
# fstrip: 0 to 31
# bstrip: 0 to 31

secant_matrices = my_load_dill( 'secant_matrices' )
print( 'secant_matrices' ) 
print( len( secant_matrices ) )
print( len( secant_matrices[0] ) )
print( len( secant_matrices[0][0] ) )
print( len( secant_matrices[0][0][0] ) )
print( secant_matrices[0][0][0] )



# tau1
# key: [ source, det, fstrip, bstrip ] 

tau1 = my_load_dill( 'tau1' )
print( '\n\ntau1' )
print( len( tau1 ) )
print( len( tau1[0] ) )
print( len( tau1[0][0] ) )
print( len( tau1[0][0][0] ) )
print( tau1[0][0][1] )
print( tau1[0][0][1].x )
print( tau1[0][0][1].dx ) 



# peak vals
# key: [ source, peaknum, det, fstrip, bstrip ]
# for params that are not independent of peak, the second index is the peak number
# starting at 0 for each source

peaks = my_load_dill( 'peaks' )
print( '\n\npeaks' ) 
print( len( peaks ) )
print( len( peaks[0] ) )
print( len( peaks[0][0] ) )
print( len( peaks[0][0][0][0] ) )
print( peaks[0][0][0][1] )


# access x and dx the same way as before 
