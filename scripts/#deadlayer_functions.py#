import array
import numpy as np
from scipy import special
import peakdetect.peakdetect as peakdetect
import sql_db_manager 

# extract next row from the file. file assumed to already be opened.
def read_doubles_from_bin( f, buf ):
    bytes = f.read(8*2 )  # since we write 2 doubles for each event: efront and eback.
    if( not bytes ): return 0;
    buf[:] = array.array('d', bytes) 
    return 1

    

# # construct histogram count arrays from the file
# def construct_histo_arrays( f, efront_histo, eback_histo ):

    
#     buf = array.array( 'd', [0.0, 0.0] )  # data gets read into here.

#     while( read_doubles_from_bin( f, buf ) ):
#         efront = buf[0]
#         eback = buf[1]
        
#         # used this to verify that read was successful.
#         # print efront
#         # print eback
#         # print int(efront)
#         # print int(eback)
    
#         efront_histo[ int( efront ) ] += 1
#         eback_histo[ int( eback ) ]  += 1


# open binary infile and populate the array efront_histo as a histogram 
def construct_histo_array( infile, efront_histo ):
    try:
        with open( infile, "rb" ) as f:
            buf = array.array( 'd', [0.0, 0.0] )  # data gets read into here.
            while( read_doubles_from_bin( f, buf ) ):
                efront_histo[ int( buf[0] ) ] += 1
                
    except IOError:
        return 0

    return 1


## write the pf parameters, chisq of the fit, etc. to a file
## also extract information from the fitfunc where relevant.
#write_info_to_file( fitfunc, pf, chisq ):
#    return
#    




# print message to the console and write it to a log file.
def log_message( msg, logfile ):
    
    if log_message.first_msg or not os.path.exists( logfile ):
        log_message.first_msg = 0
        mode = 'w'
    else:
        mode = 'a'
        
    with open( logfile, mode ) as log:
        log.write( msg + '\n' )
    
    print msg
log_message.first_msg = 1







# detect the positions of the 5 highest peaks. peak_positions is an array of tuples (peakpos, value)
# https://stackoverflow.com/questions/6910641/how-to-get-indices-of-n-maximum-values-in-a-numpy-array
# even if we found all the fits, for simplicity we still do the peak detection since it is the easiest
# way to deteremine a good location for the plot boundaries. 

def get_5_peak_positions( our_peaks, efront_histo, sql_conn=None, pixel_coords=None, check_warnings=0 ):

    # peakdetect returns 2 tuples: positions and counts of peaks
    peak_positions = peakdetect.peakdetect( efront_histo, lookahead=10 )[0]

    # if 5 peaks are not there then we have a problem. this will have to be handled 
    # eventually, for now we just log it and worry about it later.
    if( len(peak_positions) < 5 ):
        return 0
    
    # now find the 5 largest and sort by x position.
    ind = np.argpartition( [z[1] for z in peak_positions ], -5 )[-5:]
    our_peaks[:] = [ peak_positions[z][0] for z in sorted(ind) ]
    
    # do a check on peak values: energy differenc for the pairs should be constant. no check necessary on the 
    # largest peak since it always dominates, the only potential problem is really the smallest peak in the 
    # lowest energy pair.
    if check_warnings:
        if( abs(23 - (our_peaks[1] - our_peaks[0] ) ) > 15 ):
            log_message( "WARNING (%d,%d,*): invalid peak suspected for pair 1. " % (pixel_coords[0], pixel_coords[1]) ) 
    
        if( abs(23 - (our_peaks[3] - our_peaks[2] ) ) > 15 ):
            log_message( "WARNING (%d,%d,*): invalid peak suspected for pair 2. " % (pixel_coords[0], pixel_coords[1]) )
    
    return 1







# in this case you supply all args as fit params, even the det ones which should be fixed. 
def fitfunc_alpha_free_det_params( p, x ): 
    return alpha_fit( p[0],p[1],p[2],p[3],p[4],p[5], x ) 

#    return np.array( map( lambda z: alpha_fit( *p, z ), x ) )
    #return np.array( map( lambda z: alpha_fit( p[0],p[1],p[2],p[3],p[4],p[5],z ), x ) )
    # return alpha_fit( p[0],p[1],p[2],p[3],p[4],p[5] , x )  # expand array as function args 
    


# 8 total params for p.
# first 4 params: A1, mu1 ,A2, u2
def fitfunc_2_alpha_peaks( p, x ):
    return map( lambda y: ( alpha_fit( p[0],p[1],p[4],p[5],p[6],p[7], y) +
                alpha_fit( p[2],p[3],p[4],p[5],p[6],p[7], y) ), x ) 


# return an inline of fitfunc_n_alpha_peaks given the number of peaks to fit.
# the first 4 params are the detector params, and the last 2*n are (A1, u1, ... )
# previously fitfunc_n_alpha_peaks_abstract
def n_fitfuncs_abstract( fitfunc, n ):
    return lambda p, x: fitfunc( n, p, x )


# fit n alpha peaks given vector p and array x 
# p format: sigma, eta, tau1, tau2, A1, mu1, ..., A_n, mu_n
def fitfunc_n_alpha_peaks( n, p, x ):
    return np.sum( [ map( lambda y: alpha_fit(p[4+2*i],p[4+2*i+1],p[0],p[1],p[2],p[3], y), x )  for i in range(n) ], axis=0 )  
    
    
# fit n alpha peaks given vector p and array x 
# p format: sigma, tau, A1, mu1, ..., A_n, mu_n
def fitfunc_n_alpha_peaks_eta1( n, p, x ):
    return np.sum( [ map( lambda y: alpha_fit_eta1(p[2+2*i],p[2+2*i+1],p[0],p[1], y), x )  for i in range(n) ], axis=0 )  


# same as alpha_fit but with eta fixed at one. this means we are taking only one of the terms. 
# experimented with this after noting that the fits seem indep. of eta, which can cause trouble.
def alpha_fit_eta1( A, mu, sigma, tau, x ):
    logtmp = (x-mu)/tau + sigma**2.0/(2*tau**2.0) + np.log( special.erfc( (x-mu)/sigma + sigma/tau) / np.sqrt(2) ) 
    return A/(2.0*tau) * np.exp( logtmp )

 
# reference: equation 10 in Bortels 1987  
# this function is meant to be applied to scalar x, not list
def alpha_fit( A, mu, sigma, eta, tau1, tau2, x ):        
    #if eta < 0 or eta > 1:
    #    return -1000
        
    # prevent overflow by computing logs and then exponentiating, at the expense of some
    # floating pt error. logtmpz is the log of the 2 analagous terms in the integrand.
    tauz = [tau1, tau2]
    logtmpz = [ (x-mu)/tau + sigma**2.0/(2*tau**2.0) + np.log( special.erfc( (x-mu)/sigma + sigma/tau) / np.sqrt(2) ) for tau in tauz ]
    return ( A/2.0 )  * ( (1-eta)/tau1 * np.exp(logtmpz[0])
                            + (eta/tau2) * np.exp(logtmpz[1])   ) 


