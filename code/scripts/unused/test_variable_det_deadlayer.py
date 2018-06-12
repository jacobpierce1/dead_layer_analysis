import numpy as np 


# return deltaX given the energy loss. this function will
# be inverted using newton's method.

def deltaE_to_deltaX( Ef, E0, S0, Sprime ) :

    deltaE = E0 - Ef 
    
    return ( deltaE * Sprime - ( S - E0 * sprime ) * np.log( S / ( S - deltaE * Sprime ) ) ) / Sprime ** 2


# def tmp( E, E_0, E, S0, Sprime ) :

#     return ( E * Sprime - np.log( ( E - E0 ) * 
