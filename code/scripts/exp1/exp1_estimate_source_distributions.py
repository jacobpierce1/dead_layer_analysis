# in this script, we will use integrated hit hit counts from each source
# on each pixel to estimate source distributions from sources placed
# below the detector. given the massive amount of information (32^2 data points )
# for each source, not to mention more data obtained from aggregating datasets,
# the source height distributions can be heavily constrained.

import jutils 
import scipy.integrate
import matplotlib.pylot as plt
