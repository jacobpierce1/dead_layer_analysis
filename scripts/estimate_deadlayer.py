# this script reads in dataframes containing: 1. the penetration angle
# of each pixel ( constructed in deadlayer_geometry.py ), 2. the mu
# values of each peak fit ( originally detected in parse_all_data.py
# and read into a matrix in deadlayer_analysis.py ) it then reads in a
# file of stopping power data and interpolates it. using all this data
# we make estimates of both the source dead layer depths and
# the pixel dead layer depth.




import libjacob.jacob_pyplot as jplt
import libjacob.error

import deadlayer_helpers.stopping_power_interpolation as stop
import deadlayer_helpers.deadlayer_geometry

import matplotlib.pyplot as plt
import pandas as pd
import numpy as np



sources = deadlayer_helpers.s


def _main():

    # read in all stopping power data.
    stop_interp_si = stop.stopping_power_interpolation.from_nist_file(
        file_si, density_si )

    

_main()
