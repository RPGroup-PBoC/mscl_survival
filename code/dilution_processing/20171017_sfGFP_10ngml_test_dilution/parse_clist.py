import glob
import numpy as np
import matplotlib.pyplot as plt
import mscl_utils as mscl
import pandas as pd
import scipy.io
import matlab.engine as matlab
try:
    eng = jmatlab.connect_matlab()
except:
    eng = matlab.start_matlab()


colors = mscl.set_plotting_style()

#
DATE = 20171017
BASENAME = 'sfGFP_10ngmL'
STRAIN = 'HG105 galK::sfGFP'
ATC_CONC = 10

# %% Read in the clist
data_dir = '../../../data/images/{0}_{1}_dilution/'.format(DATE, BASENAME)
clist_mat = scipy.io.loadmat('{}growth/xy01/clist.mat'.format(data_dir),
                             char_as_strings=True, squeeze_me=True)
clist_mat
scipy.io.load
clist_mat['def']
scipy.io.loadmat
