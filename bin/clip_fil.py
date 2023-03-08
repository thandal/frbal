#!/opt/local/bin/python3.7
import logging
import os
import sys

import blimpy as bl
import numpy as np
import pylab as pl

logger = logging.getLogger('CLIPFIL')
logger.setLevel(logging.WARN)
#logger.setLevel(logging.DEBUG)

fil_filename = sys.argv[1]
fil_prefix = fil_filename[:-9]
fil_fullname = '/users/nfairfie/scratch/fil/' + fil_filename

# Read headers
obs = bl.Waterfall(fil_fullname, load_data=False)
file_shape = obs.file_shape
obs.info()

#START_i_s = 0
#START_i_s = 67844096 # Skynet_59902_FRB20220912A_86271_36512_0001.fil
#START_i_s = 28262400  # Skynet_59905_FRB20220912A_86271_36575_0001.fil
#START_i_s = 28286976  # Skynet_59905_FRB20220912A_86271_36575_0001.fil
#START_i_s = 125976576 # Skynet_59905_FRB20220912A_86271_36575_0001.fil
START_i_s = 125998080

#END_i_s = file_shape[0]
END_i_s = START_i_s + 2**14

obs = bl.Waterfall(fil_fullname, t_start=START_i_s, t_stop=END_i_s, max_load=2.0)
D = np.squeeze(obs.data).T

import pickle
pickle.dump(D, open("fil_clip_D.pkl", "wb"))
