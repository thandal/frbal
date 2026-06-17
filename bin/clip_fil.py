#!/opt/local/bin/python3.7
import pickle
import sys

import blimpy as bl
import numpy as np

fil_filename = sys.argv[1]

# Read headers
obs = bl.Waterfall(fil_filename, load_data=False)
file_shape = obs.file_shape
obs.info()

START_i_s = int(sys.argv[2])
DURATION_i_s = int(sys.argv[3])
END_i_s = START_i_s + DURATION_i_s

obs = bl.Waterfall(fil_filename, t_start=START_i_s, t_stop=END_i_s, max_load=2.0)
D = np.squeeze(obs.data).T
pickle.dump(D, open(f"{fil_filename}_{START_i_s}_{END_i_s}.pkl", "wb"))
