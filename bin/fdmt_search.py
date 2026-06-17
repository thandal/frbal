#!/opt/local/bin/python3.7
import os
import sys

import blimpy as bl
import numpy as np
import matplotlib
matplotlib.use('Agg')  # cluster nodes run this headless over ssh -- no X display
import matplotlib.pyplot as plt

# FDMT kernel lives in fdmt.py; per-channel normalization in preprocess.py.
# Both must be deployed alongside this script (same dir / ~/bin).
from fdmt import FDMT
from preprocess import normalize_robust
from detect import boxcar_search

fil_filename = sys.argv[1]
fil_prefix = fil_filename[:-9]
fil_fullname = '/users/nfairfie/scratch/fil/' + fil_filename
# Read headers
obs = bl.Waterfall(fil_fullname, load_data=False)
file_shape = obs.file_shape
obs.info()

dt = obs.header['tsamp']
freqs = obs.get_freqs()
f_min = min(freqs)
f_max = max(freqs)
f_start = np.argmin(np.abs(freqs - f_min))
f_end = np.argmin(np.abs(freqs - f_max))
if f_start > f_end: f_start, f_end = f_end, f_start
N_f = f_end - f_start + 1
print('frequency selection:', f_min, f_max, f_start, f_end, N_f)

DM_max = 800   # Maximum dispersion to compute
DM_min = 10.0  # Minimum dispersion of interest
ds_max = int(DM_max * 4148.808 * (f_min**-2 - f_max**-2) / dt)  # Max sample bin shift, from DM_max
ds_min = int(DM_min * 4148.808 * (f_min**-2 - f_max**-2) / dt)  # Min sample bin shift, from DM_min

N_s = 2**15       # Number of samples per chunk
overlap_s = 2**13  # Overlap chunks by this amount (must be >= ds_max for gap-free coverage)
if ds_max > overlap_s:
  print(f'WARNING: ds_max ({ds_max}) > overlap_s ({overlap_s}): high-DM bursts '
        f'near chunk boundaries may be missed. Increase overlap_s.')

os.makedirs(f'/users/nfairfie/scratch/results/{fil_prefix}', exist_ok=True)

for i_s in range(0, file_shape[0], N_s - overlap_s):
  result_filename = f'/users/nfairfie/scratch/results/{fil_prefix}/{fil_filename}_{i_s:010}.png'
  if os.path.exists(result_filename): continue  # Note: race condition with other workers right here.
  open(result_filename, 'w')  # write a placeholder to claim this chunk

  print(f'processing chunk starting at sample {i_s}...')
  obs = bl.Waterfall(fil_fullname, t_start=i_s, t_stop=(i_s + N_s), max_load=2.0)
  D = np.squeeze(obs.data).T
  D = D[f_start:(f_end + 1)]
  chan_freqs = freqs[f_start:(f_end + 1)]
  print(D.shape)

  # Robust per-channel normalization: subtract median, divide by MAD-based noise
  # sigma, then clip impulsive RFI. This whitens the band so the FDMT's sum across
  # channels is the matched-filter-optimal (max-S/N) combination, and -- unlike the
  # old min-max -- is immune to a bright burst or RFI spike compressing the channel.
  D = normalize_robust(D, clip_sigma=5.0)
  #D[-5:, :] = 0  # GB-20 256-ch data: these edge channels are supposed to be ~0
  #               # (bandpass) but sometimes really aren't.

  # FDMT requires channels in ascending frequency order (channel 0 == f_min),
  # validated in fdmt_validate.py. Sort explicitly instead of assuming the file's
  # channel order -- orientation-proof; replaces the old, ambiguous D[::-1].
  order = np.argsort(chan_freqs)
  DMT = FDMT(D[order], f_min, f_max, ds_max, 'float32')  # Compute the DMT
  print(DMT.shape)
  DMT = DMT[ds_min:, ds_max:]  # Crop off low DMs, and the first ds_max (edge-contaminated) samples

  # Boxcar matched-filter width search: per (DM, time) cell, the best S/N across
  # boxcar widths (robust median/MAD z-score per row+width). Recovers ~sqrt(W) of
  # S/N for width-W pulses that a single-sample statistic would miss, and gives a
  # robust detection threshold. `best` is the peak candidate.
  detect, best = boxcar_search(DMT, max_width=64)

  if (best['snr'] > 6.0):
    DM_best = DM_min + (DM_max - DM_min) * best['i_dm'] / max(detect.shape[0] - 1, 1)
    t_best = (i_s + ds_max + best['i_t']) * dt
    fig = plt.figure(figsize=(16, 12))
    plt.imshow(detect, cmap='magma', origin='lower', aspect='auto',
        extent=((i_s + ds_max) * dt, (i_s + ds_max + detect.shape[1]) * dt, DM_min, DM_max))
    plt.colorbar(label='matched-filter S/N')
    plt.plot(t_best, DM_best, 'c+', ms=20, mew=2)  # mark the peak candidate
    plt.ylabel('DM (pc/cm^3)')
    plt.xlabel('time (s)')
    plt.title(f"{fil_filename}_{i_s:010} S/N={best['snr']:.1f} "
              f"DM={DM_best:.0f} t={t_best:.2f}s width={best['width']}")
    plt.savefig(result_filename)
    plt.close(fig)
  else:
    open(result_filename, 'w').write(f"snr={best['snr']:.2f} width={best['width']} i_dm={best['i_dm']}")
