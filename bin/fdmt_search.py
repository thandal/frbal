#!/opt/local/bin/python3.7
import logging
import os
import sys

import blimpy as bl
import numpy as np
import pylab as pl

logger = logging.getLogger('FDMT')
logger.setLevel(logging.WARN)
#logger.setLevel(logging.DEBUG)

# Original code from https://github.com/iansbrown/FRB-FDMT-Search, mostly written by Bruce Wu with modifications by Ian Brown
# See also the fdmt.m codes under https://webhome.weizmann.ac.il/home/eofek/matlab/FunList.html

def FDMT_initialization(Image, f_min, f_max, maxDT, dtype):
   """
   Data initialization is done prior to the first FDMT iteration.
   Input:
     Image - power matrix I(f,t)
     f_min, f_max - are the base-band begin and end frequencies.
       The frequencies can be entered in both MHz and GHz, units are factored out in all uses.
     ds_max - the maximal delay (in sample bins) of the maximal dispersion.
       Appears in the paper as N_{\Delta}
       A typical input is maxDT = N_f
     dtype - To naively use FFT, one must use floating point types.  # Use float64...
       Due to casting, use either complex64 or complex128.
   Output:
     3d array, with dimensions [N_f, N_d0, Nt]
   where N_d0 is the maximal number of bins the dispersion curve travels at one frequency bin
  
   For details, see Algorithm 1 in Zackay & Ofek (2014)
   """
   # See Equations 17 and 19 in Zackay & Ofek (2014)
   N_f, N_s = Image.shape
   deltaF = (f_max - f_min) / float(N_f)
   deltaT = int(np.ceil((ds_max - 1) * (f_min**-2 - (f_min + deltaF)**-2) / (f_min**-2 - f_max**-2)))
   Output = np.zeros([N_f, deltaT+1, N_s], dtype)
   Output[:, 0, :] = Image
   for i_dT in range(1, deltaT + 1):
     Output[:, i_dT, i_dT:] = Output[:, i_dT - 1, i_dT:] + Image[:, :-i_dT]
   return Output
    
def FDMT_iteration(Input, maxDT, F, f_min, f_max, iteration_num, dtype):
  """
  Input: 
    Input - 3d array, with dimensions [N_f, N_d0, Nt]
    f_min,f_max - are the base-band begin and end frequencies.
      The frequencies can be entered in both MHz and GHz, units are factored out in all uses.
    maxDT - the maximal delay (in time bins) of the maximal dispersion.
      Appears in the paper as N_{\Delta}
      A typical input is maxDT = N_f
    dtype - To naively use FFT, one must use floating point types.  # Use float64...
      Due to casting, use either complex64 or complex128.
    iteration num - Algorithm works in log2(Nf) iterations, each iteration changes all the 
      sizes (like in FFT)
  Output: 
    3d array, with dimensions [N_f/2, N_d1, Nt]
      where N_d1 is the maximal number of bins the dispersion curve travels at one output frequency band
  
  For details, see Algorithm 1 in Zackay & Ofek (2014)
  """
  input_dims = Input.shape
  output_dims = list(input_dims)
  deltaF = 2**(iteration_num) * (f_max - f_min) / float(F)
  dF = (f_max - f_min)/float(F)
  # the maximum deltaT needed to calculate at the i'th iteration
  deltaT = int(np.ceil((maxDT-1) *(1./f_min**2 - 1./(f_min + deltaF)**2) / (1./f_min**2 - 1./f_max**2)))
  logger.debug(f"deltaT = {deltaT}")
  logger.debug(f"N_f = {F/2.**(iteration_num)}")
  logger.debug(f'input_dims {input_dims}')
  output_dims[0] = int(output_dims[0] / 2)
  output_dims[1] = deltaT + 1
  logger.debug(f'output_dims {output_dims}')
  Output = np.zeros(output_dims, dtype)
  # No negative D's are calculated => no shift is needed
  # If you want negative dispersions, this will have to change to 1+deltaT,1+deltaTOld
  # Might want to calculate negative dispersions when using coherent dedispersion, to reduce the number of 
  # trial dispersions by a factor of 2 (reducing the complexity of the coherent part of the hybrid)
  ShiftOutput = 0
  ShiftInput = 0
  T = output_dims[2]
  F_jumps = output_dims[0]
  # For some situations, it is beneficial to play with this correction.
  # When applied to real data, one should carefully analyze and understand the effect of this
  # correction on the pulse one is looking for (especially if convolving with a specific pulse profile)
  if iteration_num > 0:
    correction = dF/2.
  else:
    correction = 0
  for i_F in range(F_jumps):
    f_start = (f_max - f_min)/float(F_jumps) * (i_F) + f_min
    f_end = (f_max - f_min)/float(F_jumps) *(i_F+1) + f_min
    f_middle = (f_end - f_start)/2. + f_start - correction
    # it turned out in the end, that putting the correction +dF to f_middle_larger 
    # (or -dF/2 to f_middle, and +dF/2 to f_middle larger)
    # is less sensitive than doing nothing when dedispersing a coherently dispersed pulse.
    # The confusing part is that the hitting efficiency is better with the corrections (!?!).
    f_middle_larger = (f_end - f_start)/2 + f_start + correction
    deltaTLocal = int(np.ceil((maxDT-1) *(1./f_start**2 - 1./(f_end)**2) / (1./f_min**2 - 1./f_max**2)))
    for i_dT in range(deltaTLocal+1):
      dT_middle = round(i_dT * (1./f_middle**2 - 1./f_start**2)/(1./f_end**2 - 1./f_start**2))
      dT_middle_index = dT_middle + ShiftInput
      dT_middle_larger = round(i_dT * (1./f_middle_larger**2 - 1./f_start**2)\
                               /(1./f_end**2 - 1./f_start**2))            
      dT_rest = i_dT - dT_middle_larger
      dT_rest_index = dT_rest + ShiftInput
      i_T_min = 0
      i_T_max = dT_middle_larger	
      Output[i_F,i_dT + ShiftOutput,i_T_min:i_T_max] = Input[2*i_F, dT_middle_index,i_T_min:i_T_max]
      i_T_min = dT_middle_larger
      i_T_max = T
      Output[i_F,i_dT + ShiftOutput,i_T_min:i_T_max] = Input[2*i_F, dT_middle_index,i_T_min:i_T_max] \
        + Input[2*i_F+1, dT_rest_index,i_T_min - dT_middle_larger:i_T_max-dT_middle_larger]
  return Output

def FDMT(Image, f_min, f_max, ds_max, dtype):
    """
    The Fast discrete Dispersion Measure Transform (FDMT) algorithm.
    Input: 
      Input power matrix I(f, s)
        dimensions (N_f, N_s), N_f must be a power of 2
      f_min,f_max - the base-band begin and end frequencies in MHz
      DM_max - the maximal dispersion measure (DM) to compute
        ds_max (the max delay in samples) is derived from DM_max
      #maxDT - the maximal delay (in time bins) of the maximal dispersion
      #  Appears in the paper as N_{\Delta}
      #  A typical input is maxDT = N_f
      dtype - a valid numpy dtype  # Use float64...
        Recommended: either int32, or int64.
    Output: 
      The dispersion measure transform of the Input matrix.
        The output dimensions are [Input.shape[1], DT_max]
    For details, see algorithm 1 in Zackay & Ofek (2014)
    """
    N_f, N_s = Image.shape
    f = int(np.log2(N_f))
    if f != np.log2(N_f):
      raise NotImplementedError(f'Input frequency channel dimension ({N_f}) must be a power of 2')
    # Initialize
    State = FDMT_initialization(Image, f_min, f_max, ds_max, dtype)
    logger.debug('FDMT initialized')
    # Iterate
    for i_t in range(1, f + 1):
      State = FDMT_iteration(State, ds_max, N_f, f_min, f_max, i_t, dtype)
    return np.squeeze(State)

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
#f_min = 400
#f_max = 600
f_start = np.argmin(np.abs(freqs - f_min))
f_end = np.argmin(np.abs(freqs - f_max))
if f_start > f_end: f_start, f_end = f_end, f_start
N_f = f_end - f_start + 1
print('frequency selection:', f_min, f_max, f_start, f_end, N_f)

DM_max = 800  # Maximum dispersion to compute
DM_min = 10.0  # Minimum dispersion of interest
ds_max = int(DM_max * 4148.808 * (f_min**-2 - f_max**-2) / dt)  # Max sample bin shift from DM_max
ds_min = int(DM_min * 4148.808 * (f_min**-2 - f_max**-2) / dt)  # Max sample bin shift from DM_max

N_s = 2**15  # Number of samples per chunk
overlap_s = 2**13  # Overlap chunks by this amount

#os.makedirs(f'/users/nfairfie/scratch/results_ft/{fil_prefix}', exist_ok=True)
os.makedirs(f'/users/nfairfie/scratch/results/{fil_prefix}', exist_ok=True)

for i_s in range(0, file_shape[0], N_s - overlap_s):
  #result_filename = f'/users/nfairfie/scratch/results_ft/{fil_prefix}/{fil_filename}_{i_s:010}.png'
  result_filename = f'/users/nfairfie/scratch/results/{fil_prefix}/{fil_filename}_{i_s:010}.png'
  if os.path.exists(result_filename): continue  # Note: race condition with other workers right here.
  open(result_filename, 'w')  # write a placeholder to claim this chunk

  print(f'processing chunk starting at sample {i_s}...')
  obs = bl.Waterfall(fil_fullname, t_start=i_s, t_stop=(i_s + N_s), max_load=2.0)
  D = np.squeeze(obs.data).T
  D = D[f_start:(f_end + 1)]
  print(D.shape)
  # Normalize (-min/max) per channel
  D = D - np.min(D, axis=1)[:, None]
  D = D / np.max(D, axis=1)[:, None]
  D = np.where(np.isnan(D), 0, D)  # Zero out NaNs
  D[-5:, :] = 0  # For GB-20 256 channel data, these channels are supposed to
                 # be 0 anyway, due to the bandpass... but sometimes it
                 # *really* isn't!
  # TODO: fix DMT flipping
  DMT = FDMT(D[::-1], f_min, f_max, ds_max, 'float32')  # Compute the DMT
  print(DMT.shape)
  DMT = DMT[ds_min:, ds_max:]  # Crop off low DMs, and the first ds_max samples. TODO: fix DMT flipping
  # Normalize (-mean/std) per DM row
  DMT = DMT - np.nanmean(DMT, axis=1)[:, None]
  DMT = DMT / np.nanstd(DMT, axis=1)[:, None]

  DMs = np.linspace(DM_min, DM_max, num=DMT.shape[0])
  DM_index = np.abs(DMs - 670).argmin()

#  fft = np.fft.fft(DMT[DM_index, :])
#  fftfreq = np.fft.fftfreq(DMT.shape[1], dt)
#  FFT_f_min = 10 # Hz
#  FFT_f_max = 100 # Hz
#  FFT_start = np.argmin(np.abs(fftfreq - FFT_f_min))
#  FFT_end = np.argmin(np.abs(fftfreq - FFT_f_max))
#  fft_max_i = FFT_start + np.abs(fft[FFT_start:FFT_end]).argmax()
#  fft_max = fft[fft_max_i]
#  fftfreq_max = fftfreq[fft_max_i]
#  if fft_max > 1000:

  DMT_max = DMT.max()
  if (DMT_max > 6.0):
    fig = pl.figure(figsize=(16, 12))
    #pl.subplot(2, 1, 1)
    pl.imshow(DMT, cmap='magma', origin='lower', aspect='auto',
        extent=((i_s + ds_max) * dt, (i_s + ds_max + DMT.shape[1]) * dt, DM_min, DM_max))
    pl.ylabel('DM (pc/cm^3')
    pl.xlabel('time (s)')
    #pl.grid()
    #pl.subplot(2, 1, 2)
    #pl.margins(x=0)  # To more closely align the plot with the image
    #pl.plot(i_s + ds_min + np.arange(DMT.shape[1]), np.sum(DMT, axis=0))
    pl.title(f'{fil_filename}_{i_s:010} DMT_max={DMT_max}')
    #pl.title(f'{fil_filename}_{i_s:010} fftfreq_max={fftfreq_max} fft_max={fft_max}')
    pl.savefig(result_filename)
    pl.close(fig)
  else:
    open(result_filename, 'w').write(f'DMT_max={DMT_max}')
#    open(result_filename, 'w').write(f'fftfreq_max={fftfreq_max} fft_max={fft_max}')
