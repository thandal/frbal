"""Per-channel preprocessing for the FRB search.

Separated from fdmt_search.py so the normalization can be unit-tested against
the same code that runs in production (see fdmt_validate.py).
"""
import numpy as np


def normalize_robust(D, clip_sigma=None):
    """Per-channel robust normalization: subtract the median, divide by a
    MAD-based noise estimate (1.4826 * MAD ~= Gaussian sigma).

    Why this and not min-max:
      * The median is a robust baseline (immune to the bandpass DC level and to
        a bright burst), unlike `min`.
      * Dividing by the per-channel noise sigma *whitens* the band: every
        channel then carries unit-variance noise, so the FDMT's sum across
        channels becomes the matched-filter-optimal (max-S/N) combination.
      * MAD is set by the bulk of the samples, not the extremes, so an RFI
        spike or a bright pulse does not blow up the scale and compress the
        real signal -- the failure mode of dividing by `max`.

    Dead/constant channels (sigma == 0) are zeroed. Optionally clip to
    +/-clip_sigma afterwards to limit the impact of impulsive RFI.

    Input/Output: D is [N_f, N_s]; returned in the same shape and orientation.
    """
    D = D.astype('float32', copy=True)
    med = np.median(D, axis=1, keepdims=True)
    mad = np.median(np.abs(D - med), axis=1, keepdims=True)
    sigma = (1.4826 * mad).astype('float32')
    D -= med
    with np.errstate(divide='ignore', invalid='ignore'):
        D /= sigma
    D[~np.isfinite(D)] = 0.0          # constant/dead channels (sigma == 0)
    if clip_sigma is not None:
        np.clip(D, -clip_sigma, clip_sigma, out=D)
    return D


def normalize_minmax(D):
    """Legacy per-channel min-max normalization (kept for comparison/regression).

    This is what fdmt_search.py used before normalize_robust. It is fragile:
    `max` is set by the brightest sample, so a single RFI spike -- or the very
    burst being searched for -- compresses the channel and costs S/N.
    """
    D = D.astype('float32', copy=True)
    D = D - np.min(D, axis=1)[:, None]
    D = D / np.max(D, axis=1)[:, None]
    D = np.where(np.isnan(D), 0, D)
    return D
