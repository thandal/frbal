"""Boxcar matched-filter width search over an FDMT dispersion-measure transform.

The FDMT gives one dedispersed sample per (DM, time) cell. Real FRBs span several
samples (intrinsic width, scattering, intra-channel smearing at high DM), so a
single-sample statistic is mismatched and loses ~sqrt(W) S/N for a width-W pulse.
This module convolves each DM row with a bank of boxcars and keeps, per cell, the
best (max) significance across widths -- a matched-filter bank over pulse width.

Separated from fdmt_search.py so it can be unit-tested (see fdmt_validate.py).
"""
import numpy as np


def default_widths(max_width):
    """Powers-of-two boxcar widths 1, 2, 4, ... up to (and including) max_width."""
    widths, w = [], 1
    while w <= max_width:
        widths.append(w)
        w *= 2
    return widths


def boxcar_search(dmt, widths=None, max_width=64):
    """Matched-filter width search over a dispersion-measure transform.

    For each trial width W, sum each DM row over a sliding window of W samples,
    then robustly z-score (subtract median, divide by 1.4826*MAD) so the response
    is in units of sigma and is directly comparable across widths. A pulse of
    intrinsic width ~W is best matched by the W-boxcar.

    Robust (median/MAD) normalization also means the burst being searched for, or
    residual RFI, does not bias its own detection threshold -- unlike mean/std.

    Args:
      dmt: [n_dm, n_t] dedispersed intensity (cropped FDMT output). Operated on
           as-is; per-row baseline is removed internally.
      widths: explicit iterable of boxcar widths (samples); default powers of 2.
      max_width: cap for auto-generated widths; also clamped to n_t // 4.

    Returns:
      detect: [n_dm, n_t] float32, the per-cell max significance across widths,
        time-aligned to the window CENTER so it shares `dmt`'s extent.
      best: dict {snr, i_dm, i_t, width} for the global peak.
    """
    n_dm, n_t = dmt.shape
    if widths is None:
        widths = default_widths(min(max_width, max(1, n_t // 4)))

    cs = np.zeros((n_dm, n_t + 1), dtype='float64')
    np.cumsum(dmt, axis=1, out=cs[:, 1:])

    detect = np.full((n_dm, n_t), -np.inf, dtype='float64')
    best = {'snr': -np.inf, 'i_dm': 0, 'i_t': 0, 'width': int(widths[0])}
    for W in widths:
        if W > n_t:
            break
        S = cs[:, W:] - cs[:, :-W]                 # S[:,k] = sum over samples [k, k+W-1]
        med = np.median(S, axis=1, keepdims=True)
        mad = np.median(np.abs(S - med), axis=1, keepdims=True)
        sigma = 1.4826 * mad
        sigma[sigma == 0] = np.inf                 # dead/constant row -> z = 0
        z = (S - med) / sigma
        off = (W - 1) // 2                          # center the window in time
        view = detect[:, off:off + z.shape[1]]
        np.maximum(view, z, out=view)
        r, k = np.unravel_index(int(np.argmax(z)), z.shape)
        if z[r, k] > best['snr']:
            best = {'snr': float(z[r, k]), 'i_dm': int(r), 'i_t': int(k + off), 'width': int(W)}

    detect[~np.isfinite(detect)] = 0.0              # untouched edges (shouldn't occur after W=1)
    return detect.astype('float32'), best
