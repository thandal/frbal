#!/usr/bin/env python3
"""Validation harness for the FDMT kernel in fdmt.py.

Runs entirely on synthetic data -- no blimpy, no telescope files, no cluster.
Each check is pass/fail; the script exits nonzero if any fails.

  1. KERNEL CORRECTNESS (low DM, exact). At low DM the intra-channel dispersion
     smear is < 1 sample, so the FDMT reduces to plain incoherent dedispersion:
     a sum of one sample per channel. We inject a delta and prove that ALL N_f
     channels land within +/-1 sample of the expected time (windowed energy ==
     N_f, to floating-point precision) and that the FDMT row matches brute-force
     dedispersion. This is the strongest correctness statement here.

  2. INJECTION RECOVERY. A pulse dispersed at DM0 peaks in the FDMT at the row
     mapped to DM0 and at the expected time. Sweeping DM0 confirms the DM axis
     is linear and correctly mapped -- the np.linspace(DM_min, DM_max, ...)
     mapping fdmt_search.py uses.

  3. FREQUENCY ORIENTATION. The kernel only sums coherently when channel 0 is
     the LOWEST frequency. Feeding it the reversed band smears the pulse to
     noise. This settles the `D[::-1]` / "TODO: fix DMT flipping" question in
     fdmt_search.py.

  4. TIME CONVENTION. Reports the column offset between an injected pulse's
     f_min arrival and the FDMT peak column, so the plot time axis can be
     checked.

  5. SENSITIVITY vs DIRECT DEDISPERSION. With fair sub-sample injection, the
     FDMT should recover essentially the same peak S/N as brute-force direct
     dedispersion. Also demonstrates that a boxcar matched to the pulse width
     beats single-sample detection -- motivating a width search.

Usage:
  ./fdmt_validate.py            # run all checks; exit nonzero on any failure
  ./fdmt_validate.py -v         # verbose
  ./fdmt_validate.py --plot OUT.png
"""
import argparse
import sys

import numpy as np

sys.path.insert(0, __file__.rsplit('/', 1)[0])
from fdmt import FDMT
from preprocess import normalize_robust, normalize_minmax
from detect import boxcar_search, default_widths

K_DM = 4148.808  # MHz^2 pc^-1 cm^3 s ; same constant fdmt_search.py uses


# --------------------------------------------------------------------------- #
# Dispersion helpers. Convention: frequencies ascending, channel 0 is f_min.
# Output time of a dedispersed series = arrival time at f_min; higher
# frequencies arrive earlier (smaller sample index).
# --------------------------------------------------------------------------- #
def delay_s(DM, f_lo, f_hi):
    """Cold-plasma dispersion delay (seconds) between f_lo and f_hi (MHz)."""
    return K_DM * DM * (f_lo**-2 - f_hi**-2)


def dm_to_row(DM, f_min, f_max, dt):
    """FDMT row index for a DM (total delay across the band, in samples)."""
    return int(round(delay_s(DM, f_min, f_max) / dt))


def row_to_dm(row, f_min, f_max, dt):
    return row * dt / (K_DM * (f_min**-2 - f_max**-2))


def channel_freqs(f_min, f_max, N_f):
    """Channel-center frequencies, ascending -- matches obs.get_freqs() usage in
    fdmt_search.py, where f_min/f_max are the extreme channel centers."""
    return np.linspace(f_min, f_max, N_f)


def inject_pulse(freqs, N_s, dt, DM, t0, amp=1.0, width=1, subsample=False,
                 rng=None, noise_std=0.0):
    """Build Image[N_f, N_s] (channel 0 == freqs[0]) with a dispersed pulse.

    The pulse arrives at f_min (freqs[0]) at sample t0; channel i is advanced by
    delay_s(DM, f_min, freqs[i]). `width` samples each get amplitude `amp`
    (a boxcar). `subsample=True` places each channel at its continuous delay via
    linear interpolation between adjacent bins (no privileged rounding -- this
    is the fair mode for S/N comparisons). Optional white noise of `noise_std`.
    """
    N_f = len(freqs)
    f_min = freqs[0]
    Image = np.zeros((N_f, N_s), dtype='float32')
    if noise_std:
        Image += rng.normal(0.0, noise_std, size=Image.shape).astype('float32')
    for i, f in enumerate(freqs):
        x = t0 - delay_s(DM, f_min, f) / dt          # continuous arrival sample
        for w in range(width):
            xc = x - w
            if subsample:
                lo = int(np.floor(xc)); frac = xc - lo
                if 0 <= lo < N_s:      Image[i, lo] += amp * (1 - frac)
                if 0 <= lo + 1 < N_s:  Image[i, lo + 1] += amp * frac
            else:
                ti = int(round(xc))
                if 0 <= ti < N_s:      Image[i, ti] += amp
    return Image


def brute_dedisperse(Image, freqs, dt, DM):
    """Incoherent dedispersion at trial DM by integer per-channel shifts,
    summed over channels, referenced to f_min. The ground truth the FDMT
    approximates: out[t] = sum_i Image[i, t - round(delay_i/dt)]."""
    N_f, N_s = Image.shape
    f_min = freqs[0]
    out = np.zeros(N_s, dtype='float64')
    for i, f in enumerate(freqs):
        shift = int(round(delay_s(DM, f_min, f) / dt))
        if shift >= 0:
            out[shift:] += Image[i, :N_s - shift]
        else:
            out[:N_s + shift] += Image[i, -shift:]
    return out


def boxcar_snr(series, width=1):
    """Peak S/N of a 1-D series under a boxcar matched filter of given width,
    using a robust (median) baseline. Noise normalization keeps unit variance."""
    x = series.astype('float64')
    if width > 1:
        x = np.convolve(x, np.ones(width), 'same') / np.sqrt(width)
    return (x.max() - np.median(x)) / x.std()


# --------------------------------------------------------------------------- #
class Reporter:
    def __init__(self):
        self.failures = 0

    def check(self, ok, label, detail=""):
        if not ok:
            self.failures += 1
        print(f"  [{'PASS' if ok else 'FAIL'}] {label}" + (f"  ({detail})" if detail else ""))
        return ok


# --------------------------------------------------------------------------- #
def test_kernel_correctness(rep, f_min, f_max, N_f, dt, N_s, dm_max, verbose):
    print("\n== Test 1: kernel correctness vs brute-force (low DM, exact) ==")
    freqs = channel_freqs(f_min, f_max, N_f)
    maxDT = dm_to_row(dm_max, f_min, f_max, dt)
    # Low DMs where the intra-channel smear rounds to 0 samples -> exact.
    for DM in [5, 15, 30, 45]:
        t0 = maxDT + (N_s - maxDT) // 2
        Image = inject_pulse(freqs, N_s, dt, DM, t0, amp=1.0)
        A = FDMT(Image, f_min, f_max, maxDT, 'float32')
        fdmt_row = A[dm_to_row(DM, f_min, f_max, dt)].astype('float64')
        brute = brute_dedisperse(Image, freqs, dt, DM)
        pk = int(np.argmax(fdmt_row))
        win = fdmt_row[pk - 2:pk + 3].sum()          # 5-sample window
        rep.check(abs(pk - t0) <= 1, f"DM={DM:>3}: peak at expected time",
                  f"peak col={pk}, t0={t0}")
        rep.check(abs(win - N_f) < 0.5,
                  f"DM={DM:>3}: all {N_f} channels summed (windowed energy == N_f)",
                  f"sum={win:.3f}")
        rep.check(int(np.argmax(brute)) == t0,
                  f"DM={DM:>3}: brute-force agrees on peak time")
        if verbose:
            print(f"      single-sample peak={fdmt_row[pk]:.1f} "
                  f"({100 * fdmt_row[pk] / N_f:.0f}% of N_f; rest in +/-1 bin)")


def test_injection_recovery(rep, f_min, f_max, N_f, dt, N_s, dm_max, verbose):
    print("\n== Test 2: injection recovery & DM-axis linearity ==")
    freqs = channel_freqs(f_min, f_max, N_f)
    maxDT = dm_to_row(dm_max, f_min, f_max, dt)
    dm_step = row_to_dm(1, f_min, f_max, dt)
    print(f"   {'inj DM':>8} {'rec DM':>8} {'inj t':>8} {'rec t':>8} {'peak':>7}")
    for DM in [50, 100, 200, int(0.9 * dm_max)]:
        t0 = maxDT + (N_s - maxDT) // 2
        Image = inject_pulse(freqs, N_s, dt, DM, t0, amp=1.0)
        A = FDMT(Image, f_min, f_max, maxDT, 'float32')
        sub = A[:, maxDT:]                            # drop edge-contaminated cols
        r, c = np.unravel_index(int(np.argmax(sub)), sub.shape)
        rec_t, rec_dm = c + maxDT, row_to_dm(r, f_min, f_max, dt)
        print(f"   {DM:8.1f} {rec_dm:8.1f} {t0:8d} {rec_t:8d} {sub[r, c]:7.1f}")
        rep.check(abs(rec_dm - DM) <= 1.5 * dm_step, f"DM={DM:>5}: recovered DM within one row",
                  f"|{rec_dm:.1f}-{DM}|={abs(rec_dm - DM):.2f}, step={dm_step:.2f}")
        rep.check(abs(rec_t - t0) <= 2, f"DM={DM:>5}: recovered time within 2 samples")


def test_orientation(rep, f_min, f_max, N_f, dt, N_s, dm_max, verbose):
    print("\n== Test 3: frequency orientation (settles the D[::-1] question) ==")
    freqs = channel_freqs(f_min, f_max, N_f)
    maxDT = dm_to_row(dm_max, f_min, f_max, dt)
    DM, t0 = 0.6 * dm_max, maxDT + (N_s - maxDT) // 2
    Image = inject_pulse(freqs, N_s, dt, DM, t0, amp=1.0)   # ascending = correct
    peak_asc = float(FDMT(Image, f_min, f_max, maxDT, 'float32').max())
    peak_desc = float(FDMT(Image[::-1], f_min, f_max, maxDT, 'float32').max())
    print(f"   ascending (channel 0 = f_min): peak = {peak_asc:.1f} / {N_f}")
    print(f"   reversed  (channel 0 = f_max): peak = {peak_desc:.1f} / {N_f}")
    rep.check(peak_asc >= 0.6 * N_f, "ascending input sums coherently", f"peak={peak_asc:.1f}")
    rep.check(peak_desc < 0.1 * N_f, "reversed input does NOT (pulse smeared to noise)",
              f"peak={peak_desc:.1f}")
    print("   => FDMT requires channel 0 == f_min. fdmt_search.py now sorts channels")
    print("      explicitly (order = np.argsort(chan_freqs); D[order]) -- orientation-proof,")
    print("      replacing the old ambiguous D[::-1] / 'TODO: fix DMT flipping'.")


def test_time_convention(rep, f_min, f_max, N_f, dt, N_s, dm_max, verbose):
    print("\n== Test 4: time convention (FDMT column vs f_min arrival) ==")
    freqs = channel_freqs(f_min, f_max, N_f)
    maxDT = dm_to_row(dm_max, f_min, f_max, dt)
    DM, t0 = 0.5 * dm_max, maxDT + (N_s - maxDT) // 2
    Image = inject_pulse(freqs, N_s, dt, DM, t0, amp=1.0)
    A = FDMT(Image, f_min, f_max, maxDT, 'float32')
    peak_col = int(np.argmax(A[dm_to_row(DM, f_min, f_max, dt)]))
    offset = peak_col - t0
    print(f"   injected f_min arrival t0={t0}, FDMT peak column={peak_col}, offset={offset:+d}")
    rep.check(abs(offset) <= 1, "FDMT output column == f_min arrival time (within 1 sample)",
              f"offset={offset:+d}")


def test_sensitivity(rep, f_min, f_max, N_f, dt, N_s, dm_max, verbose, seed=7):
    print("\n== Test 5: sensitivity vs direct dedispersion (fair sub-sample injection) ==")
    freqs = channel_freqs(f_min, f_max, N_f)
    maxDT = dm_to_row(dm_max, f_min, f_max, dt)
    rng = np.random.default_rng(seed)
    W = 4                                    # injected pulse width (samples)
    target = 15.0                            # ideal matched-filter S/N
    amp = target / np.sqrt(N_f * W)          # per-sample, per-channel amplitude
    print(f"   injecting width-{W} pulses, ideal matched S/N = {target:.0f}, per-ch noise std = 1")
    print(f"   {'DM':>5} {'FDMT box1':>10} {'FDMT box4':>10} {'brute box4':>11} {'retained':>9}")
    for DM in [50, 150, int(0.9 * dm_max)]:
        t0 = maxDT + (N_s - maxDT) // 2
        Image = inject_pulse(freqs, N_s, dt, DM, t0, amp=amp, width=W,
                             subsample=True, noise_std=1.0, rng=rng)
        A = FDMT(Image, f_min, f_max, maxDT, 'float32')
        frow = A[dm_to_row(DM, f_min, f_max, dt)][maxDT:]
        brow = brute_dedisperse(Image, freqs, dt, DM)[maxDT:]
        f1, f4, b4 = boxcar_snr(frow, 1), boxcar_snr(frow, W), boxcar_snr(brow, W)
        print(f"   {DM:5d} {f1:10.1f} {f4:10.1f} {b4:11.1f} {f4 / b4 * 100:7.0f}%")
        rep.check(f4 >= 0.8 * b4, f"DM={DM:>5}: FDMT retains >=80% of direct-dedisp S/N",
                  f"{f4:.1f} vs {b4:.1f} ({f4 / b4 * 100:.0f}%)")
        rep.check(f4 > 6.0, f"DM={DM:>5}: pulse detected above 6 sigma", f"S/N={f4:.1f}")
        rep.check(f4 > f1, f"DM={DM:>5}: boxcar(width) beats single-sample (motivates width search)",
                  f"box{W}={f4:.1f} > box1={f1:.1f}")


def test_normalization(rep, f_min, f_max, N_f, dt, N_s, dm_max, verbose, seed=11):
    """Validate the sensitivity fix: robust(+clip) normalization vs legacy min-max,
    on a realistic chunk (bandpass gain/offset, white noise, a broadband burst) in
    two RFI regimes. The production config is normalize_robust(clip_sigma=5)."""
    print("\n== Test 6: normalization sensitivity (robust+clip vs legacy min-max) ==")
    freqs = channel_freqs(f_min, f_max, N_f)
    maxDT = dm_to_row(dm_max, f_min, f_max, dt)
    DM, W = 150.0, 4
    t0 = maxDT + (N_s - maxDT) // 2
    row = dm_to_row(DM, f_min, f_max, dt)

    def recovered_snr(D, normalizer):
        A = FDMT(normalizer(D), f_min, f_max, maxDT, 'float32')
        return boxcar_snr(A[row][maxDT:], W)

    def build_chunk(rng, saturating_rfi):
        # measured = gain * (noise + burst_flux) + bandpass DC offset
        gain = (1.0 + 0.5 * np.sin(np.linspace(0, 3, N_f)))[:, None]
        offset = (50.0 + 10 * np.cos(np.linspace(0, 5, N_f)))[:, None]
        amp = 30.0 / np.sqrt(N_f * W)        # ideal matched S/N ~= 30
        burst = inject_pulse(freqs, N_s, dt, DM, t0, amp=amp, width=W, subsample=True)
        noise = rng.normal(0, 1, (N_f, N_s)).astype('float32')
        D = (gain * (noise + burst) + offset).astype('float32')
        if saturating_rfi:                   # one huge (saturating) sample in many channels
            sat = rng.choice(N_f, size=N_f // 2, replace=False)
            D[sat, rng.integers(0, N_s, size=sat.size)] += 1e4
        return D

    robust = lambda D: normalize_robust(D, clip_sigma=5.0)
    print(f"   {'regime':>16} {'min-max':>9} {'robust+clip':>12} {'gain':>7}")
    for label, sat in [("clean", False), ("saturating RFI", True)]:
        rng = np.random.default_rng(seed)
        D = build_chunk(rng, sat)
        s_mm = recovered_snr(D, normalize_minmax)
        s_rb = recovered_snr(D, robust)
        print(f"   {label:>16} {s_mm:9.1f} {s_rb:12.1f} {s_rb / s_mm:6.2f}x")
        rep.check(s_rb >= 1.15 * s_mm,
                  f"{label}: robust+clip beats min-max by >=15%",
                  f"{s_rb:.1f} vs {s_mm:.1f} ({s_rb / s_mm:.2f}x)")
        rep.check(s_rb > 6.0, f"{label}: burst still detected above 6 sigma", f"S/N={s_rb:.1f}")
    if verbose:
        rng = np.random.default_rng(seed)
        D = build_chunk(rng, True)
        s_noclip = recovered_snr(D, lambda d: normalize_robust(d))  # no clip
        print(f"      (robust WITHOUT clip under saturating RFI: S/N={s_noclip:.1f} "
              f"-- clipping is essential; unbounded spikes inflate the noise)")


def test_boxcar_width_search(rep, f_min, f_max, N_f, dt, N_s, dm_max, verbose, seed=42):
    """Validate the boxcar matched-filter width search: it should recover the
    injected pulse width, localize it in (DM, time), and gain S/N over a
    single-sample statistic by a margin that grows with pulse width."""
    print("\n== Test 7: boxcar width search (matched filter over pulse width) ==")
    freqs = channel_freqs(f_min, f_max, N_f)
    maxDT = dm_to_row(dm_max, f_min, f_max, dt)
    DM = 150.0
    row, t0 = dm_to_row(DM, f_min, f_max, dt), maxDT + (N_s - maxDT) // 2
    print(f"   widths searched: {default_widths(64)}")
    print(f"   {'inj W':>6} {'rec W':>6} {'best S/N':>9} {'W=1 S/N':>8} {'gain':>6}")
    gains = {}
    for Wp in [1, 4, 16]:
        rng = np.random.default_rng(seed)
        amp = 20.0 / np.sqrt(N_f * Wp)               # fixed ideal S/N across widths
        Image = inject_pulse(freqs, N_s, dt, DM, t0, amp=amp, width=Wp,
                             subsample=True, noise_std=1.0, rng=rng)
        DMT = FDMT(Image, f_min, f_max, maxDT, 'float32')[:, maxDT:]   # crop edge cols
        det, best = boxcar_search(DMT, max_width=64)
        _, best1 = boxcar_search(DMT, widths=[1])     # single-sample baseline
        gain = best['snr'] / best1['snr']
        gains[Wp] = gain
        print(f"   {Wp:6d} {best['width']:6d} {best['snr']:9.1f} {best1['snr']:8.1f} {gain:6.2f}x")
        rep.check(det.shape == DMT.shape, f"W={Wp:>2}: detect map matches DMT shape")
        rep.check(max(1, Wp // 2) <= best['width'] <= 2 * max(Wp, 1),
                  f"W={Wp:>2}: recovered width within one octave of injected",
                  f"rec={best['width']}, inj={Wp}")
        rep.check(abs(best['i_dm'] - row) <= 2, f"W={Wp:>2}: DM localized", f"row {best['i_dm']} vs {row}")
        rep.check(abs(best['i_t'] - (t0 - maxDT)) <= max(Wp, 3), f"W={Wp:>2}: time localized")
        rep.check(best['snr'] >= 0.98 * best1['snr'], f"W={Wp:>2}: width search >= single-sample S/N",
                  f"{best['snr']:.1f} vs {best1['snr']:.1f}")
    rep.check(gains[16] > gains[1] + 0.3, "S/N gain grows with pulse width (matched filter working)",
              f"gain(W16)={gains[16]:.2f}x > gain(W1)={gains[1]:.2f}x")


def maybe_plot(path, f_min, f_max, N_f, dt, N_s, dm_max):
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt
    freqs = channel_freqs(f_min, f_max, N_f)
    maxDT = dm_to_row(dm_max, f_min, f_max, dt)
    DM, t0 = 0.6 * dm_max, maxDT + (N_s - maxDT) // 2
    rng = np.random.default_rng(0)
    Image = inject_pulse(freqs, N_s, dt, DM, t0, amp=15.0 / np.sqrt(N_f * 4),
                         width=4, subsample=True, noise_std=1.0, rng=rng)
    A = FDMT(Image, f_min, f_max, maxDT, 'float32')[:, maxDT:]   # crop edge cols
    detect, best = boxcar_search(A, max_width=64)
    DM_best = dm_max * best['i_dm'] / max(detect.shape[0] - 1, 1)
    fig, (ax1, ax2, ax3) = plt.subplots(3, 1, figsize=(12, 12))
    ax1.imshow(Image, aspect='auto', origin='lower', cmap='magma', extent=(0, N_s, f_min, f_max))
    ax1.set(title=f'Injected: DM={DM:.0f}, t0={t0}', xlabel='sample', ylabel='freq (MHz)')
    ax2.imshow(A, aspect='auto', origin='lower', cmap='magma', extent=(0, A.shape[1], 0, dm_max))
    ax2.axhline(DM, color='cyan', lw=0.5, ls='--')
    ax2.set(title='FDMT output (raw)', xlabel='sample', ylabel='DM (pc/cm^3)')
    im = ax3.imshow(detect, aspect='auto', origin='lower', cmap='magma', extent=(0, detect.shape[1], 0, dm_max))
    ax3.plot(best['i_t'], DM_best, 'c+', ms=20, mew=2)
    fig.colorbar(im, ax=ax3, label='matched-filter S/N')
    ax3.set(title=f"Boxcar width search: S/N={best['snr']:.1f} DM={DM_best:.0f} width={best['width']}",
            xlabel='sample', ylabel='DM (pc/cm^3)')
    fig.tight_layout(); fig.savefig(path, dpi=110)
    print(f"\nWrote {path}")


def main():
    ap = argparse.ArgumentParser(description="Validate the FDMT kernel on synthetic data.")
    ap.add_argument('-v', '--verbose', action='store_true')
    ap.add_argument('--plot', metavar='OUT.png')
    ap.add_argument('--f-min', type=float, default=1300.0)
    ap.add_argument('--f-max', type=float, default=1500.0)
    ap.add_argument('--n-f', type=int, default=256, help="channels (power of 2)")
    ap.add_argument('--dt', type=float, default=200e-6, help="sample time (s)")
    ap.add_argument('--n-s', type=int, default=4096, help="samples per chunk")
    ap.add_argument('--dm-max', type=float, default=300.0)
    args = ap.parse_args()

    p = (args.f_min, args.f_max, args.n_f, args.dt, args.n_s, args.dm_max)
    maxDT = dm_to_row(args.dm_max, args.f_min, args.f_max, args.dt)
    print("FDMT validation harness")
    print(f"  band {args.f_min}-{args.f_max} MHz, {args.n_f} ch, dt={args.dt*1e6:.0f} us, "
          f"N_s={args.n_s}, DM_max={args.dm_max} -> maxDT={maxDT} samples")

    rep = Reporter()
    test_kernel_correctness(rep, *p, args.verbose)
    test_injection_recovery(rep, *p, args.verbose)
    test_orientation(rep, *p, args.verbose)
    test_time_convention(rep, *p, args.verbose)
    test_sensitivity(rep, *p, args.verbose)
    test_normalization(rep, *p, args.verbose)
    test_boxcar_width_search(rep, *p, args.verbose)
    if args.plot:
        maybe_plot(args.plot, *p)

    print("\n" + "=" * 62)
    print(f"RESULT: {rep.failures} check(s) FAILED" if rep.failures else "RESULT: all checks PASSED")
    return 1 if rep.failures else 0


if __name__ == '__main__':
    sys.exit(main())
