"""
period_id_exp.py -- identify the model's fundamental OOP period the way the
                    EXPERIMENT does, not the way the driver currently does.

WHY THIS EXISTS
---------------
Moshfeghi, Smyrou, Arslan & Bal (2024), Structures 66:106815, section 3.2,
states the measurement protocol verbatim:

  "For each test run, a silent period of post-shaking recording was used, as
   long as 45 to 60s, for collecting clear free vibration signals from the
   wall. ... it is important to find the exact moment the shake table actually
   stops. The displacement sensor attached to the shake table is used for
   finding that exact moment and the part of the acceleration signal from that
   point on is cropped for further analysis."

  "Data from all three accelerometers, placed at the bottom quarter,
   mid-height and top-quarter positions, are used for estimating the
   periodogram of the signals, calculating the Power Spectral Densities (PSD),
   and Singular Values ... Peak-picking is applied for finding the fundamental
   periods. ... the PSD results are compared to the simple Fast Fourier
   Transform (FFT) outcomes and identical results were obtained."

The driver's `identify_Tend_from_csv` does none of that. It reads ONE channel,
`Channel_19_DispTopQRight`, which is (a) a DISPLACEMENT, not an acceleration,
and (b) despite its name, recorded at x = bLength/2 = 0.646 m, i.e. the
CENTRELINE, not the right quarter point.

The differences are not cosmetic:

  1. DISPLACEMENT vs ACCELERATION.  A displacement spectrum is weighted 1/w^2
     and an acceleration spectrum by w^2. A rocking run leaves a slow settling
     transient and a permanent offset in the displacement record. Removing the
     mean (which the driver does) kills a constant but not a settling ramp,
     and that ramp carries enormous power at low frequency. The displacement
     periodogram peak is then pulled toward it. This is the most likely origin
     of run 12's T_end = 0.21848 s -- 2.09x the previous run in a single step,
     with the autocorrelation and FFT estimators AGREEING to within tolerance,
     because both were looking at the same low-frequency ramp. Double
     differentiation annihilates a ramp exactly.

  2. ONE CHANNEL vs THREE.  A single centreline channel can lock onto a local
     mode of whatever course it sits on. The singular-value decomposition of
     the three-channel PSD matrix is the standard remedy (this is Frequency
     Domain Decomposition); it is what the paper did.

  3. WINDOW LENGTH.  This one CANNOT be matched and the gap must be stated in
     the paper, not buried. The experiment used 45-60 s of silence. The model
     has tail_sec (2.5 s) + inter_run_gap (0.5 s) = 3.0 s, of which the driver
     currently uses at most MAX_RD_WINDOW = 1.5 s. Over 45 s the free-vibration
     amplitude decays to almost nothing, so the experimental estimate is
     dominated by the LOW-AMPLITUDE tail -- and the paper's own abstract notes
     that "even serious cracks caused by OOP response close when the shaking
     stops". The experimental period is therefore a CLOSED-CRACK, small-
     amplitude property. A 1.5 s model window is dominated by the first ~15
     cycles, while cracks are still swinging open, which measures a different
     and much longer period. Matching the estimator removes a large part of the
     discrepancy but not this part.

WHAT THIS MODULE DOES NOT DO
----------------------------
It does not change the experiment. It changes how the MODEL is measured so the
two are the same quantity.

CHANNELS USED (already recorded by instrument_history_new.dat, lines 159-164 --
no new instrumentation and no model rerun is required):

    Channel_12_AccTable   velocity-z at (0.646, 0.00, 0.21)   table reference
    Channel_15_AccBot     velocity-z at (0.646, 0.66, 0.21)   bottom quarter
    Channel_16_AccMid     velocity-z at (0.646, 1.26, 0.21)   mid-height
    Channel_17_AccTop     velocity-z at (0.646, 2.06, 0.21)   top quarter

Note these are VELOCITY histories despite the "Acc" names; the source file says
so in its own comment ("VELOCITY CHANNELS (for acceleration derivation)"). They
are differentiated once here. The sensor elevations 0.66 / 1.26 / 2.06 m match
the paper's bottom-quarter / mid-height / top-quarter accelerometer positions.
"""

import os
import math
import numpy as np

try:
    from scipy.signal import csd as _scipy_csd
    _HAVE_SCIPY = True
except Exception:                                    # pragma: no cover
    _HAVE_SCIPY = False


# ---------------------------------------------------------------------
# Configuration
# ---------------------------------------------------------------------
CH_TABLE = "Channel_12_AccTable"
CH_WALL = ("Channel_15_AccBot", "Channel_16_AccMid", "Channel_17_AccTop")

T_MIN_PHYS = 0.02     # s   physical band for peak-picking
T_MAX_PHYS = 1.00     # s

ZPAD = 16             # zero-padding factor for peak localisation
                      # (does not add resolution, only sub-bin peak location --
                      #  the parabolic refine below does the same job as the
                      #  driver's existing 16x pad, kept identical on purpose)

TABLE_STOP_TOL = 5e-2  # Fraction of the loudest chunk RMS below which the
                       # table counts as stopped.
                       #
                       # Was 1e-3, which never fired: on NODAMP_v7_NEW run 01
                       # the channel decays from 3.002e-03 to 4.848e-05 m/s,
                       # i.e. to 1.6% of peak, and sits there. 0.1% was an
                       # arbitrary number with no basis in the data.
                       #
                       # It does not decay further because Channel_12 is NOT
                       # ON THE TABLE. instrument_history_new.dat line 159 puts
                       # it at ([bLength/2], 0.0, 0.21); y = 0.0 is the
                       # wall/table interface and block.gp.near() resolves it
                       # onto a masonry gridpoint at the wall base. Group 'S'
                       # is block fix'ed and the applied velocity is removed
                       # after the pulse, so a genuine table gridpoint would
                       # read exactly zero. The 1.6% residual is the wall base
                       # still ringing.
                       #
                       # NOTE the same coordinate is used by
                       # Channel_5_DispTable (line 150), which postprocess_
                       # stratC.extract_rel_to_table subtracts to form the EDP.
                       # If that gridpoint is also on the wall, every
                       # "relative to table" displacement is relative to the
                       # wall base instead. Verify before the next write-up.
                       #
                       # This tolerance is only the FALLBACK rule. Prefer
                       # `pulse_end_s` below, which is exact.
SETTLE_SKIP_S = 0.0    # the paper crops from the stop instant itself; expose
                       # this only so the sensitivity can be reported.
MAX_WINDOW_S = None    # None = use the entire available quiet record. Set to a
                       # number only to reproduce the driver's truncation.


# ---------------------------------------------------------------------
# I/O
# ---------------------------------------------------------------------
def find_channel_file(run_dir, run_label, channel):
    """Locate one channel's exported CSV.

    Tries the exact name the FISH export builds, then falls back to scanning
    the folder for any CSV containing the channel name. The driver already
    learned this lesson once -- see _find_channel_csv in
    strategy_C_3dec_nodamp.py, whose docstring records that locating the file
    by scan "is what previously broke period ID across OS conventions". The
    first version of this module did not do that, which is why a whole run set
    can come back empty with nothing but a note to show for it.
    """
    exact = os.path.join(run_dir, "{}_{}.csv".format(run_label, channel))
    if os.path.exists(exact):
        return exact
    if not os.path.isdir(run_dir):
        return None
    needle = channel.lower()
    hits = [os.path.join(run_dir, fn) for fn in os.listdir(run_dir)
            if fn.lower().endswith(".csv") and needle in fn.lower()]
    return sorted(hits)[0] if hits else None


def load_channel(run_dir, run_label, channel):
    """Read one exported 3DEC history CSV -> (t, y). Returns (None, None)
    if the file is absent, so a caller can degrade gracefully."""
    path = find_channel_file(run_dir, run_label, channel)
    if path is None:
        return None, None
    try:
        data = np.genfromtxt(path, delimiter=',', skip_header=1)
        if data.ndim == 1 or data.shape[1] < 2:
            raise ValueError
    except (ValueError, IndexError):
        data = np.genfromtxt(path, skip_header=2)
    if data.ndim == 1 or data.shape[1] < 2:
        return None, None
    return data[:, 0], data[:, 1]


def last_segment(t, y):
    """3DEC appends successive solves to one history file; keep the last
    monotonic segment. Same rule the driver already uses."""
    d = np.diff(t)
    neg = np.where(d < -0.5)[0]
    start = neg[-1] + 1 if len(neg) else 0
    return t[start:], y[start:]


# ---------------------------------------------------------------------
# Windowing -- "find the exact moment the shake table actually stops"
# ---------------------------------------------------------------------
def table_stop_index(v_table, tol=TABLE_STOP_TOL):
    """Index of the first sample of the post-shaking quiet record.

    The experiment located this from the shake-table displacement sensor. Here
    the table velocity channel is the direct equivalent.

    WHY NOT A SIMPLE THRESHOLD ON A SINGLE SAMPLE
    -----------------------------------------------
    The first version of this took the LAST sample above tol*peak. That is
    fragile for two reasons and it silently emptied a whole 25-run set:

      - Channel_12 samples a gridpoint at y = 0.0, the wall/table interface.
        If block.gp.near() resolves to a MASONRY gridpoint rather than the
        fixed table block, the channel keeps ringing after the table stops and
        the "last sample above threshold" lands at the end of the record.
      - Even on a genuine table channel, one numerical spike anywhere in the
        tail moves the answer to that spike.

    A chunked-RMS criterion is used instead: find the FIRST chunk whose RMS
    drops below `tol` of the loudest chunk and never recovers. That is a
    sustained-quiet test, which is what "the table has stopped" actually means
    and what the experiment's displacement-sensor reading amounts to.

    Returns (index, source_string) so the caller can report which rule fired.
    """
    a = np.abs(np.asarray(v_table, dtype=float))
    n = len(a)
    if n < 64:
        return 0, "record too short"
    peak = float(np.max(a))
    if peak <= 0.0:
        return 0, "table channel is identically zero"

    n_chunk = max(8, n // 200)
    n_full = n // n_chunk
    rms = np.sqrt(np.mean(
        a[:n_full * n_chunk].reshape(n_full, n_chunk) ** 2, axis=1))
    loud = float(np.max(rms))
    if loud <= 0.0:
        return 0, "table channel is identically zero"

    quiet = rms < tol * loud
    i_stop = None
    for k in range(n_full):
        if quiet[k] and quiet[k:].all():
            i_stop = k * n_chunk
            break
    if i_stop is not None:
        return int(i_stop), "sustained-quiet RMS"

    # Never went sustainably quiet. Rather than return an unusable window,
    # fall back to the loudest-chunk rule and SAY SO, so a run identified this
    # way is visibly different from one where the table demonstrably stopped.
    k_loud = int(np.argmax(rms))
    return int(min((k_loud + 1) * n_chunk, n - 1)), \
        "FALLBACK: table never went quiet (max chunk RMS {:.3e} m/s, " \
        "final chunk {:.3e} m/s)".format(loud, float(rms[-1]))


def quiet_window(t, channels, v_table,
                 settle_skip=SETTLE_SKIP_S, max_window=MAX_WINDOW_S,
                 pulse_end_s=None):
    """Crop every channel to the post-shaking free-vibration record.

    Returns (t, channels, window_s, source). `source` is never discarded by
    the caller -- a window found by the fallback rule is reported in the log
    and the CSV, because it is a weaker claim than one found by the table
    demonstrably going quiet.
    """
    dt = float(np.median(np.diff(t))) if len(t) > 2 else 0.0
    if pulse_end_s is not None and np.isfinite(pulse_end_s) and pulse_end_s > 0:
        # EXACT rule. The driver commands the table through a velocity table
        # that returns to zero at t = n_cycles * T_excitation, both of which
        # are logged every run. This is the direct equivalent of the
        # experiment reading the stop instant off its table displacement
        # sensor, and it does not care whether Channel_12 is on the table.
        idx = np.searchsorted(np.asarray(t, dtype=float), float(pulse_end_s))
        i0, src = int(min(idx, len(t) - 1)), \
            "known pulse end t = {:.4f} s (driver log)".format(pulse_end_s)
    else:
        i0, src = table_stop_index(v_table)
    if dt <= 0:
        return None, None, 0.0, "non-monotonic time column"
    i0 += int(round(settle_skip / dt))
    i1 = len(t)
    if max_window is not None:
        i1 = min(i1, i0 + int(round(max_window / dt)))
    if i1 - i0 < 32:
        return None, None, 0.0, \
            "quiet window only {} samples (need 32); stop rule was '{}'".format(
                max(0, i1 - i0), src)
    return t[i0:i1], [c[i0:i1] for c in channels], float(t[i1 - 1] - t[i0]), src


# ---------------------------------------------------------------------
# Signal preparation -- velocity -> acceleration, detrend, taper
# ---------------------------------------------------------------------
def to_acceleration(t, v):
    """Differentiate the recorded velocity once.

    This is the step that matters. A settling ramp in displacement survives
    mean-removal; in acceleration it is annihilated. The experiment measured
    acceleration directly, so this makes the two the same quantity.
    """
    return np.gradient(np.asarray(v, dtype=float), np.asarray(t, dtype=float))


def prepare(t, a):
    """Linear detrend + Hann taper, matching a standard periodogram."""
    a = np.asarray(a, dtype=float)
    n = len(a)
    x = np.arange(n, dtype=float)
    # least-squares linear detrend (removes any residual drift the
    # differentiation did not already kill)
    A = np.vstack([x, np.ones(n)]).T
    coef, *_ = np.linalg.lstsq(A, a, rcond=None)
    a = a - A.dot(coef)
    return a * np.hanning(n)


# ---------------------------------------------------------------------
# Estimators
# ---------------------------------------------------------------------
def _parabolic_peak(freqs, amps, i):
    """Sub-bin peak refinement on log-amplitude. Identical in form to the
    refinement already used in strategy_C_3dec_nodamp.identify_Tend_from_csv,
    so any difference between the two methods is the SIGNAL and the CHANNEL
    SET, not the peak-picker."""
    if i <= 0 or i >= len(amps) - 1:
        return float(freqs[i])
    al = math.log(amps[i - 1] + 1e-300)
    be = math.log(amps[i] + 1e-300)
    ga = math.log(amps[i + 1] + 1e-300)
    den = al - 2 * be + ga
    if abs(den) < 1e-12:
        return float(freqs[i])
    df = float(freqs[1] - freqs[0])
    return float(freqs[i]) + 0.5 * (al - ga) / den * df


def _pick(freqs, spec):
    """Peak-pick inside the physical band."""
    f_lo, f_hi = 1.0 / T_MAX_PHYS, 1.0 / T_MIN_PHYS
    band = np.where((freqs >= f_lo) & (freqs <= f_hi))[0]
    if band.size == 0:
        return float('nan')
    i = int(band[np.argmax(spec[band])])
    f = _parabolic_peak(freqs, spec, i)
    return 1.0 / f if f > 0 else float('nan')


def _local_maxima(spec):
    """Indices of interior local maxima."""
    a = np.asarray(spec)
    return np.where((a[1:-1] > a[:-2]) & (a[1:-1] > a[2:]))[0] + 1


def mode_shape_at(A_list, i):
    """Complex mode shape across the three sensors at spectral bin i.

    For a single-block periodogram the cross-spectral matrix is
    G = A A^H, rank one, so the first singular VECTOR is A itself (up to a
    scalar). The relative phases of A_bot, A_mid, A_top at a candidate peak
    therefore ARE the mode shape at that frequency -- no extra machinery.
    """
    return np.array([A[i] for A in A_list])


def is_fundamental_shape(shape, max_phase_dev_deg=60.0):
    """Does this look like the FUNDAMENTAL one-way OOP bending mode?

    The specimen spans one way between the fixed base and the horizontally
    restrained top beam, so the first mode is a half sine over 0 -> 2.58 m.
    All three sensors (0.66, 1.26, 2.06 m) sit inside that single lobe and
    therefore move IN PHASE.

    Every higher mode and every impact harmonic puts at least one node between
    the sensors, which flips the relative phase of the channels either side of
    it by 180 degrees. Phase coherence across the three sensors is thus an
    exact discriminator between the fundamental and everything else, and it
    needs no arbitrary frequency band -- which matters, because a band chosen
    to exclude today's harmonic is a band that will silently exclude a
    genuinely elongated period tomorrow. That is the same failure mode as
    TID_MAX_JUMP freezing the loop for fourteen runs.

    Returns (ok, max_phase_deviation_deg).
    """
    mag = np.abs(shape)
    if not np.all(np.isfinite(mag)) or mag.max() <= 0:
        return False, 180.0
    ref = shape[int(np.argmax(mag))]
    dev = 0.0
    for c, m in zip(shape, mag):
        if m < 0.05 * mag.max():
            continue          # near a node: phase is meaningless, skip
        ang = abs(math.degrees(np.angle(c / ref)))
        dev = max(dev, min(ang, 360.0 - ang))
    return dev <= max_phase_dev_deg, dev


def pick_fundamental(freqs, spec, A_list, n_peaks=8, zpad=None,
                     max_phase_dev_deg=60.0):
    """Peak-pick the FUNDAMENTAL, not merely the tallest peak.

    Candidates are the `n_peaks` strongest local maxima in the physical band,
    forced at least one TRUE resolution bin apart. The TALLEST candidate whose
    mode shape is in phase across all three sensors wins.

    A relative power floor was tried first and abandoned. `spec` is power, so a
    5% floor is a 22% AMPLITUDE floor, and in the harmonic regression test that
    excluded the fundamental itself -- the run then returned NaN rather than
    the right answer, which is a worse failure than picking the wrong peak
    because it looks like missing data. Lowering the floor instead runs it into
    the Hann sidelobe level (-31 dB in power), so there is no safe value. Top-N
    with a resolution-bin separation needs no threshold at all: spurious peaks
    admitted this way are killed by the phase test, because uncorrelated
    content is not phase-coherent across three sensors.

    WHY TALLEST AND NOT LOWEST-FREQUENCY
    ------------------------------------
    Lowest-frequency-first was tried first, on the reasoning that an impact
    harmonic always sits ABOVE the fundamental so scanning upward could not be
    fooled by a taller harmonic. It regressed the self-test: the rigid-body
    rocking transient sits BELOW the fundamental and, being a rotation of the
    whole wall about its base, is ALSO in phase across all three sensors. The
    ascending rule therefore picked 0.250 s where the true elastic period was
    0.104 s -- reintroducing exactly the failure this module exists to remove.

    "Tallest phase-coherent" handles both, for different reasons:
      - impact harmonics may be TALLER than the fundamental but are never in
        phase, because every one of them puts a node between the sensors;
      - the rocking transient IS in phase but is not the tallest peak once the
        signal is differentiated to acceleration, because w^2 weighting
        suppresses it.
    Neither escape route is shared, so the two failure modes cannot combine.

    When every candidate is phase-coherent this degrades exactly to the naive
    tallest-peak rule, which is the desired behaviour on a clean elastic
    ring-down.

    Returns (T, info_dict).
    """
    f_lo, f_hi = 1.0 / T_MAX_PHYS, 1.0 / T_MIN_PHYS
    band = np.where((freqs >= f_lo) & (freqs <= f_hi))[0]
    if band.size == 0:
        return float('nan'), {"candidates": [], "reason": "empty band"}

    peak_val = float(np.max(spec[band]))
    df_true = float(freqs[1] - freqs[0]) * (zpad if zpad else ZPAD)

    cand_idx = [i for i in _local_maxima(spec) if f_lo <= freqs[i] <= f_hi]
    if not cand_idx:
        cand_idx = [int(band[np.argmax(spec[band])])]
    # strongest first, then greedily drop anything within one true bin of an
    # already-kept peak (zero-padding creates many maxima per real peak)
    cand_idx.sort(key=lambda i: -spec[i])
    maxima = []
    for i in cand_idx:
        if all(abs(freqs[i] - freqs[j]) >= 0.5 * df_true for j in maxima):
            maxima.append(i)
        if len(maxima) >= n_peaks:
            break
    maxima.sort(key=lambda i: freqs[i])

    cands = []
    for i in maxima:
        shape = mode_shape_at(A_list, i)
        ok, dev = is_fundamental_shape(shape, max_phase_dev_deg)
        f_ref = _parabolic_peak(freqs, spec, i)
        cands.append({
            "f_hz": f_ref, "T_s": 1.0 / f_ref if f_ref > 0 else float('nan'),
            "rel_power": spec[i] / peak_val if peak_val > 0 else 0.0,
            "phase_dev_deg": dev, "in_phase": bool(ok),
            "is_global_max": bool(abs(spec[i] - peak_val) < 1e-12)})
    coherent = [c for c in cands if c["in_phase"]]

    # --- harmonic-family guard ---------------------------------------
    # Phase coherence across three sensors is NECESSARY but not SUFFICIENT.
    # A harmonic whose nodes happen to fall outside the sensor span reads as
    # in-phase: in the regression test the 4x mode gives phase_dev = 0.1 deg,
    # because sin(4*pi*y/2.58) has the same sign at 0.66, 1.26 and 2.06 m.
    # This is not hypothetical -- NODAMP_v7_NEW returns 0.027 s on runs 10, 12
    # and 15, which is 4x the fundamental.
    #
    # So: if a coherent candidate sits within `tol` of an integer multiple
    # (2..6) of a LOWER coherent candidate, it is that candidate's harmonic and
    # is demoted, however tall it is. Physical period elongation is continuous
    # and never lands on an exact integer ratio to another live peak, so this
    # cannot suppress a genuine measurement.
    # A base must itself carry real energy. Without this floor a near-zero
    # sidelobe becomes an eligible "fundamental" and demotes the true peak:
    # in the 4x-dominant regression case a 1.619 Hz candidate at ~0% power
    # demoted the genuine 9.615 Hz fundamental, because 9.615/1.619 = 5.94 ~ 6.
    tol = 0.03
    BASE_MIN_POWER = 1e-2
    for c in coherent:
        c["harmonic_of"] = None
        for b in coherent:
            if b["rel_power"] < BASE_MIN_POWER:
                continue
            if b["f_hz"] >= c["f_hz"] * (1.0 - tol):
                continue
            k = c["f_hz"] / b["f_hz"]
            k_int = round(k)
            if 2 <= k_int <= 6 and abs(k / k_int - 1.0) < tol:
                c["harmonic_of"] = round(b["f_hz"], 3)
                break

    primary = [c for c in coherent if c["harmonic_of"] is None]
    chosen = max(primary, key=lambda c: c["rel_power"]) if primary else \
        (max(coherent, key=lambda c: c["rel_power"]) if coherent else None)

    if chosen is None:
        return float('nan'), {"candidates": cands,
                              "reason": "no phase-coherent peak in band"}
    return chosen["T_s"], {"candidates": cands, "reason": "",
                           "picked_global_max": chosen["is_global_max"]}


def periodogram_svd(t, accels):
    """PSD matrix -> singular values -> peak-pick.  (The paper's method.)

    With a single data block the 3x3 cross-spectral matrix G(f) = A(f)A(f)^H is
    rank one, so its only non-zero singular value is sum_i |A_i(f)|^2. That is
    stated here rather than hidden: for a single-segment periodogram the SVD
    step reduces exactly to the multi-channel power sum. It is still the right
    estimator -- it uses all three sensors and suppresses a node in any one of
    them -- but it is not doing modal separation, and the write-up should not
    claim it is. Use `welch_svd` below for a genuinely rank>1 estimate.
    """
    dt = float(np.median(np.diff(t)))
    n = len(t)
    nfft = n * ZPAD
    freqs = np.fft.rfftfreq(nfft, d=dt)
    s1 = np.zeros(len(freqs))
    A_list = []
    for a in accels:
        A = np.fft.rfft(prepare(t, a), n=nfft)
        A_list.append(A)
        s1 += np.abs(A) ** 2
    T_shape, info = pick_fundamental(freqs, s1, A_list, zpad=ZPAD)
    # _pick is the naive tallest-peak rule, kept so the two can be compared
    return freqs, s1, (T_shape if np.isfinite(T_shape) else _pick(freqs, s1)), info


def welch_svd(t, accels, nperseg_frac=0.5):
    """Genuine FDD: Welch-averaged 3x3 CPSD matrix, SVD at each frequency.

    Averaging buys a rank>1 matrix at the cost of frequency resolution, which
    on a 1.5-3.0 s record is already the binding constraint. Reported as a
    cross-check, not as the primary number.
    """
    if not _HAVE_SCIPY:
        return None, None, float('nan')
    dt = float(np.median(np.diff(t)))
    fs = 1.0 / dt
    n = len(t)
    nperseg = max(64, int(n * nperseg_frac))
    nfft = nperseg * ZPAD
    m = len(accels)
    G = None
    freqs = None
    for i in range(m):
        for j in range(m):
            f, Gij = _scipy_csd(accels[i], accels[j], fs=fs,
                                nperseg=nperseg, nfft=nfft, detrend='linear')
            if G is None:
                freqs = f
                G = np.zeros((len(f), m, m), dtype=complex)
            G[:, i, j] = Gij
    s1 = np.linalg.svd(G, compute_uv=False)[:, 0]
    return freqs, s1, _pick(freqs, s1)


def fft_per_channel(t, accels):
    """Plain single-channel FFT peak for each accelerometer.

    The paper cross-validated PSD against plain FFT and reported identical
    results. Reproducing that check is the acceptance test: if the three
    channels disagree, the identification is not trustworthy for that run and
    should be flagged rather than propagated.
    """
    dt = float(np.median(np.diff(t)))
    n = len(t)
    nfft = n * ZPAD
    freqs = np.fft.rfftfreq(nfft, d=dt)
    out = []
    for a in accels:
        A = np.abs(np.fft.rfft(prepare(t, a), n=nfft))
        out.append(_pick(freqs, A))
    return out


# ---------------------------------------------------------------------
# Top level
# ---------------------------------------------------------------------
def identify_period_experiment(run_dir, run_label,
                               settle_skip=SETTLE_SKIP_S,
                               max_window=MAX_WINDOW_S,
                               pulse_end_s=None,
                               verbose=True):
    """Identify the fundamental OOP period from one run's exported histories,
    following Moshfeghi et al. (2024) section 3.2.

    Returns a dict. `T` is the primary estimate (periodogram + singular value
    peak-pick over the three accelerometers). `ok` is False when the run should
    be flagged rather than fed back into the driver.
    """
    res = {"run_label": run_label, "T": float('nan'), "ok": False,
           "T_welch": float('nan'), "T_ch": [], "spread": float('nan'),
           "window_s": 0.0, "n_cycles": 0.0, "note": "", "stop_rule": ""}

    t_tab, v_tab = load_channel(run_dir, run_label, CH_TABLE)
    if t_tab is None:
        res["note"] = "missing {}".format(CH_TABLE)
        return res
    t_tab, v_tab = last_segment(t_tab, v_tab)

    accels_raw, missing = [], []
    for ch in CH_WALL:
        tc, vc = load_channel(run_dir, run_label, ch)
        if tc is None:
            missing.append(ch)
            continue
        tc, vc = last_segment(tc, vc)
        accels_raw.append(np.interp(t_tab, tc, vc))
    if missing:
        res["note"] = "missing " + ", ".join(missing)
        if not accels_raw:
            return res

    t_q, v_q, win, stop_src = quiet_window(
        t_tab, accels_raw, v_tab,
        settle_skip=settle_skip, max_window=max_window,
        pulse_end_s=pulse_end_s)
    res["stop_rule"] = stop_src
    if t_q is None:
        res["note"] = stop_src
        return res

    accels = [to_acceleration(t_q, v) for v in v_q]

    _, _, T_svd, pick_info = periodogram_svd(t_q, accels)
    res["candidates"] = pick_info.get("candidates", [])
    res["picked_global_max"] = pick_info.get("picked_global_max", None)
    if pick_info.get("reason"):
        res["note"] = pick_info["reason"]
    _, _, T_wel = welch_svd(t_q, accels)
    T_ch = fft_per_channel(t_q, accels)

    finite = [v for v in T_ch if np.isfinite(v)]
    spread = (max(finite) - min(finite)) / min(finite) if len(finite) > 1 else 0.0

    res.update({"T": T_svd, "T_welch": T_wel, "T_ch": T_ch,
                "spread": spread, "window_s": win,
                "n_cycles": win / T_svd if T_svd and np.isfinite(T_svd) else 0.0})

    # Acceptance test, mirroring the paper's PSD-vs-FFT cross-validation.
    res["ok"] = bool(np.isfinite(T_svd) and spread <= 0.10
                     and not stop_src.startswith("FALLBACK"))
    if stop_src.startswith("FALLBACK") and np.isfinite(T_svd):
        res["note"] = stop_src
    if not res["ok"] and np.isfinite(T_svd):
        res["note"] = "channels disagree by {:.0f}%".format(spread * 100)

    if verbose:
        print("  [exp-method] T = {:.5f} s   window {:.2f} s ({:.0f} cycles)"
              "   [stop rule: {}]".format(T_svd, win, res["n_cycles"], stop_src))
        print("               per-channel FFT: " +
              ", ".join("{:.5f}".format(v) for v in T_ch) +
              "   spread {:.1f}%".format(spread * 100))
        if np.isfinite(T_wel):
            print("               Welch-FDD cross-check: {:.5f} s".format(T_wel))
        cands = res.get("candidates") or []
        if len(cands) > 1 or (cands and not cands[0]["in_phase"]):
            print("               peak candidates (ascending frequency):")
            for c in cands:
                print("                 T={:.5f}s  f={:6.2f}Hz  power={:7.3%}"
                      "  phase_dev={:5.1f}deg  {}{}{}".format(
                          c["T_s"], c["f_hz"], c["rel_power"],
                          c["phase_dev_deg"],
                          "IN-PHASE" if c["in_phase"] else "out-of-phase",
                          "  [tallest]" if c["is_global_max"] else "",
                          "  [harmonic of {} Hz]".format(c["harmonic_of"])
                          if c.get("harmonic_of") else ""))
            if res.get("picked_global_max") is False:
                print("               -> the tallest peak was NOT the "
                      "fundamental; it was rejected on mode shape.")
        if res["note"]:
            print("               NOTE: " + res["note"])
    return res


# ---------------------------------------------------------------------
# Diagnosis -- say WHY a folder produced nothing, do not just return NaN
# ---------------------------------------------------------------------
def diagnose_run_folder(run_dir, run_label):
    """Human-readable report on why identification failed for one run.

    Called automatically when a whole run set comes back empty. Silently
    returning NaN for 25 runs and then plotting the remaining series is how a
    figure ends up looking finished while containing no new information.
    """
    lines = ["  folder: {}".format(run_dir)]
    if not os.path.isdir(run_dir):
        lines.append("  DOES NOT EXIST")
        return "\n".join(lines)

    csvs = sorted(fn for fn in os.listdir(run_dir) if fn.lower().endswith(".csv"))
    lines.append("  {} CSV file(s) present".format(len(csvs)))
    for fn in csvs[:40]:
        lines.append("      {}".format(fn))
    if len(csvs) > 40:
        lines.append("      ... and {} more".format(len(csvs) - 40))

    lines.append("  required channels:")
    for ch in (CH_TABLE,) + CH_WALL:
        hit = find_channel_file(run_dir, run_label, ch)
        lines.append("      {:<24s} {}".format(
            ch, os.path.basename(hit) if hit else "*** NOT FOUND ***"))

    t, v = load_channel(run_dir, run_label, CH_TABLE)
    if t is not None:
        t, v = last_segment(t, v)
        dt = float(np.median(np.diff(t))) if len(t) > 2 else 0.0
        i0, stop_src = table_stop_index(v)
        lines.append("  table channel: {} samples, dt = {:.3e} s, "
                     "record {:.3f} s".format(len(t), dt,
                                              t[-1] - t[0] if len(t) else 0.0))
        lines.append("  stop rule    : {}".format(stop_src))
        lines.append("  table stops at sample {} (t = {:.4f} s) -> quiet window "
                     "{:.4f} s".format(i0, t[i0] if i0 < len(t) else float('nan'),
                                       (t[-1] - t[i0]) if i0 < len(t) else 0.0))
        lines.append("  |v_table| peak {:.4e} m/s, final 5% RMS {:.4e} m/s"
                     .format(float(np.max(np.abs(v))),
                             float(np.sqrt(np.mean(v[int(0.95*len(v)):]**2)))))
        lines.append("  need >= 32 samples in the quiet window; have {}".format(
            max(0, len(t) - i0)))
    return "\n".join(lines)


# ---------------------------------------------------------------------
# Self-test: does this actually fix the failure mode it claims to fix?
# ---------------------------------------------------------------------
def _selftest():
    """Reproduce the run-12 failure on a synthetic signal, then show what
    fixes it. Run this before trusting anything above.

    A post-rocking ring-down is NOT one decaying sinusoid. It is two
    superposed components:

      - a large, heavily damped RIGID-BODY ROCKING transient. Impact at each
        half-cycle dissipates hard, so xi ~ 0.20 and it is gone within a
        second, but while it lasts it is an order of magnitude larger in
        DISPLACEMENT than anything else.
      - a small, lightly damped ELASTIC vibration at the wall's actual
        fundamental period, xi ~ 0.03, which persists for the whole record.

    A displacement periodogram over a short window is dominated by the first.
    That is what T_end = 0.21848 s at run 12 was: not the wall's period, but
    the rocking transient's half-cycle rate, correctly measured off the wrong
    quantity. Note the consequence for TID_MAX_JUMP -- it rejected that value,
    which was the right call for entirely the wrong reason, and it then went
    on to freeze the loop for fourteen runs.

    An earlier version of this self-test used a settling RAMP instead, and both
    methods passed it. Mean removal plus a Hann taper handle a ramp. The ramp
    hypothesis is therefore refuted and is not the mechanism; the two-component
    superposition below is. Kept as a record so it is not re-proposed.
    """
    T_e, T_r = 0.104, 0.250          # elastic period, rocking transient period
    xi_e, xi_r = 0.03, 0.20
    dt, dur = 2.0e-4, 1.5            # 1.5 s = the driver's MAX_RD_WINDOW
    t = np.arange(0.0, dur, dt)
    we, wr = 2 * math.pi / T_e, 2 * math.pi / T_r

    shape = [math.sin(math.pi * y / 2.58) for y in (0.66, 1.26, 2.06)]
    base = (20.0e-3 * np.exp(-xi_r * wr * t) * np.sin(wr * t) +
             1.5e-3 * np.exp(-xi_e * we * t) * np.sin(we * t))
    disp = [s * base for s in shape]
    vel = [np.gradient(d, t) for d in disp]

    def peak_of(tt, yy):
        nfft = len(yy) * ZPAD
        f = np.fft.rfftfreq(nfft, d=float(np.median(np.diff(tt))))
        return _pick(f, np.abs(np.fft.rfft(prepare(tt, yy), n=nfft)))

    print("Self-test -- truth: elastic {:.3f} s, rocking transient {:.3f} s"
          .format(T_e, T_r))
    print()
    acc = [to_acceleration(t, v) for v in vel]
    half = len(t) // 2

    rows = [
        ("displacement, full 1.5 s  (CURRENT DRIVER)", peak_of(t, disp[2])),
        ("displacement, last 50%",                     peak_of(t[half:], disp[2][half:])),
        ("acceleration, full 1.5 s  (THIS MODULE)",    periodogram_svd(t, acc)[2]),
        ("acceleration, last 50%",                     periodogram_svd(t[half:], [a[half:] for a in acc])[2]),
    ]
    for name, T in rows:
        tag = "  <-- locks onto the rocking transient" if abs(T - T_r) < 0.02 else ""
        print("  {:<44s} {:.4f} s{}".format(name, T, tag))
    print()
    print("  Switching to acceleration is sufficient on its own: the w^2")
    print("  weighting suppresses the low-frequency transient without")
    print("  discarding any of the record. Cropping to the late window also")
    print("  works and is what a 45-60 s experimental record does implicitly,")
    print("  but it throws away data the model can barely afford to lose.")
    return rows


if __name__ == "__main__":
    _selftest()


# ---------------------------------------------------------------------
# Modal tracking -- the selection rule, VALIDATED against the experiment
# ---------------------------------------------------------------------
# Taking the tallest phase-coherent peak is wrong, and the experiment's own
# records prove it. Run 1 of Test 9 carries a coherent peak at 10.97 Hz
# (T = 0.09116 s, 74% power) -- Table 5's 0.091 s -- and a TALLER coherent peak
# at 18.76 Hz. The ratio 1.71 is not an integer, so it is a genuine second mode
# (almost certainly the timber floor: channels 13/14 are the slab
# accelerometers and the joists are embedded in the wall), not a harmonic any
# integer-ratio guard can catch. "Tallest" returns 0.0533 s, off by 41%.
#
# Seeding on run 1 and tracking the same mode forward reproduces Table 5:
#
#     run  1   0.09116 s   vs 0.091   (+0.2%)
#     run 24   0.10844 s   vs 0.107   (+1.3%)
#     elongation +19.0%    vs +17.6%
#
# Two independent ground-truth points, both within 1.3%. That is the closest
# this project has come to validating a measurement against the source paper,
# and it is why this rule -- not the tallest-peak rule -- is the one to use.
#
# LIMITATION, state it in the write-up: tracking is sequential and needs a
# seed, so an error in the seed propagates. Guard by checking the seed run
# against a known initial period, and by requiring the tracked series to stay
# phase-coherent every run (both are reported below).

TRACK_SEED_MIN_POWER = 0.50   # seed = lowest coherent peak carrying >= this
TRACK_MIN_POWER      = 0.05   # candidates considered when tracking
TRACK_SHRINK_PENALTY = 0.35   # log-distance penalty for a candidate BELOW the
                              # running value. Damage lengthens a period; it
                              # does not shorten one. This biases against
                              # dropping onto a higher mode without forbidding
                              # a genuine decrease outright.


def seed_period(candidates, min_power=TRACK_SEED_MIN_POWER):
    """Lowest phase-coherent peak carrying at least `min_power`."""
    c = [x for x in candidates if x["in_phase"] and x["rel_power"] >= min_power]
    if not c:
        c = [x for x in candidates if x["in_phase"]]
    return min(c, key=lambda x: x["f_hz"])["T_s"] if c else float("nan")


def track_step(candidates, T_prev, min_power=TRACK_MIN_POWER,
               shrink_penalty=TRACK_SHRINK_PENALTY):
    """Follow the same mode into the next run. Returns (T, record or None)."""
    c = [x for x in candidates if x["in_phase"] and x["rel_power"] >= min_power]
    if not c or not (T_prev and np.isfinite(T_prev)):
        return float("nan"), None
    best = min(c, key=lambda x: abs(math.log(x["T_s"] / T_prev))
               + (shrink_penalty if x["T_s"] < T_prev * 0.97 else 0.0))
    return best["T_s"], best