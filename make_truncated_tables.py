# -*- coding: ascii -*-
"""
make_truncated_tables.py -- shortened record velocity tables for Strategy F: the quiet lead-in is cut, the cut start
is faded in, the end is faded out to zero, so 3DEC spends its time on the strong motion. The free-vibration tail
after each record is the driver's job (TAIL_SEC_RECORD), not this file's.

    vel_HU.txt -> vel_HU_tr.txt      cut at CUT_S["HU12"]   (35.9 s -> ~3.3 s)
    vel_EC.txt -> vel_EC_tr.txt      cut at CUT_S["EC40"]   (39.6 s -> ~26 s)
    vel_FR.txt -> vel_FR_tr.txt      cut at CUT_S["FR76"]   (17.3 s -> ~9 s)

Per record: keep t >= CUT, shift to start at 0, cosine fade-in over FADE_IN_S, cosine fade-out over FADE_OUT_S (the
originals end with a one-step drop to zero velocity, i.e. an acceleration spike), same two-line header as the
originals (name, "N<tab>0"). The original files are only read.
Check printed per record: duration saved, the largest velocity and the table displacement in the dropped part, and
the 1.5 %-damped Sd of original and truncated at T1 (0.089, 0.103 s) and T_eff (0.36 s) -- they should agree to ~0.5 %.

USAGE   python make_truncated_tables.py            (in the folder with vel_HU.txt / vel_EC.txt / vel_FR.txt)
"""
import os
import numpy as np

RECORDS = {"HU12": ("vel_HU.txt", "vel_HU_tr.txt", "vel_HU_tr"),
           "EC40": ("vel_EC.txt", "vel_EC_tr.txt", "vel_EC_tr"),
           "FR76": ("vel_FR.txt", "vel_FR_tr.txt", "vel_FR_tr")}
CUT_S      = {"HU12": 32.6, "EC40": 13.5, "FR76": 8.0}   # FR76: 8 s keeps more lead-in than the 9.4 s minimum
FADE_IN_S  = 0.5
FADE_OUT_S = 0.25
CHECK_T    = (0.089, 0.103, 0.36)
CHECK_XI   = 0.015


def read(path):
    with open(path) as f:
        f.readline(); f.readline()
        d = np.loadtxt(f)
    return d[:, 0], d[:, 1]


def write(path, name, t, v):
    with open(path, "w", newline="\n") as f:
        f.write("{}\n{}\t0\n".format(name, len(t)))
        for ti, vi in zip(t, v):
            f.write("{:.6f}\t{:.9e}\n".format(ti, vi))


def sd(t, v, T, xi=CHECK_XI, tail=3.0, dt=1e-3):
    tt = np.arange(t[0], t[-1] + tail, dt)
    a = np.interp(tt, t, np.gradient(v, t), right=0.0)
    w = 2 * np.pi / T; k = w * w; c = 2 * xi * w; keff = k + 2 * c / dt + 4 / dt ** 2
    u = vv = um = 0.0; ac = -a[0]
    for i in range(1, len(a)):
        p = -a[i] + (4 / dt ** 2 * u + 4 / dt * vv + ac) + c * (2 / dt * u + vv)
        un = p / keff; vn = 2 / dt * (un - u) - vv; ac = 4 / dt ** 2 * (un - u) - 4 / dt * vv - ac
        u, vv = un, vn; um = max(um, abs(u))
    return um * 1e3


def main():
    for rec, (src, dst, name) in RECORDS.items():
        if not os.path.isfile(src):
            print("! {}: {} missing".format(rec, src)); continue
        if os.path.abspath(src) == os.path.abspath(dst):
            raise SystemExit("refusing to overwrite the original " + src)
        t, v = read(src)
        dt = float(np.median(np.diff(t)))
        m = t >= CUT_S[rec]
        tc = t[m] - t[m][0]; vc = v[m].copy()
        n_in, n_out = int(round(FADE_IN_S / dt)), int(round(FADE_OUT_S / dt))
        vc[:n_in] *= 0.5 * (1 - np.cos(np.pi * np.arange(n_in) / n_in))
        vc[-n_out:] *= 0.5 * (1 + np.cos(np.pi * np.arange(n_out) / n_out))
        vc[-1] = 0.0
        write(dst, name, tc, vc)
        drop_v = np.max(np.abs(v[~m])) / np.max(np.abs(v)) * 100 if np.any(~m) else 0.0
        drop_d = np.trapezoid(v[~m], t[~m]) * 1e3 if np.any(~m) else 0.0
        chk = ", ".join("T {:.3f}: {:.3f} / {:.3f} mm".format(T, sd(t, v, T), sd(tc, vc, T)) for T in CHECK_T)
        print("{}: {} -> {}  {:.1f} -> {:.1f} s (-{:.0f} %), dropped part: max |v| {:.2f} % of PGV, table disp {:.3f} mm"
              .format(rec, src, dst, t[-1], tc[-1], 100 * (1 - tc[-1] / t[-1]), drop_v, drop_d))
        print("      Sd {:.1f} % original / truncated: {}".format(CHECK_XI * 100, chk))


if __name__ == "__main__":
    main()
