"""
reid_periods_exp.py -- re-identify the fundamental period of an ALREADY
                       COMPLETED run set using the experiment's method.

No 3DEC rerun. Reads only the exported history CSVs that every run already
wrote (Channel_12/15/16/17), so any of NODAMP_v6_NEW, v7_NEW, SECANT, stratD,
G00_LS can be re-read today.

Run it the same way as postprocess_stratC.py -- edit the constants and
`python reid_periods_exp.py`, or `call 'reid_periods_exp.py'` from 3DEC (no
argv is passed in that case, which is why these are constants and not
argparse).

Outputs, into <SIM_DIR>/postproc/:
  period_comparison.csv   one row per run
  period_comparison.png   T vs run number, both estimators, with Table 5
                          of Moshfeghi et al. (2024) overlaid for US-1
"""

import os
import csv
import numpy as np

import period_id_exp as P

# ---------------------------------------------------------------------
SIM_DIR = "stratF_full_results_US1"   # <-- edit
LABEL = "Full_results_noDamp"                    # <-- edit
DRIVER_LOG = "strategy_C_log.csv"          # written by the driver, optional

# MEASURED experimental period, per run. This file already exists in the repo
# and is what postprocess_stratC.py plots -- do not regenerate it, and do not
# substitute a straight line between Table 5's two endpoints, which is what an
# earlier version of this script drew and which made both estimators look wrong
# for reasons that were an artefact of the reference.
# Fig 14 of the paper, digitised from its own gridlines (digitize_fig14.py).
# Calibration validates against Table 5 at BOTH ends: run 1 = 0.09110 vs 0.091,
# run 24 = 0.10749 vs 0.107, elongation +18.0% vs +17.6%. This is the published
# experiment and is the reference.
#
# The two RE-DERIVED series are kept for context but are not the benchmark:
#   exp_Test9_period_psd.csv  mean |err| 1.51% vs Fig 14, but runs 12-15 are
#                             frozen at 0.10059 and runs 23-24 overshoot by
#                             ~5% (0.11210 at run 24, above anything in the
#                             paper).
#   US1_experimental_periods.csv (this project's tracking rule)
#                             mean |err| 3.26%, worst +8.1% at run 9. WORSE
#                             than the file already in the repo. Do not
#                             substitute it for either.
# Both unstrengthened specimens. US-1 is the specimen being modelled; US-2 is
# the companion reference wall and is plotted for context, because the two
# bracket the scatter that "nominally identical masonry" actually produces.
#
#   US-1   0.0911 -> 0.1075   +18.0%   (paper states +17.5%)
#   US-2   0.0875 -> 0.0994   +13.7%   (paper states +13.9%)
#
# Plotted as DOTS WITH A REGRESSION LINE, the way Fig 14 itself presents them.
# Joining consecutive points implies a run-to-run trajectory the measurement
# does not support: the scatter between adjacent runs is comparable to the
# total elongation over several runs, and a few markers are recovered from
# occlusions (flagged px_area = -1, drawn as open symbols).
EXP_SERIES = [
    ("US-1", "US1_fig14_digitised.csv", "k",    "o"),
    ("US-2", "US2_fig14_digitised.csv", "0.45", "s"),
]
EXP_T_INITIAL = 0.091                             # Table 5, endpoints only
EXP_T_FINAL = 0.107
# ---------------------------------------------------------------------


def run_folders(sim_dir):
    if not os.path.isdir(sim_dir):
        raise SystemExit("No such directory: {}".format(sim_dir))
    out = []
    for d in sorted(os.listdir(sim_dir)):
        p = os.path.join(sim_dir, d)
        if os.path.isdir(p) and d.startswith("Run"):
            try:
                n = int(d[3:5])
            except ValueError:
                continue
            out.append((n, d, p))
    return sorted(out)


def driver_log(sim_dir):
    """Read strategy_C_log.csv -> {run: {"T_end":..., "pulse_end":...}}.

    pulse_end = n_cycles_used * T_excitation is the EXACT instant the driver's
    velocity table returns to zero, i.e. the instant the table stops. Using it
    is the direct equivalent of the experiment reading that instant off its
    table displacement sensor, and unlike a threshold on Channel_12 it does not
    depend on which gridpoint that history landed on.
    """
    path = os.path.join(sim_dir, DRIVER_LOG)
    if not os.path.exists(path):
        return {}
    got = {}
    with open(path) as f:
        for row in csv.DictReader(f):
            try:
                n = int(row["run"])
            except (KeyError, ValueError, TypeError):
                continue
            rec = {}
            try:
                rec["T_end"] = float(row["T_end"])
            except (KeyError, ValueError, TypeError):
                pass
            try:
                rec["pulse_end"] = (float(row["n_cycles_used"]) *
                                    float(row["T_excitation"]))
            except (KeyError, ValueError, TypeError):
                pass
            got[n] = rec
    return got


def main():
    folders = run_folders(SIM_DIR)
    if not folders:
        raise SystemExit("No Run** folders in {}".format(SIM_DIR))
    logged = driver_log(SIM_DIR)
    if not logged:
        print("NOTE: {} not found. Falling back to the RMS table-stop rule, "
              "which is weaker -- runs found that way are marked in the "
              "stop_rule column.".format(DRIVER_LOG))

    rows = []
    T_track = None      # running fundamental, seeded on the first good run
    print("Re-identifying periods in {}  ({} runs)".format(SIM_DIR, len(folders)))
    print("Method: Moshfeghi et al. (2024) sec. 3.2 -- three accelerometers,")
    print("        cropped at table stop, PSD -> singular values -> peak-pick.")
    print()
    for n, label, path in folders:
        print("Run {:02d}  {}".format(n, label))
        res = P.identify_period_experiment(
            path, label, pulse_end_s=logged.get(n, {}).get("pulse_end"),
            verbose=True)
        # --- selection rule: seed once, then TRACK the same mode ---------
        # Validated against the experiment's own records: seeding on run 1 and
        # tracking forward reproduces Table 5 to +0.2% (run 1) and +1.3%
        # (run 24). Taking the tallest peak returns 0.0533 s on run 1, off by
        # 41%, because the timber floor mode at 18.76 Hz is taller than the
        # wall's own 10.97 Hz. See period_id_exp.track_step.
        cl = res.get("candidates") or []
        if cl:
            if T_track is None:
                T_track = P.seed_period(cl)
                res["T"] = T_track
                res["ok"] = bool(np.isfinite(T_track))
                res["note"] = (res["note"] + " | seed").strip(" |")
            else:
                T_new, rec = P.track_step(cl, T_track)
                if np.isfinite(T_new):
                    res["T"] = T_new
                    res["note"] = (res["note"] + " | tracked p={:.0%}".format(
                        rec["rel_power"])).strip(" |")
                    # res["ok"] keeps the channel-agreement verdict; the anchor
                    # only advances on an accepted run.
                    if res["ok"]:
                        T_track = T_new
        rows.append({
            "run": n,
            "folder": label,
            "T_exp_psd_s": "" if not np.isfinite(res["T"]) else round(res["T"], 6),
            "T_exp_welch_s": "" if not np.isfinite(res["T_welch"]) else round(res["T_welch"], 6),
            "T_ch_bot_s": "" if not res["T_ch"] else round(res["T_ch"][0], 6),
            "T_ch_mid_s": "" if len(res["T_ch"]) < 2 else round(res["T_ch"][1], 6),
            "T_ch_top_s": "" if len(res["T_ch"]) < 3 else round(res["T_ch"][2], 6),
            "channel_spread": "" if not np.isfinite(res["spread"]) else round(res["spread"], 4),
            "window_s": round(res["window_s"], 4),
            "cycles_in_window": round(res["n_cycles"], 1),
            "accepted": int(res["ok"]),
            "stop_rule": res.get("stop_rule", ""),
            "note": res["note"],
            "pulse_end_s": round(logged.get(n, {}).get("pulse_end", float("nan")), 5)
                           if "pulse_end" in logged.get(n, {}) else "",
            "T_driver_logged_s": logged.get(n, {}).get("T_end", ""),
        })
        print()

    # -----------------------------------------------------------------
    # Refuse to produce a figure that contains no new information.
    #
    # If every run returned NaN, the plot below would still draw the driver's
    # logged series and the Table 5 reference line, and would look like a
    # finished result. It is not one. Stop here and say why.
    # -----------------------------------------------------------------
    produced = [r for r in rows if r["T_exp_psd_s"] != ""]
    if not produced:
        print()
        print("=" * 70)
        print("NO RUN PRODUCED AN EXPERIMENT-METHOD PERIOD.")
        print("Nothing has been plotted. The reasons reported were:")
        print("=" * 70)
        tally = {}
        for r in rows:
            tally[r["note"] or "(no note)"] = tally.get(r["note"] or "(no note)", 0) + 1
        for note, k in sorted(tally.items(), key=lambda kv: -kv[1]):
            print("  {:3d} run(s): {}".format(k, note))
        print()
        print("Folder inspection for the first run:")
        print(P.diagnose_run_folder(folders[0][2], folders[0][1]))
        print()
        print("Most likely causes, in order:")
        print("  1. The accelerometer channels were never exported for this run")
        print("     set. instrument_history_export_new.dat lines 22-27 export")
        print("     Channel_12/13/14/15/16/17; instrument_history_export_v2.dat")
        print("     may not. Check which export file export_all_histories calls.")
        print("  2. The run label used in the filenames differs from the folder")
        print("     name. This version falls back to a folder scan, so that")
        print("     should no longer bite -- if it still does, the files are")
        print("     genuinely absent.")
        print("  3. The quiet window is shorter than 32 samples, i.e. the table")
        print("     never stopped inside the recorded history.")
        raise SystemExit(1)

    outdir = os.path.join(SIM_DIR, "postproc")
    if not os.path.isdir(outdir):
        os.makedirs(outdir)
    csv_path = os.path.join(outdir, "period_comparison.csv")
    with open(csv_path, "w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=list(rows[0].keys()))
        w.writeheader()
        w.writerows(rows)
    print("wrote {}".format(csv_path))

    # ---- summary against Table 5 -------------------------------------
    acc = [r for r in rows if r["accepted"] and r["T_exp_psd_s"] != ""]
    if acc:
        T0, T1 = acc[0]["T_exp_psd_s"], acc[-1]["T_exp_psd_s"]
        print()
        print("  {:<28s} {:>8s} {:>8s} {:>10s}".format(
            "", "initial", "final", "elongation"))
        print("  {:<28s} {:8.4f} {:8.4f} {:9.1f}%".format(
            "model (experiment method)", T0, T1, 100 * (T1 / T0 - 1)))
        print("  {:<28s} {:8.4f} {:8.4f} {:9.1f}%".format(
            "US-1, Table 5", EXP_T_INITIAL, EXP_T_FINAL,
            100 * (EXP_T_FINAL / EXP_T_INITIAL - 1)))
        rejected = [r["run"] for r in rows if not r["accepted"]]
        if rejected:
            print()
            print("  {} run(s) failed the PSD-vs-FFT cross-check and are "
                  "excluded above: {}".format(len(rejected), rejected))
            print("  That is the paper's own acceptance test, not an extra "
                  "filter. Do not quietly average them back in.")

    # ---- figure ------------------------------------------------------
    try:
        import matplotlib
        matplotlib.use("Agg")
        import matplotlib.pyplot as plt
    except Exception:
        print("matplotlib unavailable; CSV written, figure skipped.")
        return

    runs = [r["run"] for r in rows]
    fig, ax = plt.subplots(figsize=(9, 5))
    fits = {}

    # ---- model, experiment method: dots + regression -------------------
    # Same presentation as the published series, for the same reason: the
    # run-to-run scatter is comparable to the elongation over several runs,
    # so a connecting line implies a trajectory the measurement does not
    # support. The comparison that carries meaning is slope against slope.
    #
    # The "driver as logged" series is deliberately NOT plotted. It is a
    # single-channel DISPLACEMENT estimate, and over runs 11-25 of the NODAMP
    # set it reported one held value to five decimals -- the TID_AGREE_TOL
    # hold, not a measurement. It stays in the CSV as T_driver_logged_s for
    # audit, but it does not belong on a figure next to a measured curve.
    ex = [(r["run"], r["T_exp_psd_s"]) for r in rows
          if r["accepted"] and r["T_exp_psd_s"] != ""]
    if ex:
        mx = np.array([a for a, _ in ex], float)
        my = np.array([b for _, b in ex], float)
        ax.plot(mx, my, "o", ls="none", color="tab:blue", ms=5, zorder=6,
                label="model, experiment method (acc, 3 ch)")
        if len(mx) >= 2:
            k, c = np.polyfit(mx, my, 1)
            xs = np.array([mx.min(), mx.max()])
            g = 100.0 * ((k * xs[1] + c) / (k * xs[0] + c) - 1.0)
            fits["model"] = (k, c, g, xs)
            ax.plot(xs, k * xs + c, "-", color="tab:blue", lw=1.6, zorder=5,
                    label="model fit  ({:+.1f}% over {:.0f} runs)".format(
                        g, xs[1] - xs[0] + 1))

    # Rejected runs are shown but excluded from the fit: they failed the
    # paper's own PSD-vs-FFT cross-check, so they are not measurements.
    bad = [(r["run"], r["T_exp_psd_s"]) for r in rows
           if not r["accepted"] and r["T_exp_psd_s"] != ""]
    if bad:
        ax.plot([a for a, _ in bad], [b for _, b in bad], "x",
                color="0.6", ms=5, zorder=6,
                label="failed cross-check (excluded from fit)")

    # ---- published experimental series: dots + regression, as in Fig 14 ----
    drew_exp = False
    for name, path, colour, marker in EXP_SERIES:
        if not os.path.exists(path):
            print("NOTE: {} not found -- {} not plotted.".format(path, name))
            continue
        try:
            er, et, flag = [], [], []
            with open(path) as f:
                for row in csv.DictReader(f):
                    er.append(int(row["run"]))
                    et.append(float(row["T_exp_s"]))
                    flag.append(float(row.get("px_area", 0) or 0) < 0)
        except Exception as e:
            print("could not read {}: {}".format(path, e))
            continue
        er = np.array(er, float); et = np.array(et, float)
        flag = np.array(flag, bool)

        # clean reads filled, occlusion-recovered points open
        ax.plot(er[~flag], et[~flag], marker, ls="none", color=colour, ms=5,
                zorder=5, label="{} published (Fig 14)".format(name))
        if flag.any():
            ax.plot(er[flag], et[flag], marker, ls="none", ms=5, zorder=5,
                    mfc="none", mec=colour,
                    label="{} recovered from occlusion".format(name))

        # least-squares trend, which is what Fig 14 draws through its dots
        k, c = np.polyfit(er, et, 1)
        xs = np.array([er.min(), er.max()])
        growth = 100.0 * ((k * xs[1] + c) / (k * xs[0] + c) - 1.0)
        ax.plot(xs, k * xs + c, "-", color=colour, lw=1.6, alpha=0.85, zorder=4,
                label="{} fit  ({:+.1f}% over {:.0f} runs)".format(
                    name, growth, xs[1] - xs[0] + 1))
        fits[name] = (k, c, growth, xs)
        drew_exp = True

    if not drew_exp:
        ax.plot([min(runs), max(runs)], [EXP_T_INITIAL, EXP_T_FINAL], "k--",
                lw=1.5, alpha=0.6,
                label="Table 5 endpoints ONLY (interpolated -- not data)")
        print("NOTE: no digitised series found. The dashed reference is a "
              "straight line between two endpoints, NOT the experimental "
              "curve. Do not read run-by-run agreement off it.")

    ax.set_xlabel("Run")
    ax.set_ylabel("Fundamental OOP period [s]")
    ax.set_title("Period identification -- {}".format(LABEL))
    ax.grid(alpha=0.3)
    ax.legend(fontsize=7, ncol=2)
    fig.tight_layout()
    png = os.path.join(outdir, "period_comparison.png")
    fig.savefig(png, dpi=150)
    print("wrote {}".format(png))

    # ---- the headline number: slope against slope ----------------------
    if fits:
        print()
        print("Period elongation, least-squares fit over the run sequence")
        print("  {:<10s} {:>10s} {:>10s} {:>12s}".format(
            "", "T at run 1", "T at last", "elongation"))
        for key in ("model", "US-1", "US-2"):
            if key not in fits:
                continue
            k, c, g, xs = fits[key]
            print("  {:<10s} {:10.5f} {:10.5f} {:11.1f}%".format(
                key, k * xs[0] + c, k * xs[1] + c, g))
        if "model" in fits and "US-1" in fits:
            print()
            print("  model elongates at {:.0f}% of the US-1 rate.".format(
                100 * fits["model"][2] / fits["US-1"][2]))
        if "US-1" in fits and "US-2" in fits:
            k1, c1, _, x1 = fits["US-1"]
            k2, c2, _, x2 = fits["US-2"]
            spread = 100 * abs((k1 * x1[0] + c1) / (k2 * x2[0] + c2) - 1)
            print("  US-1 and US-2 differ by {:.1f}% in fitted initial period."
                  .format(spread))
            print("  That is the scatter two nominally identical unstrengthened")
            print("  walls actually produce. Read the model's initial-period")
            print("  offset against that band, not against zero.")


if __name__ == "__main__":
    main()