"""compare_processed_old_new.py -- old vs new TestXXRunYY_processed_globalzero.xlsx, run by run.

Old = the files made with 2626 kg (midpoint split x 2.317) and ch12 added to the wall channels.
New = make_processed_globalzero_v2.py (group's masses 455/670/606 kg, absolute channels, no factor).

Checks and shows:
  (a) peak base shear per run, old vs new, and the ratio old/new  (should be ~constant if only the
      masses changed; it is NOT constant when the double-counted table term matters, i.e. on strong runs)
  (b) U_avg must be identical in both (the displacement recipe did not change) -- max |difference| printed
  (c) base shear vs displacement loops, old and new overlaid, for the runs you name
Usage:
    python compare_processed_old_new.py <old_dir> <new_dir> Test9 [22 23 24]
Output: <new_dir>/compare_old_new_TestXX.png and compare_old_new_TestXX.csv
"""
import sys, os, re, glob, csv
import numpy as np
import openpyxl
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

def read_proc(path):
    ws = openpyxl.load_workbook(path, read_only=True).worksheets[0]
    it = ws.iter_rows(values_only=True); hdr = [str(h).strip() if h is not None else "" for h in next(it)]
    ci = {h: j for j, h in enumerate(hdr)}
    rows = [r for r in it if r[ci["Time"]] is not None]
    g = lambda k: np.array([r[ci[k]] for r in rows], float)
    return g("Time"), g("U_avg") * 1e3, g("base_shear_kN")

def run_no(p):
    m = re.search(r"Run0*(\d+)_", os.path.basename(p)); return int(m.group(1)) if m else None

def main():
    if len(sys.argv) < 4:
        print(__doc__); sys.exit(1)
    old_dir, new_dir, test = sys.argv[1:4]; loop_runs = [int(x) for x in sys.argv[4:]] or [22, 23, 24]
    new_files = {run_no(p): p for p in glob.glob(os.path.join(new_dir, test + "Run*_processed_globalzero.xlsx")) if run_no(p)}
    old_files = {run_no(p): p for p in glob.glob(os.path.join(old_dir, test + "Run*_processed_globalzero.xlsx")) if run_no(p)}
    runs = sorted(set(new_files) & set(old_files))
    if not runs: raise FileNotFoundError("no common runs between {} and {}".format(old_dir, new_dir))
    rows = []; data = {}
    print("{:>4} {:>10} {:>10} {:>8} {:>12}".format("run", "old kN", "new kN", "old/new", "max|dU| mm"))
    for r in runs:
        to, uo, vo = read_proc(old_files[r]); tn, un, vn = read_proc(new_files[r])
        n = min(len(uo), len(un)); dU = float(np.nanmax(np.abs(uo[:n] - un[:n])))
        po, pn = float(np.nanmax(np.abs(vo))), float(np.nanmax(np.abs(vn)))
        rows.append([r, po, pn, po / pn if pn else np.nan, dU]); data[r] = (un, vn, uo, vo)
        print("{:4d} {:10.2f} {:10.2f} {:8.2f} {:12.2e}".format(r, po, pn, po / pn if pn else np.nan, dU))
    with open(os.path.join(new_dir, "compare_old_new_{}.csv".format(test)), "w", newline="") as f:
        w = csv.writer(f); w.writerow(["run", "old_peak_kN", "new_peak_kN", "old_over_new", "max_abs_dU_mm"]); w.writerows(rows)
    R = np.array(rows)
    fig = plt.figure(figsize=(16, 8)); gs = fig.add_gridspec(2, max(len(loop_runs), 2))
    a = fig.add_subplot(gs[0, :]); x = R[:, 0]
    a.plot(x, R[:, 1], "s-", color="#8a8f97", label="old: 2626 kg x (ch + ch12)"); a.plot(x, R[:, 2], "o-", color="#b4491a", label="new: 1732 kg, absolute channels")
    a2 = a.twinx(); a2.plot(x, R[:, 3], "^--", color="#2b5f8a", alpha=0.8, label="old / new"); a2.set_ylabel("old / new", color="#2b5f8a")
    a.set_xticks(x); a.set_xlabel("run"); a.set_ylabel("peak base shear (kN)"); a.grid(alpha=0.25)
    a.set_title("{}: peak experimental base shear per run, old vs new processing".format(test), loc="left", fontweight="bold")
    h1, l1 = a.get_legend_handles_labels(); h2, l2 = a2.get_legend_handles_labels(); a.legend(h1 + h2, l1 + l2, fontsize=8, frameon=False, loc="upper left")
    for k, r in enumerate(loop_runs):
        b = fig.add_subplot(gs[1, k])
        if r not in data:
            b.text(0.5, 0.5, "run {} not in both folders".format(r), ha="center", transform=b.transAxes); continue
        un, vn, uo, vo = data[r]
        b.plot(uo, vo, lw=0.6, color="#8a8f97", alpha=0.8, label="old"); b.plot(un, vn, lw=0.8, color="#b4491a", label="new")
        b.axhline(0, color="k", lw=0.5); b.axvline(0, color="k", lw=0.5); b.grid(alpha=0.25)
        b.set_title("run {}".format(r), loc="left", fontweight="bold"); b.set_xlabel("U_avg (mm)"); b.set_ylabel("base shear (kN)"); b.legend(fontsize=8, frameon=False)
    fig.tight_layout(); out = os.path.join(new_dir, "compare_old_new_{}.png".format(test)); fig.savefig(out, dpi=150); print("->", out)

if __name__ == "__main__":
    main()