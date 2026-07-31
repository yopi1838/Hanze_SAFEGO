# -*- coding: ascii -*-
"""
Batch plot export from saved states (run inside 3DEC)
=====================================================
Restores every .sav in a folder, and for each one exports the named plot
views to high-res bitmaps into a "plots" subfolder of that same folder.

Run inside 3DEC (from anywhere):
    call 'export_plots.py'
Set SAV_DIR below to the folder that holds the .sav files (or leave "."
and run from inside that folder).

The plot views named in PLOTS must already exist in the restored state
(or be (re)created by PLOT_DEF_FILE). Each export uses:
    plot '<view>' export bitmap filename '<out>' size W H dpi DPI
"""
import itasca as it
import os, glob, re

# ===================== CONFIG =====================
SAV_DIR      = "stratC_results_RATCHETING"                       # folder with the .sav files ("." = cwd)
SAV_GLOB     = "stratC_run_*.sav"        # which saves to process ("*.sav" for all)
OUT_SUBDIR   = "plots"                   # created inside SAV_DIR
# (3DEC plot-view name, output-file label):
PLOTS        = [("Master_", "Master"), ("Plane", "Plane"),("History_", "History")]
IMG_W        = 3840
IMG_H        = 2160
IMG_DPI      = 600
PLOT_DEF_FILE = ""                       # optional .dat that (re)builds the views; "" to skip
# ==================================================

def cmd_path(p):
    # forward slashes work on Windows and Linux; avoids backslash issues
    return os.path.normpath(p).replace("\\", "/")

def run_tag(fname):
    # use the last run-number in the filename, zero-padded to 2 digits
    nums = re.findall(r'\d+', os.path.basename(fname))
    return nums[-1].zfill(2) if nums else os.path.splitext(os.path.basename(fname))[0]

def main():
    sav_dir = os.path.abspath(SAV_DIR)
    out_dir = os.path.join(sav_dir, OUT_SUBDIR)
    os.makedirs(out_dir, exist_ok=True)

    savs = sorted(glob.glob(os.path.join(sav_dir, SAV_GLOB)))
    if not savs:
        print("No saves matching '{}' in {}".format(SAV_GLOB, sav_dir))
        return

    print("Found {} save(s). Exporting {}x{} @ {} dpi to: {}".format(
        len(savs), IMG_W, IMG_H, IMG_DPI, out_dir))

    ok, fail = 0, 0
    for sav in savs:
        tag = run_tag(sav)
        try:
            it.command("model restore '{}'".format(cmd_path(sav)))
        except Exception as e:
            print("  [SKIP] cannot restore {}: {}".format(os.path.basename(sav), e))
            fail += 1
            continue
        if PLOT_DEF_FILE:
            it.command("call '{}'".format(cmd_path(PLOT_DEF_FILE)))
        for view, label in PLOTS:
            out_png = cmd_path(os.path.join(out_dir, "{}_Run{}.png".format(label, tag)))
            try:
                it.command("plot '{}' export bitmap filename '{}' size {} {} dpi {}".format(
                    view, out_png, IMG_W, IMG_H, IMG_DPI))
                print("  Run {}: {}".format(tag, os.path.basename(out_png)))
                ok += 1
            except Exception as e:
                print("  [WARN] Run {} plot '{}' failed: {}".format(tag, view, e))
                fail += 1

    print("Done. {} image(s) exported, {} skipped/failed. Folder: {}".format(ok, fail, out_dir))

main()
