# Hanze_SAFEGO

3DEC (Distinct Element Method) modelling of an out-of-plane loaded unreinforced
masonry wall with timber floor joists, cross-validated against shake-table tests
US-1 (Test 9) and US-2 (Test 12).

Full project context — pipeline, conventions, model parameters, known issues — lives
in **[`CLAUDE.md`](CLAUDE.md)**. Read that first.

---

## Setting this up on a new laptop

### 1. Clone

```bash
git clone <your-remote-url> Hanze_SAFEGO
cd Hanze_SAFEGO
```

### 2. Python environment (for post-processing, outside 3DEC)

```bash
python -m venv .venv
# Windows:  .venv\Scripts\activate
# Linux/macOS:  source .venv/bin/activate
pip install -r requirements.txt
```

### 3. 3DEC

Install **Itasca 3DEC 9.1**. The custom joint constitutive model `mason_v6` is a
compiled plugin — `libjmodelmason009.so` / `libjmodelmasonv6009.so` in this repo are
the **Linux** builds. On Windows you need the equivalent `.dll` on 3DEC's plugin
search path, otherwise `block contact jmodel assign mason_v6` fails.

Launch 3DEC with this folder as the working directory — every driver resolves its
inputs relative to cwd.

### 4. What is *not* in this repo

~10 GB of simulation state was excluded to keep the repo usable:

- `*.sav` — 3DEC saved states (`Part_I_MASON*.sav`, `GEO_WALL.sav`,
  `stratC_results_*/stratC_run_NN.sav`, 47–101 MB each)
- `stratC_results_*/RunNN_*/` — raw per-channel history CSVs (~2.5 MB × ~23 per run)
- `stratC_results_*/plots/` — 3DEC bitmap dumps (regenerate with `export_plots.py`)
- `P.prj`, `P.temp`, `P.backup` — machine-specific 3DEC project files
- The experimental `EXP_DATA` xlsx files (never lived in this folder)

**What *is* here:** all code, the 3DEC model definition, `Groningen.dec`, the target
spectra, the velocity tables actually used (`vel_run_NN.txt`), every
`strategy_C_summary.csv` / `strategy_C_log.csv` / `stratC_checkpoint.json`, and the
derived `postproc/` CSVs and figures. That is enough to reproduce or re-analyse
without re-running, for most purposes.

To restore the heavy data, either copy those folders across manually or re-run:
`Geo_Prep.dat` → `ANALYSIS_PART_I_MASON.dat` → `strategy_C_3dec_nodamp.py`.

---

## Running things

### Full IDA (inside 3DEC)

```
python-reset-state false
call 'strategy_C_3dec_nodamp.py'
```

Writes to `stratC_results_NODAMP_v6_NEW/`. Resumable via `stratC_checkpoint.json`.

Variants: `strategy_C_3dec_damp1p5.py` (1.5% Rayleigh),
`strategy_C_3dec_ratcheting.py` (3% Rayleigh),
`sinus_wave_IDA_strategyC_Maxwell.py` (original, 6% mass-proportional).
**See `CLAUDE.md` §3 before running the older two — they have Windows-only path
handling and restore a base save that is no longer produced.**

### Post-processing (outside 3DEC)

```bash
python postprocess_stratC.py stratC_results_NODAMP_v6_NEW \
    --label "No viscous damping (hysteretic only), asym. pulse" \
    --compare stratC_results_DAMP1p5 stratC_results_NODAMP_v6 \
    --compare-labels "1.5% Rayleigh" "NODAMP v6 (previous)"
```

Run this **first** — it creates `postproc/`, which every figure script writes into.

Then any of:

```bash
python exp_vs_sim_figs.py                       # figP1–figP5
python make_presentation_figs.py [SIM_DIR] [OUT_DIR]
python fig_sd_period.py                         # figP2_nodamp, figP2_ratcheting
python ch19_xval.py sim  &&  python ch19_xval.py finish
python profile_hysteresis.py sim  &&  python profile_hysteresis.py fig
```

Several of these are **staged**, caching intermediate state in `/tmp/*.json`. That
cache does not survive a reboot or a machine move — re-run the earlier stages before
`finish` / `finalize` / `fig`.

### Experimental data processing

```bash
python exp_tilt_from_raw.py  EXP_DIR --runs 1 24 --test 9
python exp_period_wrapper.py 9  EXP_DIR 1 25
python exp_period_wrapper.py 12 EXP_DIR 1 25
python exp_period_wrapper.py finalize
python process_tiltmeter.py  path/to/US1_Tilt_values.csv
```

`exp_tilt_from_raw.py` must process **Run 1 first** (it sets the baselines).

### 3DEC plot export (inside 3DEC)

```
call 'export_plots.py'
```

Edit the config block at the top for `SAV_DIR`, resolution, and which views to export.
