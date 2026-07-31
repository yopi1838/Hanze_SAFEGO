# Hanze_SAFEGO — project context

> This file is read automatically by Claude (Cowork / Claude Code) when this folder is
> opened. It is the project's working memory. Keep it current — if you change a
> convention, a path, or a driver, update the relevant section here in the same commit.

---

## 1. What this project is

Distinct Element Method (DEM) modelling in **Itasca 3DEC 9.1** of a two-leaf,
out-of-plane loaded **unreinforced masonry wall** with timber floor joists, tested on a
shake table. The numerical model is cross-validated against two physical specimens:

- **US-1 = Test 9**
- **US-2 = Test 12**

The central numerical experiment is an **adaptive sequential IDA, "Strategy C"**: a
fixed 25-run sequence of spectrum-matched sinusoidal pulses applied to the *same,
progressively damaged* model. After each run the fundamental period `T1` is
re-identified from that run's free-vibration ring-down, and the next run's excitation
is re-tuned to that elongated period. A second research thread ("ratcheting",
Modification A) makes the pulse **asymmetric** so the wall accumulates one-directional
lean, reproducing the residual tilt seen experimentally.

Research questions being answered:

1. Does an adaptive (period-tracking) IDA reproduce the experimental damage sequence
   better than a fixed-period IDA?
2. How much viscous damping is defensible? The comparison runs 0% / 1.5% / 3% / 6%.
3. Does an asymmetric pulse reproduce the observed residual tilt / ratcheting, and is
   the single free parameter `ASYM_K` transferable out of sample?

---

## 2. Folder layout

```
Hanze_SAFEGO/
├── CLAUDE.md, README.md          this file + new-machine setup
├── *.dat                         3DEC FISH / command files (model build + instruments)
├── Groningen.dec                 base brick layout, called by Geo_Prep.dat
├── strategy_C_3dec_*.py          the IDA drivers — run INSIDE 3DEC
├── ratcheting_pulse.py           asymmetric-pulse module (pure numpy, imported by drivers)
├── spectrum_{HU12,EC40,FR76}.csv target displacement spectra
├── libjmodelmason*.so            compiled `mason_v6` joint constitutive model (Linux)
├── postprocess_stratC.py         main postprocessor — run OUTSIDE 3DEC
├── exp_*.py, process_tiltmeter.py experimental-data processing
├── *_figs.py, fig_*.py, ch19_xval.py, profile_hysteresis.py   figure scripts
└── stratC_results_<VARIANT>/
    ├── strategy_C_summary.csv    ← per-run summary (TRACKED in git)
    ├── strategy_C_log.csv        ← run log (TRACKED)
    ├── stratC_checkpoint.json    ← resume state (TRACKED)
    ├── vel_run_NN.txt            ← generated velocity tables (TRACKED)
    ├── stratC_run_NN.sav         ← 3DEC state, 47–101 MB each (NOT in git)
    ├── RunNN_REC_sXpYY/*.csv     ← raw channel histories, ~2.5 MB each (NOT in git)
    ├── plots/*.png               ← 3DEC bitmap dumps (NOT in git, regenerable)
    └── postproc/                 ← derived CSVs + paper figures (TRACKED)
```

**What is deliberately not in git** (~10 GB): all `*.sav`, every `RunNN_*/` raw-CSV
folder, `plots/`, and the 3DEC project files `P.prj` / `P.temp` / `P.backup`.
See §9 for how to get them onto a new machine.

---

## 3. The four driver variants — READ THIS BEFORE RUNNING ANYTHING

All four are ~90% duplicated code with no shared module. **A fix must be applied to
each one separately.** `strategy_C_3dec_nodamp.py` is the only one that is both
current and cross-platform.

| Script | `OUT_DIR` | Damping actually applied | Base save restored | Portable paths? | Preflight? |
|---|---|---|---|---|---|
| `strategy_C_3dec_nodamp.py` **(current)** | `stratC_results_NODAMP_v6_NEW` | **none** — contact dissipation only | `Part_I_MASON_v6.sav` | yes | yes |
| `strategy_C_3dec_damp1p5.py` | `stratC_results_DAMP1p5` | `rayleigh 0.015 10.8696` (full) | `Part_I_MASON.sav` | yes | yes |
| `strategy_C_3dec_ratcheting.py` | `stratC_results_RATCHETING` | `rayleigh 0.03 10.8696` (full, hard-coded) | `Part_I_MASON.sav` | **no** (Windows-only) | no |
| `sinus_wave_IDA_strategyC_Maxwell.py` (ancestor) | `stratC_results_MP_NEW` | `rayleigh 0.06 10.8696 mass` | `Part_I_MASON.sav` | **no** | no |

### Traps to know about

- **`DAMP_RATIO` in `strategy_C_3dec_nodamp.py` is a lie if you change it.** Around
  lines 440–446 the Rayleigh command inside `if DAMP_RATIO > 0:` is commented out and
  the body is `pass`, while `preflight_checks()` still *prints* that damping is active.
  Setting it non-zero silently gives you an undamped run. Fix the branch, don't trust
  the printout.
- **`_damp1p5` and `_ratcheting` restore `Part_I_MASON.sav`**, which the current
  `ANALYSIS_PART_I_MASON.dat` no longer produces (it now writes `Part_I_MASON_v6`).
  Those two cannot be re-run from scratch without the old save file.
- **`cmd_path()` is broken in the two older drivers**: `path.replace("/", "\\")`.
  On Linux this corrupts paths and was the documented cause of period-ID failures
  (the driver could not locate the Channel-19 CSV). The newer two use
  `os.path.normpath(path).replace("\\", "/")`.
- **The folder name `stratC_results_RATCHETING` is misleading.** All three
  `strategy_C_3dec_*` drivers have `USE_RATCHETING = True`. What actually
  distinguishes that folder is **3% Rayleigh damping** and the non-`_v6` base save.
- **`stratC_results_NODAMP_v6` and `stratC_results_NODAMP_v6_NEW` are different runs.**
  The former is used as a `--compare` baseline; the latter is current.

---

## 4. Pipeline

```
Groningen.dec
   └─ call ─> Geo_Prep.dat ─────────────────────> model save 'GEO_WALL'
                 └─ restore ─> ANALYSIS_PART_I_MASON.dat ─> 'Part_I_MASON_v6'
                                  (gravity, precompression, jmodel mason_v6)
                                        │
        ┌───────────────────────────────┘   inside 3DEC:
        │                                     python-reset-state false
        │                                     call 'strategy_C_3dec_nodamp.py'
        ▼
  + ratcheting_pulse.py, spectrum_*.csv,
    instrument_history_new.dat, instrument_history_export_new.dat
        │
        │  per run: build vel_run_NN.txt → table import → apply velocity-z to
        │  groups 'S' and 'T_B' → solve pulse → remove BC → solve 0.5 s gap →
        │  model save → export histories → identify T_end from Ch19 ring-down →
        │  checkpoint → feed T_end into the next run
        ▼
  stratC_results_NODAMP_v6_NEW/
        │
        ▼  python postprocess_stratC.py            (outside 3DEC)
  postproc/postproc_all_channels.csv, postproc_summary.csv, fig1–fig5.png
        │
        ├─> exp_vs_sim_figs.py         → figP1…figP5
        ├─> make_presentation_figs.py  → figP1…figP3 into presentation_figs/
        ├─> fig_sd_period.py           → figP2_nodamp.png, figP2_ratcheting.png
        ├─> ch19_xval.py finish        → xval_ch19_peaks.csv, figP6, figP7
        └─> profile_hysteresis.py fig  → figP8, figP8b, figP9, figP10

Experimental branch (feeds the same postproc/ folder):
  Test{9,12}Run{N}.xlsx ─exp_tilt_from_raw.py──> exp_Test{9,12}_tilt.csv
  Test{9,12}Run{N}.xlsx ─exp_period_wrapper.py─> exp_Test{9,12}_period_psd.csv
  US1_Tilt_values.csv   ─process_tiltmeter.py──> exp_US{1,2}_tiltmeter.csv
  (exp_Test9_metrics.csv — NO PRODUCER SCRIPT EXISTS, see §8)

Side branch: export_plots.py (inside 3DEC) reads stratC_run_*.sav → plots/*.png
```

**Ordering constraints that actually bite:**

- `postprocess_stratC.py` must run first — it creates `postproc/`, which every
  downstream figure script writes into.
- `exp_period_wrapper.py TEST_NO …` must run for **both** tests before `finalize`,
  and `finalize` must precede `fig_sd_period.py`.
- `exp_tilt_from_raw.py` must process **Run 1 first** — it raises
  `SystemExit("Run 1 must be processed first (baselines).")` otherwise.
- Several scripts are **staged with a JSON cache in `/tmp`** (`/tmp/ch19_state.json`,
  `/tmp/prof_state.json`, `/tmp/period_state.json`, `/tmp/exp_tilt_state.json`).
  These do not survive a reboot or a new machine — the intermediate stages must be
  re-run before the `finish` / `finalize` / `fig` stage.

---

## 5. Conventions

### Run naming
`Run{NN:02d}_{RECORD}_s{scale, '.'→'p'}` — e.g. `Run01_HU12_s0p50`, `Run25_FR76_s2p00`.
Built by `"{:.2f}".format(scale).replace(".", "p")`; parsed back by
`RUN_RE = ^Run(\d+)_([A-Za-z0-9]+)_s(\d+)p(\d+)$` in `postprocess_stratC.py`.

### Records
Three ground-motion records used only as spectrum-lookup and asymmetry keys:

| Key | Role | `s` (dominant-peak sign) | `beta` (reverse/dominant) |
|---|---|---|---|
| `HU12` | workhorse, runs 1–21, scale 0.50→6.00 | +1 | 0.54 |
| `EC40` | interleaved low scales, runs 3/6/8/11 | +1 | 0.70 |
| `FR76` | final runs 22–25, scale 1.00→2.00 | −1 | 0.78 |

(`RECORD_ASYM`, `ratcheting_pulse.py` ~line 39. The acronyms are never expanded in code.)

### Scale factor `s`
Third element of each `PROTOCOL` tuple. It multiplies the record's spectral
displacement **at the wall's current period**:
`Sd_target = scale * interpolate_sd(record, T_current)`. Because the period elongates
with damage, the same nominal scale gives a different absolute Sd than fixed-T1
scaling — the ratio is logged as `amplification = Sd_target / Sd_fixedT1`.

### The 25-run protocol (identical in all drivers)
```
 1 HU12 0.50 |  2 HU12 0.75 |  3 EC40 0.20 |  4 HU12 1.00 |  5 HU12 1.25
 6 EC40 0.30 |  7 HU12 1.50 |  8 EC40 0.40 |  9 HU12 1.75 | 10 HU12 2.00
11 EC40 0.50 | 12 HU12 2.25 | 13 HU12 2.50 | 14 HU12 2.75 | 15 HU12 3.00
16 HU12 3.50 | 17 HU12 4.00 | 18 HU12 4.50 | 19 HU12 5.00 | 20 HU12 5.50
21 HU12 6.00 | 22 FR76 1.00 | 23 FR76 1.50 | 24 FR76 1.75 | 25 FR76 2.00
```

### History-CSV format
3DEC `history export … vs "time" file …` writes **2 header lines, then
whitespace-separated 2 columns**: col 0 = dynamic time-total (s), col 1 = the
quantity. Read with `np.genfromtxt(path, skip_header=2)`. Filename is
`<run_label>_<channel>.csv`.

### Units
SI throughout: m, kg/m³, Pa, m/s, m/s², s. **y is vertical** (`gravity 0 -9.81 0`),
z is out-of-plane, x is along the wall. Derived exports: `rel_disp_top_mm` in **mm**,
all `tilt_*` in **degrees**, `cstav` / `joist_s*_shear` / `total_shear` in **kN**.
Summary CSV: `Sd_*_mm` mm, `A_mps2` m/s², `PGA_g` g, `V_peak_mps` m/s, `T_*` s.
Experimental xlsx channels are in **metres** and get `*1000.0`.

### Instrument channels (block histories, positions in m)

| Channel | Quantity | Position (x, y, z) |
|---|---|---|
| `Channel_1_DispBot` | disp-z | 0.646, 0.66, 0.21 |
| `Channel_2_DispMid` | disp-z | 0.646, 1.26, 0.21 |
| `Channel_3_DispTopQLeft` | disp-z | 0.106, 2.06, 0.21 |
| `Channel_4_DispTopQRight` | disp-z | 1.186, 2.06, 0.21 |
| `Channel_5_DispTable` | disp-z | 0.646, 0.00, 0.21 |
| `Channel_8/9_LeftJoist…Slab1` | disp-z | 0.106, 0.56 / 0.66, 0.22 |
| `Channel_10/11_LeftJoist…Slab2` | disp-z | 0.106, 2.35 / 2.45, 0.22 |
| `Channel_12_AccTable` | **velocity-z** | 0.646, 0.00, 0.21 |
| `Channel_13/14_AccSlab1/2` | velocity-z | 0.646, 0.49 / 2.28, 0.73 |
| `Channel_15/16/17_AccBot/Mid/Top` | velocity-z | 0.646, 0.66 / 1.26 / 2.06, 0.21 |
| `Channel_19_DispTopQRight` | disp-z | 0.646, 2.06, 0.21 |
| `Channel_20_TopBeam` | disp-z | 0.646, 2.68, 0.21 |

Two gotchas: the `Acc*` channels are **velocity** histories (acceleration is derived
downstream), and **`Channel_19` is at the wall centreline despite the name
"TopQRight"**. Channel 19 is the period-identification channel.

FISH tilt chords, `atan(Δdz/dy)` in degrees:
`tilt_bot_seg` (0→0.66), `tilt_low_seg` (0.66→1.26), `tilt_up_seg` (1.26→2.06),
`tilt_beam_seg` (2.06→2.68), `tilt_full_wall` (0→2.06).
`rel_disp_top_mm = (dz_top − dz_table) * 1000`.

### Postproc metric definitions (`TAIL_FRAC = 0.05`)
- `peak` = max |x(t) − x(t₀)| within the run
- `residual` = mean of the last 5% of x(t) minus x at the **start of Run 1** — i.e.
  **cumulative across the sequence**, not per-run.

Key channels: `KEY_TILT = "tilt_full_wall"`, `KEY_DISP = "rel_disp_top_mm"`.

### Figure convention (all presentation scripts)
Red squares = simulation (NUM) · blue circles = US-1 (Test 9) · orange circles =
US-2 (Test 12) · purple diamonds = 3% Rayleigh / RATCHETING run.
Record colours: HU12 `tab:blue`, EC40 `tab:orange`, FR76 `tab:red`.

---

## 6. Model parameters

**Geometry.** Wall 1.292 (x) × 2.58 (y) × 0.21 (z) m, two leaves with a collar joint
at z = 0.105. Steel support `'S'` at y ∈ [−0.05, 0], top beam `'T_B'` at
y ∈ [2.58, 2.68]. Joist pockets at x ∈ [0.3230, 0.4038] and [0.8883, 0.9690], at
y ∈ [0.29318, 0.46909] and [2.05227, 2.22818]. `bLength = 1.292`.
`model large-strain off`.

**Densities.** Masonry 1885 kg/m³ · steel parts 7850 · timber joists `[5510*1.034]`
≈ 5697.

**Joint model `mason_v6`** (`block contact jmodel assign mason_v6`) — a **user-written
plugin**, not a built-in 3DEC jmodel; the compiled binary is `libjmodelmason*.so`
(see §7). Parameters from `fish define properties` in `ANALYSIS_PART_I_MASON.dat`:

```
kn = ks   = 2.5e9 Pa/m
tension   = 0.2e6 Pa
cohesion  = 0.3e6 Pa   (cohesion-residual = 0)
friction  = 36.9 deg   (friction-residual = same)
fc_comp   = 20e6 Pa    (comp-residual = 0.2*fc = 4e6)
Gt = (0.025*(2*tension*1e-6)^0.7)*1e3      → G_I
Gs = 10*Gt                                 → G_II
Gc = (15 + 0.43*fc_MPa − 0.036*fc_MPa^2)*1e3 → G_c
n  = 2.0 (peak_ratio) ; Css = 9.0 ; Cnn = 1.0 ; Cn = 0.0
```

Joist-to-support contacts use the built-in **Mohr** model with
`stiffness-normal [kn*0.1]`, `stiffness-shear 1`, `tension 1`, `cohesion 1`,
`friction 1` — deliberately near-frictionless sliding supports.

**Precompression.** `block face apply stress-yy [-10.35e3/(0.21*1.292)]` ≈ 38.2 kPa
on group `"Comp_apply"` (a 10.35 kN vertical load).

**Timestep.** Never set explicitly — 3DEC's automatic dynamic timestep is used.
`delta_t = 0.005` in the drivers is the **velocity-table sample interval**, not the
solve timestep.

**Excitation parameters** (identical in all drivers):
```
T1_init = 0.092 s   (→ 1/T1_init = 10.87 Hz, the Rayleigh centre frequency)
xi = 0.05           (only for the Newmark SDOF spectrum calibration)
delta_t = 0.005 s ; n_cycles = 1.5 ; tail_sec = 2.5 s ; inter_run_gap = 0.5 s
FFT_F_MIN/MAX = 2.0 / 50.0 Hz ; T_MIN/MAX_PHYS = 0.02 / 1.0 s
RMS_ALIVE_MM = 0.05 ; RMS_CHUNK_S = 0.15 ; MAX_RD_WINDOW = 0.50 s
phase_accel = pi/2  (declared, never used)
```

**Period identification.** Last free-vibration segment → RMS-based live window →
de-mean + Hanning → **normalised autocorrelation (primary)** with parabolic lag
refinement, **16× zero-padded rFFT (fallback)** with log-parabolic peak refinement.
If both fail the function returns `T1_init` — downstream scripts detect this as
`T_end_over_Tinit == 1.0` and flag it "period-ID fallback (excluded)".

**Dynamic BCs** (re-applied on every restore; not saved with the model):
```
block contact group 'Joist_S1_contact' range pos-y 0.25 0.3 pos-z 0.9 1.5
block contact group 'Joist_S2_contact' range pos-y 2.0 2.5 pos-z 0.9 1.5
block free velocity-z range group 'S'
block free rotation-y / rotation-z range group 'T_B'
block mech damp local 0.0 ; block mech damp global 0.0
```
Excitation: `block apply velocity-z 1.0 table '<run>' range group '<S|T_B>'` — **both
the base and the top beam are driven by the same table motion** (rigid-diaphragm-like
restraint) — then removed with `block gridpoint apply-remove velocity-z`.

**Ratcheting.** `USE_RATCHETING = True`, `ASYM_K = 1.0` in all three
`strategy_C_3dec_*` drivers; `False` reproduces the original symmetric cosine exactly.
`beta_eff()` is two-regime: β = 1 while `T_end/T1 ≤ 1.05`, ramping linearly to
`beta_record` by `T_end/T1 ≥ 1.20`. **`ASYM_K` is the only free parameter and is
documented as still needing out-of-sample validation.** The sign convention (`s` is in
the channel-12 table-accelerometer positive direction, and must be confirmed to map to
3DEC +z) also has an open check noted in `ratcheting_pulse.py`.

---

## 7. Environment and dependencies

- **3DEC 9.1** (`3dec910` appears in a hard-coded path). Syntax throughout is
  3DEC 7+/9 style.
- **`mason_v6` jmodel plugin** — `libjmodelmason009.so` / `libjmodelmasonv6009.so`
  are in this folder. Note these are **Linux** shared objects. Nothing in the code
  loads them explicitly; 3DEC must find them on its plugin search path, or a startup
  file not in this repo loads them. **On a Windows machine you need the corresponding
  `.dll` instead** — `block contact jmodel assign mason_v6` will fail without it.
- **Python** — the drivers run in 3DEC's embedded interpreter, so they are written in
  a Python-2/3-compatible style (`.format()`, no f-strings, ASCII coding headers).
  Preserve that style when editing `strategy_C_3dec_*.py`, `export_plots.py`,
  `sinus_wave_IDA_*.py`.
- **Packages**: `numpy` (all), `matplotlib` (always `matplotlib.use("Agg")`),
  `scipy.signal` (`exp_period_evolution.py`), `openpyxl` (all xlsx readers,
  `read_only=True, data_only=True`), `pandas` (`process_tiltmeter.py` only),
  `itasca` (3DEC-embedded only). See `requirements.txt`.

---

## 8. Known issues, WIP, and things that will break on a new machine

### Hard-coded absolute paths — fix these first

| File | Path | Overridable? |
|---|---|---|
| `exp_tilt_from_raw.py` (~L39) | `C:\Users\yopi1\Documents\Itasca\3dec910\Hanze\Wall_Floor_Interaction\EXP_DATA` | yes, CLI arg |
| `process_tiltmeter.py` (~L22) | `C:\Users\yopi1\Documents\000.DISSERTATION_FILES\00.Postprocess\dynamic\Hanze\EXP_DATA\US1_Tilt_values.csv` | yes, `sys.argv[1]` |
| `xval_exp_vs_sim.py` (~L6) | `/sessions/zen-brave-heisenberg/mnt/Hanze_SAFEGO/stratC_results_NODAMP_v6_NEW` | **no** |
| `xval_exp_vs_sim.py` (~L14) | `/tmp/exp_metrics.csv` | **no** |
| `profile_hysteresis.py` (~L28) | `/sessions/epic-stoic-galileo/mnt/Hanze--EXP_DATA/processed_globalzero` | **no** |
| `profile_hysteresis.py` (~L189) | `/sessions/epic-stoic-galileo/mnt/EXP_DATA` | **no** |
| `exp_period_evolution.py` (L20–21) | `dynamic\Hanze\EXP_DATA_2`, `dynamic\results\Hanze` (Windows-relative, currently unused) | n/a |

The `/sessions/...` paths belong to two dead cloud sandbox sessions and will never
resolve anywhere. Many other scripts hard-wire the **relative** folder name
`stratC_results_NODAMP_v6_NEW` next to the script (`exp_vs_sim_figs.py`,
`ch19_xval.py`, `exp_period_wrapper.py`, `exp_tilt_from_raw.py`,
`process_tiltmeter.py`, `fig_sd_period.py`, `profile_hysteresis.py`).

The 3DEC drivers resolve **everything** relative to 3DEC's current working directory,
so 3DEC must be launched with this folder as cwd.

### Data not in this repo

| File | Needed by | Where it comes from |
|---|---|---|
| `Test{9,12}Run{1..25}.xlsx` | `exp_tilt_from_raw.py`, `ch19_xval.py`, `profile_hysteresis.py`, `exp_period_evolution.py` | raw shake-table exports (`EXP_DATA` folder) |
| `Test9Run{N}_processed_globalzero.xlsx` | `profile_hysteresis.py` | processed exp data (`U_avg`, `base_shear_kN`) |
| `Test{9,12}_Info.xlsx` | channel maps | `EXP_DATA` |
| `US1_Tilt_values.csv`, `US2_…` | `process_tiltmeter.py` | 5 s tiltmeter logger; `;`-separated, decimal comma, millidegrees |
| `exp_Test9_metrics.csv` | `exp_vs_sim_figs.py`, `make_presentation_figs.py` | **no producer script exists.** Needs columns `run, peak_rel_mm, resid_rel_mm, raw_base4_mm, peak_tilt_deg, resid_tilt_deg`. If a copy survives in `postproc/` it is tracked; otherwise it must be rebuilt. |
| `ratcheting_pulse_spec.md` | referenced in `ratcheting_pulse.py` docstring | never written |
| `compare_tilt_robust.py` | referenced in `process_tiltmeter.py` | superseded; used the **wrong tilt axis** (`tvalx`; correct is `tvaly`) |
| `Part_I.sav`, `V_L.tab`, `V_.tab` | `Dynamic_Analysis.dat` | legacy, not needed |

### Open / half-finished

- **Four near-duplicate drivers** with no shared module (§3). Consolidating them into
  one parameterised driver is the obvious refactor and has not been done.
- **`exp_period_evolution.py` has no entry point.** `exp_period_wrapper.py` refers to
  `exp_period_evolution.run_all`, which does not exist. The module is a library lifted
  out of a larger script.
- **Three overlapping "presentation figure" scripts** producing same-named figures in
  different folders with different fonts and slightly different residual/tilt
  definitions: `exp_vs_sim_figs.py` (→ `postproc/`, DejaVu Serif),
  `make_presentation_figs.py` (→ `presentation_figs/`, STIXGeneral), `fig_sd_period.py`
  (→ `figP2_nodamp.png` / `figP2_ratcheting.png`). **Not documented which is canonical.**
- **`xval_exp_vs_sim.py` is throwaway** — no argparse, no `__main__` guard (it runs on
  import), two dead absolute paths, a duplicated record map, a figure title with a
  pre-baked conclusion, and it `KeyError`s on `sim_resid[25]` if run 25 is missing.
- **`tilt_angles` (FISH index 10)** is registered but deliberately never exported
  (duplicate of `tilt_bot_seg`). `instrument_history_export_new.dat` skips index 10.
- **`ncstav` was dropped** — still defined and historied in
  `ANALYSIS_PART_I_MASON.dat`, but the rewritten `cstav` in
  `instrument_history_new.dat` no longer computes it. The two `cstav` definitions
  differ substantively (the `.dat` version area-weights and divides by `jarea`).
- **`FISH_HISTORIES` index maps disagree** between the newer (16-entry) and older
  (13-entry) drivers. Harmless in Python — the drivers read CSVs by filename — but
  `instrument_history_export_new.dat` **does** export by numeric index, so pairing it
  with the older instrument file would silently mislabel exports.
- **Dead variables**: `phase_accel`, `bLength` (in Python), `check_run_complete()`
  (defined, never called, all four drivers), `INSTRUMENT_CHANNELS`, `HAS_PLOTLIB`,
  `PLOT_DEF_FILE`.
- **US-2 (Test 12)** experimental comparison is flagged "will be added when the data
  is available" in `exp_vs_sim_figs.py`.
- **`Dynamic_Analysis.dat`** is legacy (restores `Part_I`, not `Part_I_MASON_v6`);
  superseded by the Strategy-C drivers.

### A methodological note worth not re-deriving
`exp_tilt_from_raw.py` documents that the wall displacement sensors are mounted on the
shake-table frame, so they already read **relative** displacement — `d_table` is
therefore treated as zero. Verified: in Run 3 the table swings ±22 mm while the wall
channels move < 0.2 mm. `d_top` = mean of channels 3 and 4.
Separately, `process_tiltmeter.py` records that the out-of-plane tiltmeter axis is
**Y** (`sensor1 tvaly`), which reproduces the published US-1 curve (2.77° at Run 21,
6.03° at Run 24).

---

## 9. Getting the heavy data onto a new machine

Git carries the code, the model definition, the run summaries and the derived
`postproc/` outputs — about 30 MB. It does **not** carry ~10 GB of `.sav` states and
raw channel CSVs. Options, in order of preference:

1. **Regenerate.** `Geo_Prep.dat` → `ANALYSIS_PART_I_MASON.dat` → the driver.
   Everything is deterministic; this is the honest route, but it is a full re-run.
2. **Copy the excluded folders manually** (external drive / OneDrive / rsync):
   `*.sav`, `stratC_results_*/Run*/`, `stratC_results_*/plots/`, `P.prj`.
3. **Postprocess only.** If you just need figures, `postproc/` is already in git —
   the figure scripts can run without the raw CSVs, except `ch19_xval.py` and
   `profile_hysteresis.py`, which read raw channel CSVs directly.

---

## 10. Working agreements for Claude in this repo

- **Preserve the 3DEC-embedded Python style** in `strategy_C_3dec_*.py`,
  `sinus_wave_IDA_*.py`, `export_plots.py`: `.format()` not f-strings, no walrus, no
  dataclasses. They run in 3DEC's interpreter.
- **Never silently "fix" a duplicate driver in isolation.** If a bug exists in
  `_nodamp` it almost certainly exists in the other three — say so, and ask which to
  patch.
- **Do not commit `.sav`, raw `RunNN_*/` CSVs, or `plots/`.** `.gitignore` handles it;
  do not `git add -f` around it.
- **Line endings are frozen** by `.gitattributes` (`* -text`). The `.dat` files are
  CRLF and 3DEC is fine with that — do not normalise them.
- When changing a convention, a channel definition, or a path, **update this file in
  the same commit.**
