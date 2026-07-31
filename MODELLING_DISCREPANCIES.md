# Model vs experiment: where the discrepancies actually come from

**Scope.** Strategy-C adaptive IDA (`stratC_results_NODAMP_v6_NEW`, `_DAMP1p5`,
`_RATCHETING`) against shake-table specimens US-1 (Test 9, runs 1–24) and US-2
(Test 12, runs 1–25). All numbers below are computed from the tracked CSVs in
`stratC_results_NODAMP_v6_NEW/postproc/` and the raw channel exports; nothing is
estimated.

---

## Summary

The two symptoms in the title — "collapse is under-predicted" and "the model ratchets
but the experiment doesn't" — are **not two facets of one modelling error**. They
decompose into four distinct problems, and only two of them are constitutive:

| # | Problem | Kind | Severity |
|---|---|---|---|
| 1 | "Tilt" names two incommensurable quantities; the model outputs only one of them | Measurement / definition | **Invalidates the comparison** |
| 2 | Period identification failed at run 23, so run 24 was excited at the undamaged period | Algorithmic bug | **Explains the run-24 under-prediction entirely** |
| 3 | The model forms a full-height rocking mechanism at run 14; the specimens never do | Constitutive / mechanism | Real, and the important one |
| 4 | No single damping value works across the amplitude range | Constitutive / energy | Real |

Problems 1 and 2 have to be cleared before problem 3 can be assessed honestly,
because at present a good deal of the apparent disagreement is bookkeeping.

---

## 1. The tilt comparison is not like-for-like

There are two experimental "tilt" series and they disagree with each other by roughly
two orders of magnitude **and in sign**:

| Run | `exp_Test9_tilt.csv` residual chord tilt | `exp_US1_tiltmeter.csv` residual |
|---|---|---|
| 10 | −0.0106° | +0.202° |
| 15 | −0.0122° | +0.566° |
| 21 | −0.0294° | +2.768° |
| 24 | +0.0482° | +6.034° |

They are not measuring the same thing.

**The displacement-derived series is a chord angle.** `atan(peak_mm / 2060)` reproduces
`peak_tilt_full_wall` to three decimals across every run (run 21: 0.4680° computed vs
0.4552° reported). It is a global measure of how far the top of the wall has moved
relative to the table.

**The inclinometer is not a chord angle, and its peak values are not usable.**
The decisive evidence is in the experimental data alone, no model required:

| Run | Wall movement measured | Inclinometer peak | Top displacement that angle *would* require |
|---|---|---|---|
| 1 | 0.327 mm | 9.92° | 360 mm |
| 3 | **0.018 mm** | **26.77°** | 985 mm |
| 10 | 2.008 mm | 24.32° | 931 mm |
| 21 | 16.825 mm | 23.73° | 905 mm |

At run 3 the wall moved eighteen microns while the device reported 26.8°, which would
need 985 mm of top displacement. Whatever that number is, it is not a rigid-body
rotation of the wall.

### What the raw logger file establishes about the instrument

`EXP_DATA/US1_Tilt_values.csv` records, every 5 s, a tilt triple (value/min/max, in
millidegrees) on two axes **plus two full 3-axis accelerometer packs** (value/min/max,
in raw counts). That structure allows the instrument's behaviour to be checked directly:

- **It is gravity-referenced.** Over the first 50 resting windows the acceleration vector
  has constant magnitude 16750 counts (sd = 3), of which 16745 lies on the z axis. The
  device is tracking the gravity vector and deriving tilt from its direction. A
  gravity-referenced device necessarily reads horizontal acceleration as apparent tilt
  during shaking — that much is unavoidable physics, not an assumption about the sensor.
- **The accelerometer is a 16-bit ±2 g digital part.** `amaxy1` and `amaxz1` both rail at
  exactly **32760 counts = 1.96 g** (2¹⁵ − 8), and the measured scale of ~16750 counts/g
  is close to the nominal 2¹⁴ = 16384 counts/g. This is characteristic of a MEMS
  accelerometer, though the data does not identify the part or the sensing technology.
- **Both channels saturate during the strong runs.** The tilt channel rails symmetrically
  at **±25.75°** (`tmaxy` max +25754 mdeg, `tminy` min −25752 mdeg; 34 of 1016 windows sit
  within 1% of the rail), and the accelerometers rail at ±1.96 g. In the strong runs the
  device is simply at its limits.

`tmaxy` correlates with peak lateral acceleration at r = 0.69 across all 1016 windows —
consistent with an acceleration-driven reading, but **a simple `atan(a/|a|)` law does not
reproduce the magnitudes** (it predicts 5.7°–62.7° in windows where the device reported
25.1°–25.8°), largely because both channels are clipping. The mechanism is established;
a quantitative transfer function is not.

**What the device is attached to is not documented in the data available here.**
`Test9_Info.xlsx` and `Test12_Info.xlsx` list only the 17 displacement/force/acceleration
channels — the tiltmeter is not in either. So while the kinematics rule out its reading
being a global chord rotation (2.77° over 2.06 m would be 99 mm of permanent top
displacement, against −1.18 mm on the absolute transducer), **what local rotation it does
represent cannot be determined without the mounting record.** That record is the single
most valuable missing piece for this comparison.

`process_tiltmeter.py` is already right to extract settled quiet-window values rather
than peaks.

**Consequence.** The model outputs only chord tilt (`tilt_full_wall` = `atan((Ch19 −
Ch5)/2.06)`, verified to 1e-6°). Plotted against the inclinometer it under-reads by ~7x;
plotted against the displacement-derived tilt it over-reads by ~12x. Both comparisons
have been made in the figure scripts. Neither is meaningful until the model emits a
local-rotation channel — see §5.

### A definitional mismatch that turned out not to matter

The FISH takes `d_top` from the **centreline** gridpoint (x = 0.646); the experiment
averages the two top quarter points, Ch3 (x = 0.106) and Ch4 (x = 1.186). Recomputing
from the raw channels:

| Run | centreline | quarter-point mean | centreline bulge | Ch3−Ch4 spread |
|---|---|---|---|---|
| 15 | 12.719 mm | 12.713 mm | 0.072 mm | 0.250 mm |
| 21 | 32.974 mm | 32.958 mm | 0.024 mm | 0.097 mm |

A 0.2% effect. **The wall spans one-way** — there is no meaningful arching or twist, so
this mismatch is not a source of error. It is worth fixing for defensibility, not for
accuracy. Reported here so it is not mistaken for a finding.

One mismatch that *does* bite: the FISH hard-codes the bottom sensor at 0.66 m, but on
US-2 that sensor sits at 0.60 m. Every model-vs-US-2 segment-tilt comparison has been
using the wrong chord length.

---

## 2. Run 24 was excited at the wrong period — this is the "under-predicted collapse"

The driver identifies `T_end` from each run's ring-down and tunes the next run to it.
When identification fails it silently returns `T1_init`, which shows up as
`T_end_over_Tinit == 1.0000` exactly. That happened on **runs 23 and 25**.

Run 23's failure set run 24's excitation period to 0.092 s — the period of an
*undamaged* wall — while the specimen model had softened to about 3.9× that:

| Run | Excitation T/T₀ | Sd target | Amplification | Peak response |
|---|---|---|---|---|
| 22 | 4.15 | 27.62 mm | 21.3 | 33.9 mm |
| 23 | 3.92 | 38.06 mm | 19.5 | 46.7 mm |
| **24** | **1.00** | **2.27 mm** | **1.0** | **9.2 mm** |
| 25 | 4.35 | 57.58 mm | 22.2 | 194.4 mm |

The target spectral displacement dropped **17-fold** between runs 23 and 24 purely
because of the failed identification. US-1's run 24 — the collapse run, where the
inclinometer jumps from 2.4° to 6.0° — is being compared against a model run that was
given almost no energy.

**So the run-24 under-prediction is not a constitutive failure.** It is an artefact, and
any conclusion drawn from that data point should be withdrawn. Note also that the model
does not fail to reach collapse-level response in general: run 25 produced 194 mm.

**Fix before re-running:** make period-ID failure loud rather than silent. Options, in
order of preference — (a) abort the sequence and require intervention; (b) carry forward
the *previous successful* `T_end` instead of `T1_init`, which is far less wrong;
(c) at minimum, write a `period_id_ok` flag column into `strategy_C_summary.csv` so
affected runs can be excluded automatically. Currently the only trace is the exact
value 1.0000, which is easy to miss and impossible to distinguish from a genuinely
unsoftened run.

---

## 3. The model forms a rocking mechanism at run 14; the specimens never do

This is the real modelling problem, and it is visible in the period evolution:

| Run | Model T/T₀ | US-1 T/T₀ | US-2 T/T₀ |
|---|---|---|---|
| 12 | 1.22 | 1.09 | 0.96 |
| 13 | 1.27 | 1.09 | 0.99 |
| **14** | **3.16** | 1.09 | 1.01 |
| 15 | 3.90 | 1.09 | 1.00 |
| 21 | 4.15 | 1.12 | 1.01 |
| max over sequence | **4.37** | **1.22** | **1.23** |

The model's stiffness collapses by a factor of ~10 (period ×3.2) **in a single run**,
and stays there. The specimens soften smoothly and never exceed ~23% period elongation
over the entire 24–25 run sequence, right through to physical collapse.

This is the mechanism behind the ratcheting discrepancy. Compare like with like — model
chord tilt against the experiment's absolute transducer drift:

| Quantity at run 21 | Value |
|---|---|
| Experiment, absolute transducer drift | −1.18 mm (−0.033°) |
| Model, cumulative residual chord tilt | +14.24 mm (+0.396°) |
| **Model over-predicts permanent lean by** | **≈12×** |

Across all 24 runs the specimen's absolute drift is **+0.21 mm**. The specimen genuinely
does not accumulate global lean — while the inclinometer shows it *did* accumulate 6° of
local unit rotation.

**Interpretation.** In the specimens, damage distributed itself as many small
bed-joint rotations that largely cancelled — units rotate, the wall stays plumb, and
global stiffness barely changes. In the model, damage localises into a single one-sided
full-height rocking mechanism: a few joints open, the period quadruples, and every
subsequent asymmetric pulse pushes the same mechanism further in the same direction.

The ratcheting is therefore **a consequence of premature mechanism formation, not an
independent defect**. `ASYM_K` is not the thing to tune first — an asymmetric pulse
applied to a wall that has already become a one-sided rocking block will ratchet
whatever value it takes. Candidates for the actual cause, in the order I would test them:

1. **Joint tensile/cohesive softening is too brittle.** `G_I` is derived as
   `0.025·(2·f_t)^0.7` ≈ 5.4 J/m² for `f_t` = 0.2 MPa. Once a course opens it cannot
   recover, so damage cannot redistribute. A more ductile mode-I response, or the
   `mason_v7` healing behaviour, should let cracking spread rather than localise.
2. **Insufficient joints available to crack.** If the collar joint and bed-joint
   discretisation give only a few candidate planes, localisation is forced by the mesh
   rather than by the mechanics. Worth checking how many distinct bed joints actually
   opened in `stratC_run_14.sav` versus run 13.
3. **The boundary condition drives both ends in phase.** `block apply velocity-z` is
   applied to group `'S'` *and* group `'T_B'` with the same table motion, which imposes a
   near-rigid diaphragm at the top. Whether the physical top restraint was that stiff
   determines whether a full-height single-curvature mechanism is even admissible.
4. **No re-contact stiffness recovery.** If closed joints do not recover normal
   stiffness on impact, each rocking cycle loses a little more, biasing accumulation.

Item 3 is the cheapest to check and the most likely to change the mechanism
qualitatively.

---

## 4. No damping value works across the amplitude range

Peak response ratio (model / US-1), run 24 excluded as invalid:

| Regime | NODAMP | Rayleigh 1.5% | Rayleigh 3% |
|---|---|---|---|
| runs 1–13 (pre-softening) | **2.10×** | 1.27× | 1.61× |
| runs 14–21 (rocking) | 2.05× | **0.64×** | 0.67× |
| runs 22–25 (FR76) | 6.40× | 0.52× | 0.89× |

Undamped over-predicts everywhere. Add 1.5% Rayleigh and the pre-softening runs come
into line (1.27×), but the rocking regime is then **under**-predicted by a third.

That crossover is diagnostic. It says the joint model dissipates **too little energy
before the mechanism forms** — which is why viscous damping appears to help early — and
that the viscous damping then does the wrong job later, because Rayleigh damping
proportional to velocity heavily penalises the slow, large-amplitude rocking that
dominates runs 14+. Tuning a single global damping ratio cannot fix both ends; the
dissipation needs to come from the joint hysteresis, which returns to §3 item 1.

The `strategy_C_3dec_nodamp.py` trap is worth restating here: `DAMP_RATIO` is
non-functional in that file (the Rayleigh command is commented out inside
`if DAMP_RATIO > 0:` with `pass` as the body) while `preflight_checks()` still prints
that damping is active. Any damping study run from that script will silently produce
undamped results.

---

## 5. The tilt reimplementation

Two new files, both additive — nothing existing is renumbered or overwritten.

### `instrument_tilt_v2.dat` — call after `instrument_history_new.dat`

- `d_top` is the **mean of Ch3 and Ch4**, matching the experiment's definition rather
  than the centreline gridpoint. (Numerically a 0.2% change, per §1 — done for
  defensibility.)
- The bottom sensor elevation is a variable `[bot_sensor_y]`: **0.66 for US-1, 0.60 for
  US-2**. This fixes the wrong-chord comparison against US-2.
- Adds **`tilt_local_incl`** — rotation of a chord *inside a single masonry unit*. This
  is the missing counterpart to the bonded inclinometer, and is what should be compared
  against `residual_deg` in `exp_US{1,2}_tiltmeter.csv`.
- Adds `bulge_top_mm` (centreline minus quarter-point mean) as a running check that the
  wall is still spanning one-way. If it grows, the single-chord tilt definition has
  stopped being adequate for either model or test.

> **`tilt_local_incl` cannot be validated until the tiltmeter's mounting is known.**
> `[incl_y_lo]` / `[incl_y_hi]` are placeholders (1.95–2.06 m) and the mounting is *not*
> recorded in `Test9_Info.xlsx` / `Test12_Info.xlsx` — the tiltmeter is absent from both.
> Recover it from the test report, the instrumentation photographs or the lab, and set
> these before using the channel. A wrong height gives a plausible but wrong number,
> which is worse than no number.

> `instrument_history_export_new.dat` exports FISH histories **by numeric index**. The
> new channels are registered at 17–24; that export file must be extended to match or
> they will simply never be written out.

### `tilt_experiment_match.py` — works on existing runs, no re-run needed

Recomputes chord tilt using the experimental definition directly from already-exported
channel CSVs, and emulates the inclinometer as
`atan(a_z/g + tan θ_true)` with a band-limit to represent the device's finite response.

Verified against the existing v6 data: run 15 → 12.713 mm / 0.3536°, run 21 →
32.958 mm / 0.9165°, both matching hand calculation, and the exported `tilt_full_wall`
reproduces `atan((Ch19−Ch5)/2.06)` to 1e-6°.

```bash
python tilt_experiment_match.py stratC_results_NODAMP_v6_NEW --test 9
python tilt_experiment_match.py stratC_results_NODAMP_v6_NEW --test 12
```

**Honest limitation on the emulation.** With no `tilt_local_incl` channel it reproduces
only the acceleration artefact, and the emulated peak is sensitive to the assumed device
bandwidth — 31°/47°/52° at 10/25/50 Hz for run 15, against 24.3° measured. It
demonstrates the mechanism and the order of magnitude; it is **not** a calibrated
prediction, and shouldn't be presented as one. The settled value is far less sensitive
and is the one worth comparing, once the model actually emits a local rotation.

---

## Recommended order of work

1. **Fix the period-ID fallback** (§2) before any further runs. Everything downstream of
   a silent fallback is uninterpretable, and it hit the two most important runs.
2. **Re-plot every tilt comparison** with the right pairing (§1): chord vs chord,
   local rotation vs inclinometer. Some of the current disagreement will disappear;
   what remains is the real signal.
3. **Recover the tiltmeter's mounting location from the test records** — it is not in
   the data here — then set `[incl_y_lo]`/`[incl_y_hi]` and get `tilt_local_incl` out of
   the v7 run. Without it there is no valid model counterpart to the headline
   experimental result. Also treat every inclinometer reading above ~25° as clipped.
4. **Then** investigate premature mechanism formation (§3) — starting with the top
   boundary condition, which is cheap to test and could change the mechanism outright.
5. Leave `ASYM_K` alone until 1–4 are done. It is currently absorbing the error from
   everything above, and any value fitted now will be fitted to artefacts.
