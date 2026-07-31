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
| 0 | The tiltmeter's tilt channel disagrees with its own two accelerometers by 20–30× | Instrumentation | **Decides whether problem 3 exists** |
| 1 | "Tilt" names two incommensurable quantities; the model outputs only one of them | Measurement / definition | **Invalidates the comparison** |
| 2 | Period identification failed at run 23, so run 24 was excited at the undamaged period | Algorithmic bug | **Explains the run-24 under-prediction entirely** |
| 3 | The model forms a full-height rocking mechanism at run 14; the specimens never do | Constitutive / mechanism | Real, and the important one |
| 4 | No single damping value works across the amplitude range | Constitutive / energy | Real |

Problems 0, 1 and 2 have to be cleared before problem 3 can be assessed honestly,
because at present a good deal of the apparent disagreement is bookkeeping. Problem 0 in
particular may remove most of the ratcheting discrepancy on its own.

---

## 0. Corrections from the source paper

Moshfeghi, Smyrou, Arslan & Bal, *Structures* **66** (2024) 106815 (the SafeGo shake-table
paper) resolves several things this document originally inferred from data alone. Where
they conflict, the paper wins.

| Earlier claim here | Corrected |
|---|---|
| US-1 "collapsed" at run 24 at ~29.7 mm | **Neither specimen collapsed.** The sequence stopped at the shake table's capacity (PTA 0.78 g). Cracks were "visible… with almost no change in length and width" at the end. 29.7 mm is simply the largest displacement the table could impose. |
| The tiltmeter's mounting is undocumented | **Table 3: X = 0, Z = 2.43 m**, "Separately Retrieved from an Online Portal" — i.e. a vendor-processed value, not a raw channel. |
| The tilt-channel / accelerometer dispute is unresolved | **Resolved in favour of the tilt channel** — see below. |
| The wall is two-leaf with a collar joint | **Solid one-brick wall**, alternating stretcher and header courses; the paper states "the two layers in the solid walls are fully engaged due to the bricklaying technique". No collar joint exists. |

**Why the tilt channel wins.** The deformed shape (Fig. 15) has US-1 reaching ~28–30 mm at
2.06 m with the top beam held at zero at 2.85 m, so the chord rotation over that upper
segment is `atan(28/790)` ≈ **2.0°** at peak — the same order as the reported residual
tilt, and an order of magnitude above the 0.27° the accelerometers give. The authors
independently read the 2.2–6° residual tilt as consistent with those deformed shapes, and
the model's own `tilt_beam_seg` gives −1.51° at run 21 and −2.37° at run 23. Three
independent lines land on the tilt channel. Since the tilt value comes from the vendor's
online portal rather than the raw device, there is no reason to expect the accelerometer
words in the same export to be DC-preserving — auto-zeroing or high-pass filtering there
would produce exactly the under-reading measured.

The §1 analysis below therefore stands in its **category** conclusion — chord tilt and the
tiltmeter measure different things and the model was compared against the wrong one — but
the accelerometer cross-check should now be read as *evidence about the export pipeline*,
not as a challenge to the published curve. The peak columns remain unusable regardless:
they are acceleration-dominated and clipped at ±25.75°, and the paper confirms the
instrument was only ever read pre- and post-run.

**Two further findings, both directly useful:**

- **Measured period elongation: +17.5% (US-1, 0.091→0.107 s) and +13.9% (US-2,
  0.087→0.099 s)** — confirming the values derived here independently, and making the
  model's +337% unambiguously the central problem.
- **The asymmetry has a physical cause the model already contains.** The paper attributes
  it to "the timber floors applying a moment towards one side of the wall section", and
  measures **8–12 mm of cumulative joist slip**. The bow-tie specimen, which eliminated
  that slip, reversed the sign of its residual tilt. `ASYM_K` imposes an asymmetric input
  pulse instead — a different mechanism, and possibly a double count.

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

**Mounting.** The tiltmeter appears in neither `Test9_Info.xlsx` nor
`Test12_Info.xlsx` (both list only the 17 displacement/force/acceleration channels). The
only description is the paper's one sentence: *pre- and post-run tilt measurements on top
of the wall, for in-plane and out-of-plane residual movements.* That settles three things:
the device sits at the **top of the wall** (so `[incl_y_lo]`/`[incl_y_hi]` ≈ 1.95–2.06 m
is the right region); its two axes are **in-plane (`tvalx`) and OOP (`tvaly`)**, matching
`process_tiltmeter.py`'s choice of `tvaly`; and it was only ever intended to be read
**statically, before and after each run**. The peak columns were never a measurement —
which independently confirms that they should not be plotted against anything.

### The tilt channel disagrees with its own accelerometers by a factor of 20–30

This is the most consequential thing in this document, and it needs resolving before the
ratcheting comparison can be trusted at all.

At rest the accelerometer vector *is* the gravity direction, so `atan2(a_y, a_z)` is an
independent measurement of exactly the rotation the tilt channel reports. Measured well
after all shaking has stopped, with both readings stable (tilt sd = 1 mdeg over the final
7 logger windows):

| Permanent rotation, final state | US-1 | US-2 |
|---|---|---|
| Tilt channel `tvaly` — **the published curve** | **+6.038°** | **+2.367°** |
| Accelerometer pack 1 | +0.275° | +0.074° |
| Accelerometer pack 2 | +0.238° | +0.073° |
| Disagreement | **22×** | **32×** |

The two accelerometer packs are independent chips, and they agree with each other to
r = 0.995 across all 24 run-by-run settled values. Where two independent sensors agree
with each other and not with a third, the third is the outlier.

The displacement transducers agree with the accelerometers, not with the tilt channel:

| Final-run residual | US-1 | US-2 |
|---|---|---|
| Upper-segment chord (1.26→2.06 m), nearest the sensor | +0.069° | −0.035° |
| Full-wall chord | +0.048° | +0.038° |
| Tiltmeter accelerometers | +0.275° | +0.074° |
| Tiltmeter tilt channel | +6.038° | +2.367° |

On US-2 the accelerometers (+0.074°) and the chord tilts (+0.038 to +0.11°) are the same
order of magnitude. On US-1 the accelerometers read ~4× the chord, which is what you would
expect for a genuinely *local* rotation at the top of a damaged wall. The tilt channel is
90–125× the chord on US-1 and 60× on US-2.

It is **not a constant scale error.** In the well-resolved region the ratio grows
monotonically with amplitude — US-1: 5.1× at run 10, 10.3× at run 15, 19.6× at run 21,
22.1× at run 24. A unit or calibration mistake would give a fixed factor. A log-log fit
gives `tilt_channel ≈ (accelerometer angle)^k` with **k = 2.04 for US-1 and 1.61 for
US-2** — superlinear on both specimens, but not a single clean power law.

Three caveats on this, from a closer pass over both logger files:

- **Accelerometer resolution is ~0.0034°** (1 count at 16750 counts/g). For roughly the
  first ten runs the accel-derived angle is only a handful of counts, so the run-by-run
  ratios there are quantisation noise and should be ignored. The growth quoted above is
  taken from the resolved region only, and is real.
- **US-2's final value is less trustworthy than US-1's.** Its log ends just 45 s after
  the last run, and the wall was still moving 25 s before the end (swing 2504 counts at
  t = 5480 s, `tvaly` swinging between −1669 and +3696 mdeg). Only the last four windows
  are settled, at ~2210 mdeg. US-1 by contrast has five quiet minutes and a tilt sd of
  1 mdeg. Treat US-2's 32× as indicative, US-1's 22× as solid.
- **The in-plane channel moved too**: `tvalx` net +1.120° on US-2 and −0.359° on US-1.
  Worth explaining under any hypothesis — a purely OOP wall rotation should leave it
  near zero.

### …but the instrumentation drawing and the model both cut the other way

The instrumentation figure places sensor 18 (tiltmeter) **at the very top of the wall**,
on a bracket above the slab-2 connection — higher than the top quarter level, around
2.45–2.58 m. That is a single small housing, so its tilt element and its accelerometers
undergo the *same* rotation; the disagreement between them cannot be a location effect
and must be instrumental.

It also identifies the right model counterpart, and the model already has one:
`tilt_beam_seg`, the 2.06 → 2.68 m chord. Comparing cumulative residuals:

| Run | Model `tilt_beam_seg` (top region) | US-1 tilt channel | Model `tilt_full_wall` | US-1 accelerometers |
|---|---|---|---|---|
| 15 | −0.98° | +0.57° | +0.26° | +0.056° |
| 21 | −1.51° | +2.77° | +0.40° | +0.141° |
| 23 | −2.37° | +2.42° | +0.62° | +0.145° |

The model's top-of-wall rotation is **within a factor of ~2 of the tilt channel, and at
run 23 essentially equal to it** — while being 10–16× the accelerometer values. (Signs
differ; the mapping between the model's +z and the device's +y is not established, and
`ratcheting_pulse.py` already flags that same unresolved sign convention.)

So the evidence is genuinely split:

- **For the accelerometers:** at rest `atan2(a_y, a_z)` *is* the tilt, two independent
  chips agree to r = 0.995, and the displacement transducers show almost no permanent
  chord rotation.
- **For the tilt channel:** the model's independent prediction of top-of-wall rotation
  lands on it, not on the accelerometers. And there is a plausible failure mode on the
  other side — if the device periodically auto-zeroes its horizontal acceleration
  channels, the accelerometer-derived angle would under-read a slowly accumulated
  permanent rotation, and increasingly so as it grows. That would produce exactly the
  growing ratio observed.

I cannot adjudicate this from the data alone; it needs the instrument's documentation or
a bench check. What I can say is that **it is a real inconsistency inside one device, it
is reproducible across two independent specimens, and it changes the headline conclusion
either way**:

- If the accelerometers are right, permanent local rotation is 0.27° / 0.074° and the
  model over-predicts top-of-wall rotation by ~10×.
- If the tilt channel is right, the model reproduces top-of-wall rotation to within a
  factor of ~2 — and **the ratcheting discrepancy is largely an artefact of having
  compared the model's full-wall chord against a top-of-wall local rotation**, which is
  precisely the §1 category error.

`check_tiltmeter_consistency.py` reproduces the comparison from the raw logger files.

### A note on the figure itself

The dimension chain reads 800 / 660 / 600 top-down, putting the bottom quarter level at
**600 mm**. That is the **US-2 (Test 12)** geometry — `Test9_Info.xlsx` puts US-1's bottom
sensor at 660 mm. If this figure is being used to describe US-1, the bottom quarter level
is wrong by 60 mm. It does confirm the ±540 mm offsets of sensors 3 and 4 and the
sensor-to-channel mapping used throughout.

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
   `0.025·(2·f_t)^0.7 × 10³` = 13.2 J/m² for `f_t` = 0.2 MPa — which is *inside* van der
   Pluijm's 5–20 J/m² band for clay-brick bed joints, so G_I is not obviously the
   culprit. The stronger suspect is `cohesion-residual 0`: once a course opens it
   retains nothing at all, so damage cannot redistribute. A more ductile mode-I response, or the
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

0. **Resolve the tiltmeter channel disagreement** (§1). Everything about ratcheting
   depends on whether the experimental permanent rotation is 6° or 0.27°. Check the
   instrument's calibration/datasheet, and if possible re-read the raw device rather than
   the exported summary. Run `check_tiltmeter_consistency.py` to see the evidence.
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
