# Strategy: matching top-of-wall rotation and OOP response

Companion to `MODELLING_DISCREPANCIES.md`. That document diagnoses; this one proposes
what to change, in what order, and how you would know it worked.

---

## Correction: neither specimen collapsed

Moshfeghi, Smyrou, Arslan & Bal, *Structures* **66** (2024) 106815 settles several things
I had inferred wrongly from the data alone. Most importantly:

**US-1 and US-2 were never taken to collapse.** The test sequence stopped because the
shake table reached its capacity — PTA 0.78 g — not because the walls failed. The paper
is explicit that the tests "could not be continued with higher accelerations", and that
the cracks visible at the start "were also visible, with almost no change in length and
width, at the end of the test sequence". Run 25 was simply not applied to US-1.

So **"OOP displacement at collapse" is not an experimental observable in this campaign.**
There is nothing to match. My earlier framing — that US-1 "collapsed at ~29.7 mm" — was
wrong; that is just the largest displacement the table could impose (Table 4: US-1
+25.0 / −30.0 mm, US-2 +10.0 / −8.0 mm).

That does not make the model look better. It makes the target different:

| | US-1 | US-2 | Model (NODAMP v6) |
|---|---|---|---|
| Period elongation over the whole sequence | **+17.5%** (0.091→0.107 s) | **+13.9%** (0.087→0.099 s) | **+337%** (4.37×) |
| Max OOP displacement at top quarter | 30 mm | 10 mm | 46.7 mm by run 23, 194 mm at run 25 |
| End state | damaged, standing, cracks barely grown | damaged, standing | full-height rocking mechanism from run 14 |

The measured period elongations confirm the numbers I derived independently, and the
model's 4.37× is the central problem. **The model is producing a collapse mechanism where
the experiment produced only moderate, distributed damage.** The goal is therefore not to
match a collapse displacement but to *stop the model collapsing at all* within this input
range — while still reproducing the ~30 mm peak displacement and the residual tilt.

### What the paper says about the wall section — and why the collar-joint idea is dead

The specimens are **solid one-brick walls**, not two-leaf cavity walls. They are built
"using the one-brick technique, where the wall is as wide as the long edge of a brick"
(207 mm brick, 210 mm wall), alternating stretcher and header courses, "The head course
is laid with the short side of the brick exposed to greatly increase the structural
integrity". The paper states plainly that **"the two layers in the solid walls are fully
engaged due to the bricklaying technique"**.

So there is no weak collar joint, and my §2.2 hypothesis — that the model's collar joint
was too strong and was suppressing leaf separation — is **withdrawn**. The cavity-wall
specimens from the same campaign are reported separately precisely because they behave
differently.

This also answers the question directly: giving the mid-thickness plane `mason` joint
properties is defensible. If anything the model errs the *other* way — `Geo_Prep.dat` cuts
a **continuous** plane at z = 0.105, whereas in the real wall every header course crosses
that plane with solid brick. The model therefore under-represents the through-thickness
bond rather than over-representing it. Two options, in order of preference:

1. Leave the plane cut but give it **higher** tension/cohesion than a bed joint, to stand
   in for the header interlock; or
2. Do not cut it at all, and treat the wall as monolithic through its thickness, which is
   what "fully engaged" implies.

Neither is likely to fix the premature mechanism — the model's problem is at the bed
joints and the boundary conditions, not through the thickness.

### The ratcheting mechanism the paper identifies, which the model may be double-counting

Two statements matter a great deal for `ratcheting_pulse.py`:

- *"The hysteretic curves are asymmetric, an effect caused by the timber floors applying
  a moment towards one side of the wall section."*
- *"In the unstrengthened specimens, the joists slipped towards the positive direction
  (i.e. far from the wall surface) **8 to 12 mm cumulatively** at the end of all the
  tests… The bow ties, however, fully couple joists to the wall limiting the slipping of
  the embedded timber joists to almost zero."*

And decisively: the specimen with bow ties (ST-2-HB), which eliminated joist slip,
developed residual tilt **in the opposite direction** to every other wall.

So the experimental asymmetry and ratcheting are attributed to **eccentric joist loading
and cumulative joist slip** — a mechanism your model already contains geometrically.
`ASYM_K` imposes an asymmetric *input pulse* derived from the ground-motion record, which
is a different mechanism entirely. There is a real risk it is double-counting an
asymmetry the model should be generating on its own.

**Test this before anything else in §2**: run the model with `USE_RATCHETING = False`
(symmetric pulse) and see how much residual tilt the joists alone produce. If the model
already ratchets without an asymmetric pulse, `ASYM_K` should be dropped, not fitted.

You also now have a direct, quantitative validation target that needs no new
instrumentation: **cumulative joist slip, 8–12 mm**. Channels 8–11 already record joist
displacements. Compare.

---

## 1. Matching the top-of-wall rotation

### 1.1 Compare the right pairs (do this now, no re-run needed)

This is bookkeeping and it costs nothing:

| Experimental series | Correct model counterpart | Currently compared against |
|---|---|---|
| `exp_US{1,2}_tiltmeter.csv` settled values | `tilt_beam_seg` (2.06→2.68 m), or `tilt_local_incl` once available | `tilt_full_wall` — **wrong, and 4× smaller and opposite in sign** |
| `exp_Test{9,12}_tilt.csv` | `tilt_full_wall` | correct |
| inclinometer **peak** columns | *nothing* | plotted — must stop |

Just re-pairing these moves the model from 7× under the tiltmeter to within ~2× of it.
That is the single highest-value change in this document per unit of effort.

### 1.2 The instrument dispute is resolved — in favour of the tilt channel

Three independent lines now agree that the tiltmeter's tilt channel is the trustworthy
one and its accelerometer channels are not usable for static tilt:

- **The deformed shape.** Fig. 15 shows US-1 reaching ~28–30 mm at 2.06 m with the top
  beam held at zero at 2.85 m. The chord rotation over that upper segment is
  `atan(28/790)` ≈ **2.0°** at peak — the same order as the reported residual tilt, and an
  order of magnitude above the 0.27° the accelerometers give.
- **The authors' own reading.** They interpret the 2.2–6° residual tilt as physically
  consistent with the deformed shapes, and note that the bow-tie specimen reversed its
  sign — a coherent physical story, not an instrument artefact.
- **Your model.** `tilt_beam_seg` (2.06→2.68 m) gives −1.51° at run 21 and −2.37° at
  run 23, right in that range.

Table 3 also explains the anomaly: the tiltmeter reading is *"Separately Retrieved from an
Online Portal"*. The tilt value is the vendor's processed output, not a raw channel, so
there is no reason to expect the raw accelerometer words in the same file to be
DC-preserving — auto-zeroing or high-pass filtering in that pipeline would produce exactly
the under-reading observed.

So: **compare against the tilt channel, use `tilt_beam_seg`, and disregard the
accelerometer-derived angles.** The earlier "0.27°–6.04° bound" is withdrawn.

The peak columns remain unusable — they are acceleration-dominated and clipped at
±25.75°, and the paper confirms the instrument was only ever read pre- and post-run.

### 1.3 Resolve the sign convention

The model's `tilt_beam_seg` is negative where the tiltmeter is positive, and
`ratcheting_pulse.py` already carries an unresolved note about whether the record's `s`
maps to 3DEC's +z. Fix both together: pick one physical direction (say, "away from the
shake table's positive drive direction"), and assert it in one place. Right now a sign
error and a real counter-rotation are indistinguishable.

### 1.4 Set the sensor height correctly

Table 3 gives the tiltmeter at **X = 0, Z = 2.43 m**. Set `[incl_y_lo]`/`[incl_y_hi]` to
bracket 2.43 m — and note that `tilt_beam_seg` (2.06→2.68) already straddles it, which is
why it is a good proxy. The per-leaf channels I proposed earlier are dropped: the wall is
solid and fully engaged through its thickness, so there are no leaves to separate.

---

## 2. Matching the OOP response and the damage mechanism

### 2.1 You cannot match a quantity the algorithm never computes

`strategy_C_3dec_nodamp.py` has **no collapse criterion at all**. It runs
`model solve dynamic time <pulse_dur>`, then the gap, saves, and proceeds. There is no
convergence limit, no displacement threshold, no abort. `check_run_complete()` is defined
at line 192 and never called. Blocks can be flying apart and the run still "completes".

So the first change is to make collapse an *outcome the model can report*:

```python
# after the pulse and the settling gap, before saving
COLLAPSE_DISP_MM   = 210.0   # = wall thickness; instability reference
COLLAPSE_TILT_DEG  = 5.0     # full-wall chord
COLLAPSE_VEL_MPS   = 0.5     # still moving fast after the gap = not settling

def assess_stability(run_no):
    """Return (state, metrics). state in {'stable','marginal','collapsed'}."""
    rel = read_channel(run_no, "rel_disp_top_mm")
    tilt = read_channel(run_no, "tilt_full_wall")
    vmax_tail = peak_abs(read_channel(run_no, "Channel_17_AccTop")[-tail:])
    resid = tail_mean(rel)
    if abs(resid) > COLLAPSE_DISP_MM or abs(tail_mean(tilt)) > COLLAPSE_TILT_DEG:
        return "collapsed", ...
    if vmax_tail > COLLAPSE_VEL_MPS:
        return "marginal", ...          # not settling within the gap
    return "stable", ...
```

Then record the state per run in `strategy_C_summary.csv`, and **stop the sequence on
`collapsed`** rather than grinding through the remaining protocol. The run at which that
first triggers, and the peak displacement in it, is your model's "OOP displacement at
collapse" — the quantity you actually want to compare.

Normalise it: report displacement as a fraction of the section thickness (210 mm
composite, 105 mm per leaf). The specimen failed at 0.14 t. That framing makes the
leaf-separation hypothesis immediately testable — if separation is real, the specimen
failed at 0.28 of a *single-leaf* thickness, still low but far more plausible.

### 2.2 Suspects, revised

With the collar joint withdrawn, the candidates for "why does the model form a full-height
mechanism at run 14 when the specimen never did" are:

1. **Bed-joint fracture energy is too brittle.** `G_I = 0.025·(2·f_t)^0.7 × 10³` ≈
   **5.4 J/m²** at `f_t` = 0.2 MPa. The paper's measured bond-wrench strength is
   **0.28 MPa** (CoV 22%) and the wallette flexural strength 0.41–0.45 MPa, so the
   strength is about right but the *ductility* is the free parameter. Too brittle a
   mode-I response forces damage to localise into one mechanism instead of distributing.
   This is my primary suspect now.
2. **The top boundary condition.** The paper: the top beam "is restrained horizontally,
   but it is free to rotate", and the horizontal displacement at the top of the wall "is
   always zero". Your driver applies the table velocity to **both** group `'S'` and group
   `'T_B'`, which imposes the table motion on the top beam rather than holding it fixed
   in space. Those are different boundary conditions. Check this — it is cheap and it
   directly controls whether a full-height single-curvature mechanism can form.
3. **Hinge locations.** The paper reports cracks at the base *and at the floor levels*,
   with a hinge at the lower slab level in the unstrengthened specimens (Fig. 15, Fig. 17).
   Check where the model's joints actually open at run 14 — if they are not at the floor
   levels, the mechanism is wrong regardless of the displacement magnitude.
4. **Joist support model.** The joist contacts currently have essentially zero strength
   (`stiffness-shear 1, tension 1, cohesion 1`), i.e. frictionless sliding. The experiment
   shows 8–12 mm of *cumulative* slip, which implies real friction with slip-dependent
   accumulation, not free sliding.

### 2.3 Do not tune damping to fix the mechanism

Rayleigh 1.5% brings the peak-response ratio from 2.05× to 0.64× in the rocking regime,
so it will superficially improve the displacement match too. Resist it. The
two-regime behaviour in §4 of the discrepancies document shows the dissipation is
mis-assigned, and matching a collapse displacement by adding velocity-proportional
damping to a rocking mechanism buys a number at the cost of the physics.

---

## 3. Fixing the adaptive algorithm itself

### 3.1 Period identification must fail loudly

Currently `identify_Tend_from_csv()` returns `T1_init` when it cannot find a period, which
is indistinguishable from a genuinely unsoftened wall and silently destroyed run 24.
Minimum change:

```python
T_end, ok, source = identify_Tend_from_csv(...)   # ok is new
if not ok:
    T_end = T_prev_successful                     # not T1_init
    print("*** PERIOD ID FAILED on run %d - carrying forward %.4f s" % (run_no, T_end))
summary_row["period_id_ok"] = "Y" if ok else "N"
summary_row["period_id_source"] = source          # 'autocorr' | 'fft' | 'carried'
```

Carrying forward the last good period is far less wrong than reverting to the undamaged
one. Writing the flag into the summary lets every downstream figure exclude affected runs
automatically instead of relying on spotting `1.0000`.

Consider also aborting outright: two ID failures in the final three runs is not a
robustness problem to paper over, it is a signal that the ring-down has stopped looking
like a decaying single mode — which is itself physically meaningful (see §3.2).

### 3.2 Period-based scaling is ill-posed once a mechanism forms

This is the deeper methodological issue. A rocking block is **not** a linear oscillator:
its equivalent period depends on amplitude, increasing as the rocking amplitude grows.
Once the model enters the rocking regime at run 14, "the identified period" is a function
of how hard the previous run hit it, and tuning the next run's `Sd` to that period chases
a moving target. The tell is in your own data: `T_end/T₀` bounces around 3.9–4.4 with no
trend across runs 15–22 while the response keeps growing, and the amplification factor
reaches 22×.

Two defensible options:

- **Switch control variable at mechanism onset.** Detect onset (a step change in `T_end`,
  or the §2.1 stability metric going marginal), then stop scaling by `Sd(T)` and scale by
  **target displacement** or **input energy** instead. This is the standard move for
  collapse-oriented IDA on rocking systems.
- **Keep period control but report it honestly.** Freeze the excitation period at its
  pre-mechanism value and let the amplitude do the work, so at least the input is a
  well-defined function of the protocol rather than of an ill-conditioned identification.

I would take the first. It also removes the failure mode in §3.1 entirely, because after
onset you no longer depend on identifying a period at all.

### 3.3 Detect mechanism onset explicitly

Right now the run-14 transition is only visible in hindsight. Make it a first-class
output: flag the run where `T_end/T_prev > 1.5` or where the number of open joints jumps.
That single number — *the run at which a mechanism forms* — is directly comparable
between model and experiment (the specimen: never, up to collapse), and is a much more
discriminating validation target than peak displacement.

---

## 4. Fixing the dissipation, not the damping — MEASURED

`hysteresis_dissipation.py` builds base-shear vs top-quarter-displacement loops for model
and US-1 using the paper's own definitions (U = mean of Ch3/Ch4 per Eq. 1, F = base shear
per Eq. 2) and measures the enclosed area. Both signals are low-passed at 50 Hz first —
the model's contact force carries impact chatter that the experimental accelerometer-derived
base shear does not, and comparing raw traces is not like-for-like.

**Result 1 — yes, the model under-dissipates before the mechanism forms, by about 2×.**

| Loop fullness (area / 4·F_pk·U_pk, per cycle) | Model | US-1 | Ratio |
|---|---|---|---|
| Runs 5–12 (pre-mechanism) | 0.036 | 0.070 | **0.51×** |
| Runs 14–21 (post-mechanism) | 0.103 | 0.136 | 0.76× |

The pre-mechanism deficit is real and robust: it survived three different formulations of
the metric (0.24× whole-run-signed raw, 0.44× whole-run-signed filtered, 0.51× per-cycle
filtered) and is insensitive to whether `total_shear` or `cstav` is used. Both model bands
sit at or below `xi_eq` ≈ 0.02–0.05, i.e. the joints are barely dissipating.

**Result 2 — and this is the bigger finding: the model has almost no strength left after
run 14.**

| Run | Model peak base shear | US-1 peak base shear | Model E_diss | US-1 E_diss |
|---|---|---|---|---|
| 14 | 4.58 kN | 18.39 kN | 17.9 J | 34.3 J |
| 17 | 5.04 kN | 23.40 kN | 62.0 J | 137.4 J |
| 20 | 2.14 kN | 28.76 kN | 47.8 J | 273.5 J |
| 21 | 2.86 kN | 31.92 kN | 53.0 J | 350.0 J |

The model's base shear **decays to 2–3 kN while the experiment's climbs to 32 kN** — a
10× strength deficit by run 21 — and it dissipates 7× less energy in absolute terms.

2–3 kN is about what **pure gravity-stabilised rigid rocking** would give: with roughly
27 kN of total vertical load (12.9 kN self-weight + 10.35 kN spring + 3.84 kN slabs), a
mid-height hinge mechanism gives on the order of `N·t/h_half` ≈ 27 × 0.21 / 1.29 ≈ 4.4 kN.
So after run 14 the model is behaving as a rigid block held up only by gravity, with the
joints contributing nothing. The specimen, still developing 32 kN, plainly retained
substantial flexural and arching capacity throughout.

**What that points at.** The bed joints lose *all* cohesive capacity once cracked:
`cohesion-residual 0` is set explicitly, and `G_I = 0.025·(2·f_t)^0.7 × 10³` ≈ **5.4 J/m²**
at `f_t` = 0.2 MPa is very brittle. Once a joint opens there is no residual tension, no
cohesion, and friction only mobilises under normal stress — which an opening joint does not
have. Priority order:

1. **Give the joints residual capacity.** A non-zero `cohesion-residual`, and a larger
   `G_I`/`G_II`, so cracked joints still carry something. The paper's measured bond-wrench
   strength is 0.28 MPa (CoV 22%) and wallette flexural strength 0.41–0.45 MPa, so the
   *peak* strength is about right — it is the post-peak behaviour that is wrong.
2. **Check the arching path.** A wall restrained top and bottom develops arching thrust,
   which is the main reason the real specimen kept gaining strength. If the top boundary
   condition is wrong (§2.2 item 2) the thrust never develops, which would produce exactly
   this collapse in capacity.
3. **Only then** revisit viscous damping — and if you keep it, stiffness-proportional only,
   so it does not penalise the low-frequency rocking.

Fix the `DAMP_RATIO` trap in `strategy_C_3dec_nodamp.py` (lines ~440–446) before any
damping study, or it will silently produce undamped results while reporting otherwise.

**A caveat on the metric.** Zero-crossing cycle segmentation degrades for strongly
one-sided (ratcheted) responses — US-1 run 18 swings −9.0 to +2.1 mm and its cumulative
energy is under-counted as a result. The fullness values are sound; treat absolute `E_diss`
for heavily one-sided runs as a lower bound.

```bash
python hysteresis_dissipation.py --runs 5 12 --test 9 --lowpass 50
python hysteresis_dissipation.py --runs 14 21 --test 9 --lowpass 50
```

---

## 5. Validation plan

You now have **both specimens in the repo**, which makes a genuine out-of-sample test
possible for the first time. Use it:

- **Calibrate on US-1** — collar-joint strength, `G_I`, `ASYM_K`.
- **Predict US-2 blind**, changing only the two things that genuinely differ: the bottom
  sensor elevation (0.60 vs 0.66 m) and the 25th run. Do not re-tune.
- US-2 collapsed later and more gently than US-1 (final tiltmeter 2.37° vs 6.04°;
  accelerometer 0.074° vs 0.275°). A model calibrated on US-1 that also reproduces that
  *ordering* is meaningfully validated. One that only reproduces US-1 is fitted.

Acceptance targets, in the order they matter:

| # | Target | Current status |
|---|---|---|
| 1 | No full-height mechanism forms within the tested input range | model forms one at run 14; neither specimen did — **fails badly** |
| 2 | Period elongation ≈ +17.5% (US-1) / +13.9% (US-2) | model +337% — **fails badly** |
| 3 | Max OOP displacement within ±30% of 30 mm (US-1) / 10 mm (US-2) | model reaches 46.7 mm by run 23 — fails |
| 4 | Top-of-wall rotation vs the tilt channel (2.77° at run 21, 6.04° at run 24) | `tilt_beam_seg` gives 1.51° / 2.37° — within ~2×, passes once re-paired |
| 5 | Peak response ratio 0.8–1.25× across all runs | 2.10×/2.05×/6.40× undamped — fails |
| 6 | Cumulative joist slip 8–12 mm | never checked; channels 8–11 already record it |

Targets 1 and 2 are the real ones. If the model stops forming a premature full-height
mechanism, 3 and 5 will move a long way on their own — and `ASYM_K`, which is currently
absorbing the error from all of the above, can finally be fitted to something meaningful.

---

## Suggested sequence

1. **Check whether the model ratchets with `USE_RATCHETING = False`.** The paper
   attributes the asymmetry to eccentric joist loading and joist slip, which your model
   already has. If it ratchets without an imposed asymmetric pulse, `ASYM_K` is
   double-counting and should be dropped rather than fitted.
2. **Re-pair the tilt comparisons** (§1.1) and set the sensor height to 2.43 m (§1.4).
   Free, immediate, removes a category error.
3. **Check the top boundary condition** (§2.2 item 2). The paper holds the top beam
   horizontally fixed but free to rotate; the driver drives it with the table motion.
   Cheap to test and it directly controls whether a full-height mechanism can form.
4. **Compare cumulative joist slip against the measured 8–12 mm** using channels 8–11.
   A quantitative target that needs no new instrumentation.
5. **Add the mechanism-onset flag and the `period_id_ok` flag** (§2.1, §3.1), then
   attack the premature mechanism via bed-joint fracture energy (§2.2 item 1).
6. **Revisit dissipation** (§4), and validate blind on US-2 (§5).
