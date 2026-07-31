# Strategy: matching top-of-wall rotation and OOP displacement at collapse

Companion to `MODELLING_DISCREPANCIES.md`. That document diagnoses; this one proposes
what to change, in what order, and how you would know it worked.

---

## The one observation that should drive everything

| | US-1 (Test 9) | Model (NODAMP v6) |
|---|---|---|
| Max period elongation before failure | **1.22×** | **4.37×** (reached at run 18) |
| Period elongation at the run before failure | 1.20× (run 23) | 3.92× |
| Peak OOP displacement at failure | **~29.7 mm** (run 24) | survived **46.7 mm** at run 23 without failing |
| Character | fails **abruptly, while still stiff** | softens **massively, then survives 11 more runs** |

The specimen collapsed at roughly **14% of the wall thickness** (29.7 mm on a 210 mm
section) having barely softened. That is far below any rocking-instability limit — a
two-leaf wall failing by rocking would need mid-height displacement on the order of the
thickness itself, and would announce itself with large period elongation first.

**The model is not failing by the mechanism the specimen failed by.** It builds a
ductile, full-height rocking mechanism that tolerates very large displacement; the
specimen suffered something abrupt and comparatively brittle. No amount of tuning
`ASYM_K`, damping, or the pulse shape will reconcile those — they are different failure
modes.

So the strategy has to start with a question that is not a modelling question:

> **What actually failed in US-1 at run 24?**

Recover this from the test report, video, or post-test photographs before committing
compute. The candidates below all predict ~30 mm collapse displacement with little prior
softening, and they need different model changes:

| Candidate mechanism | Why it fits | What it would need in the model |
|---|---|---|
| **Collar-joint delamination, outer leaf peels** | Explains low capacity (a 105 mm leaf is half as stable) *and* why the tiltmeter on the wall face reads a large rotation while the chord stays plumb | Weaken the collar joint — see §2 |
| **Joist / slab connection loss** | Removes intermediate restraint abruptly; wall had been relying on it, so no prior softening | Joist contacts currently have near-zero strength *and cannot detach* — needs a real pull-out path |
| **Local crushing at a support or joist pocket** | Brittle, localised, little global softening | `comp-residual` and `G_c` govern; check whether crushing is even being triggered |
| **In-plane / OOP interaction** | The tiltmeter records both axes; in-plane residual was −0.36° | Model is driven OOP only |

Everything below assumes you will identify this first. §1 and §3 are worth doing
regardless, because they are correctness fixes.

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

### 1.2 Bound the answer while the instrument dispute is open

The tilt channel and its own accelerometers disagree by 20–30× and I could not
adjudicate. Until a bench check settles it, **report the experimental rotation as a
range**, not a point: 0.27°–6.04° for US-1, 0.074°–2.37° for US-2. A model result inside
that band is not yet falsified; one outside it is. This is honest and it keeps the
comparison usable in the meantime.

### 1.3 Resolve the sign convention

The model's `tilt_beam_seg` is negative where the tiltmeter is positive, and
`ratcheting_pulse.py` already carries an unresolved note about whether the record's `s`
maps to 3DEC's +z. Fix both together: pick one physical direction (say, "away from the
shake table's positive drive direction"), and assert it in one place. Right now a sign
error and a real counter-rotation are indistinguishable.

### 1.4 Add per-leaf channels

If the leaf-separation hypothesis is live, the model must be able to see it. Add
displacement and rotation histories on **each leaf separately** at the top of the wall —
the collar joint is at z = 0.105, so gridpoints near z = 0.05 (inner leaf) and z = 0.19
(outer leaf) at the tiltmeter height. Their difference is the delamination signal, and
the outer-leaf rotation is what a face-mounted tiltmeter would actually see.

This is the test that would unify the two anomalies: if the outer leaf rotates several
degrees while the wall chord returns to plumb, both the tiltmeter reading *and* the low
collapse displacement are explained by one mechanism.

---

## 2. Matching OOP displacement at collapse

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

### 2.2 The collar joint is the prime suspect

`ANALYSIS_PART_I_MASON.dat` line 46 assigns `mason_v7` to **every** contact, and the
collar joint cut at z = 0.105 in `Geo_Prep.dat` is a contact like any other. It therefore
carries `tension = 0.2 MPa`, `cohesion = 0.3 MPa` — the same as a bed joint.

Collar joints in two-leaf masonry are the classic weak plane; they are often barely
filled. Giving one bed-joint strength makes the two leaves act compositely, which turns
the wall into a monolithic 210 mm rocker with a large, ductile displacement capacity.
That is precisely the behaviour the model shows and the specimen did not.

**Proposed parametric study** — cheap, well-posed, and it tests a real hypothesis:

| Case | Collar-joint tension | Collar-joint cohesion | Rationale |
|---|---|---|---|
| C0 | 0.2 MPa (as now) | 0.3 MPa | baseline, composite action |
| C1 | 0.05 MPa | 0.075 MPa | weak but bonded |
| C2 | 0.01 MPa | 0.02 MPa | nominally unbonded |
| C3 | 0 | 0 | frictional contact only |

Implement by assigning a separate contact group at the collar plane before the global
property assignment, e.g. group the contacts in a thin `pos-z` band around 0.105 and give
them their own `block contact property ... range group 'collar'`. Run each to collapse
using the §2.1 criterion and plot collapse displacement against collar strength. If the
curve passes through ~30 mm at a physically plausible collar strength, you have your
mechanism — and a defensible calibration.

Watch the **period elongation** in these runs as much as the displacement. The target is a
case that collapses at ~30 mm *while staying near T/T₀ ≈ 1.2*. A case that collapses at
30 mm only after softening 4× has matched the number for the wrong reason.

### 2.3 Do not tune damping to fix collapse

Rayleigh 1.5% brings the peak-response ratio from 2.05× to 0.64× in the rocking regime,
so it will superficially improve the collapse displacement too. Resist it. The
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

## 4. Fixing the dissipation, not the damping

From the discrepancies document: undamped over-predicts peak response 2.10× before
softening, while 1.5% Rayleigh under-predicts 0.64× during rocking. One global ratio
cannot serve both, because the deficiency is in the joint's hysteretic dissipation before
the mechanism forms, and Rayleigh damping then over-penalises the slow large-amplitude
rocking afterwards.

Priority order:

1. Check whether the pre-mechanism response is dissipating anything at all — plot a
   `Channel_3` force–displacement loop for runs 5–12 and measure the enclosed area. If
   the loops are nearly closed, `G_I`/`G_II` are too small and the model is behaving
   near-elastically until it suddenly is not.
2. `G_I = 0.025·(2·f_t)^0.7 × 10³` ≈ **5.4 J/m²** at `f_t` = 0.2 MPa. That is at the
   brittle end. A more ductile mode-I response would let cracking distribute instead of
   localising into one full-height mechanism at run 14 — which is the §3 problem in the
   discrepancies document, and plausibly the same root cause as the premature softening.
3. Only after 1–2, revisit viscous damping — and if you keep it, use stiffness-proportional
   only, so it does not penalise the low-frequency rocking.

And fix the `DAMP_RATIO` trap in `strategy_C_3dec_nodamp.py` (lines ~440–446) before any
damping study, or it will silently produce undamped results while reporting otherwise.

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
| 1 | Run at which a mechanism forms: model ≈ experiment (i.e. **not before collapse**) | model run 14, experiment never — **fails badly** |
| 2 | Peak period elongation ≤ ~1.3× | model 4.37× — **fails badly** |
| 3 | OOP displacement at collapse within ±30% | model survives 1.6× the collapse displacement — fails |
| 4 | Top-of-wall rotation within the instrument bound | plausibly already passes once re-paired (§1.1) |
| 5 | Peak response ratio 0.8–1.25× across all runs | 2.10×/2.05×/6.40× undamped — fails |

Targets 1 and 2 are the real ones. If the model stops forming a premature full-height
mechanism, 3 and 5 will move a long way on their own — and `ASYM_K`, which is currently
absorbing the error from all of the above, can finally be fitted to something meaningful.

---

## Suggested sequence

1. **Identify the experimental failure mode** from the test records. Nothing else is
   well-posed until this is known.
2. **Re-pair the tilt comparisons** (§1.1) — free, immediate, and removes a category error.
3. **Add the collapse criterion and the `period_id_ok` flag** (§2.1, §3.1) — small code
   changes, and they make every subsequent run interpretable.
4. **Run the collar-joint parametric study** (§2.2) against collapse displacement *and*
   period elongation.
5. **Switch the control variable at mechanism onset** (§3.2) if the sequence still depends
   on post-mechanism period identification.
6. **Revisit dissipation** (§4), then re-fit `ASYM_K` last, and validate blind on US-2 (§5).
