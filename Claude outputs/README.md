# Route 2 — adaptive sinusoidal spectral matching on the SAFEGo US-1 wall

**Handover note, 2 October 2026.** Everything done on Route 2 to date, in the order it happened, with the figures, the numbers and the files. Written so that someone who has not followed the discussion can pick it up and run the next step.

---

## Background — why this work exists

**The experiment.** SAFEGo is a shake-table campaign on unreinforced-masonry (URM) walls representative of Dutch housing exposed to the induced seismicity of the Groningen gas field (Moshfeghi et al. 2024, *Structures* 66:106815). The walls are 1.292 × 2.58 × 0.21 m, built with timber floor joists bearing on them at two levels and a loading beam applying ~0.1 MPa of vertical stress, tested out-of-plane in one-way bending. Two unstrengthened specimens, US-1 (Test 9) and US-2 (Test 12), went through the same 25-run incremental protocol: a Groningen-type record (HU12) scaled from ×0.5 to ×6.0, with the 1940 El Centro record (EC40) interleaved at four runs, then the 1976 Friuli record (FR76) at ×1.0, 1.5, 1.75 and 2.0. Both reached rocking and large displacements in the FR76 runs; US-1 stopped after run 24 (29.5 mm at 2.06 m), US-2 after run 25 (31 mm). The campaign also tested strengthened walls, which are outside this note. Instrumentation: potentiometers at the bottom-quarter, mid and top-quarter levels (0.60 / 1.26 / 2.06 m, the top ones left and right at ±0.54 m), accelerometers at the same levels, a tiltmeter on the top courses, and load cells on the overburden rods.

**The model.** A rigid-block distinct-element model of the wall in 3DEC 9: every brick a block, joints with normal and shear stiffness, tensile strength, cohesion and friction, the joists and the loading beam as blocks with their own contacts. It was first calibrated against US-1 by running the actual 25 records in sequence ("Strategy F", the control: the damage state carried from run to run, the model's own response compared run by run with the test). That is the reference everything else is measured against, and it costs ~824 s of simulated time — hours to days of wall-clock on a DEM model — for one wall and one damping/geometry setting.

**The research question.** The cost of running real records through a DEM model is what stands between this kind of model and routine use: parametric studies, fragility, second specimens, damping variants all multiply it. The question is whether a *cheap surrogate input* — a short synthetic signal — can put the model through the same damage sequence and reach the same peaks as the records do, at a small fraction of the cost. Two ways of building such a signal were set up as competing "routes":

- **Route 1 (record-based, Eleni):** the surrogate takes its properties from the earthquake only. In its first form, one sine pulse per run with amplitude from the record's PGV and period from the record's dominant frequency. It delivered 10–27% of the record's spectral demand at the wall's effective period, because these records are not pulse-like and a single pulse at the record's period carries little energy where the wall responds. Its strongest form turned out to be *extraction*: cutting the strongest 1.6–4 s window out of the record itself (98% of the spectral demand at a tenth of the cost; the model's response was half the full record's because the wall's peak came after the cut — a longer window is the fix). That thread has its own driver (`route1_extract_3dec.py`) and is not covered further here.
- **Route 2 (structure-informed, Ihsan): adaptive sinusoidal spectral matching.** The surrogate takes its periods from the wall and its amplitudes from the record's response spectrum: a cycle at the wall's first period T1 and a cycle at its effective (cracked) period T_eff, each scaled so that a 5%-damped SDOF at that period sees the record's spectral displacement for that run. "Adaptive" because the signal is meant to follow the wall's state as it cracks; "spectral matching" because the amplitude rule is a two-point match to the record's Sd spectrum. The claim to be tested: that the damage sequence and the peaks of the 25-run protocol can be reproduced with ~1/10 of the record's simulated time.

**Where Route 2 stood when this note starts.** The two-wave version had been run over the full protocol. On the displacement plot it looked encouraging — run 24 within 10% of US-1 — and that is where the discussion with Ihsan and Eleni began: the early runs were far too high, the loops at run 22 looked soft, and the run-21–23 over-prediction meant a truncated protocol would have read as a failure. The sections below follow that discussion to its resolution: what was actually wrong (the target, not the model and not the wave's period), what the corrected target looks like, and what remains to be run.

**People.** Ihsan proposed Route 2 and the "periods from the record" refinement discussed in §5; Eleni holds the record-based route; the model, the drivers and the analysis in this note are Yopi's. The results feed a progress deck for the Hanze meeting and, in parallel, a second specimen (the EUCENTRE single-leaf wall of Graziotti et al. 2016) on which the method is to be generalised — that work is separate and not in this note.

---

## 0. One-paragraph summary

Route 2 replaces each shake-table record of the 25-run SAFEGo protocol with a two-cycle sinusoid whose amplitudes are set so that a linear SDOF at the wall's periods sees the same spectral displacement as it would under the record. It was run over the full protocol on the 3DEC model (small strain, no damping). It reaches the US-1 peak at run 24 (32 vs 29.5 mm) but over-predicts runs 1–13 by 3–5× and runs 21–23 by 2–5×, and — the decisive comparison — it over-drives the *same model* by 2–5× relative to the actual record. A like-for-like base-shear comparison showed the model's stiffness and capacity are the test's; the over-driving comes from the input. The cause is the target: Route 2 asks every run for the demand at the cracked-state period T_eff, which the wall only has near its 21 mm bond-failure displacement. A capacity-spectrum fixed point (period read from the pushover at the displacement the run will actually reach) reproduces the record-on-model and US-2 peaks within the scatter between the two test walls, and identifies the ~21 mm tipping point beyond which the wall runs away. The refined Route 2 keeps the two-wave sinusoid and the Sd amplitude rule, replaces the frozen T_eff with the per-run fixed point, keeps the original wave on runs that can tip, and solves the two amplitudes jointly. It has not been run yet. The next step is a 3-signal test on FR76 runs 22–25 from the saved run-21 state.

---

## 1. What Route 2 is (as run)

**Model.** 3DEC rigid-block model of the SAFEGo US-1 wall (1.292 × 2.58 × 0.21 m, 0.1 MPa precompression as an equivalent top-beam density). Out-of-plane axis = model z; height = y; wall centreline x = 0.646. Base save `Part_I_MASON_v8_SmallStrain.sav`, small strain, no viscous damping (contact dissipation only). Drive groups `S` and `T_B` by `block apply velocity-z table`. Displacement channel: `rel_disp_top_exp_mm` = ½(Ch3+Ch4) − Ch5 at 2.06 m, the experiment's definition.

**Protocol.** 25 runs, damage carried between runs: HU12 ×0.5 → ×6.0 (runs 1–21, with EC40 at 3, 6, 8, 11), then FR76 ×1.0, 1.5, 1.75, 2.0 (runs 22–25). US-1 (Test 9) ran 1–24; US-2 (Test 12) ran 1–25.

**The signal, every run.** One velocity cycle at T1, then one at T_eff, no gap, same phase, then a 2.5 s ring-down tail.

- T1 = 0.112 s and T_eff = 0.241 s are fixed properties of the wall from the bilinear idealisation of the model's pushover (`bilinear_periods.json`, block `model`): Vu 29.2 kN, K1 5.16 kN/mm, dy 5.7 mm, du 21.1 mm, μ 3.72, m_eff 1635 kg. T_eff is the secant period at 0.8·Vu on the descending branch, i.e. at du.
- Each cycle's amplitude is calibrated alone so that its own 5%-damped Sd at its own period equals the run scale × the record's Sd at that period (`spectrum_HU12/EC40/FR76.csv`).
- Cost: 0.35 s of table per run + tail → ~71 s of model time for 25 runs, against ~824 s for the records (21 × 35.9 s HU12 + 4 × 17.3 s FR76).

Driver: `scripts/route2_full25_revised_3dec.py` (the version in this folder has the partial-re-run switch of §4 added; set `RESUME_AFTER = 0` for the original behaviour). Log: `data/route2_full_log.csv`. Postprocess output: `data/route2_full_metrics.csv` (US-1 = `exp_*`, US-2 = `exp2_*`, record on the model = `ctrl_*`, Route 2 = `r2_*`).

---

## 2. Result of the 25-run sequence

![Route 2 full sequence](img/01_route2_full25_vs_record_vs_tests.png)

*Left: peak OOP displacement at 2.06 m per run. Middle: peak base shear (model = `cstav`, support reaction only — see §3). Right: full-wall chord tilt.*

Peak displacement at 2.06 m, mm (exact values from `route2_full_metrics.csv`):

| run | input | US-1 | US-2 | record on model | Route 2 |
|---|---|---|---|---|---|
| 1 | HU12 ×0.5 | 0.30 | 0.11 | — | 0.79 |
| 2 | HU12 ×0.75 | 0.41 | 0.39 | — | 2.47 |
| 7 | HU12 ×1.5 | 1.24 | 0.85 | 2.83 | 2.25 |
| 13 | HU12 ×2.5 | 2.57 | 1.93 | 2.29 | 4.05 |
| 18 | HU12 ×4.5 | 8.48 | 5.19 | 3.97 | 6.32 |
| 21 | HU12 ×6.0 | 16.37 | 10.11 | 16.97 | 21.54 |
| 22 | FR76 ×1.0 | 3.65 | 2.38 | 3.72 | 8.11 |
| 23 | FR76 ×1.5 | 13.15 | 5.26 | 5.89 | 28.57 |
| 24 | FR76 ×1.75 | 29.48 | 7.58 | 8.29 | 32.09 |
| 25 | FR76 ×2.0 | — | 31.14 | 41.63 | 38.62 |

What it says:

1. **Route 2 matches US-1 at run 24 by coincidence.** The record on the *same model* gives 8.3 mm at run 24, next to US-2's 7.6. Route 2 is not reproducing the model's response to the record; it is over-driving the model by 2.2× (run 22), 4.9× (run 23) and 3.9× (run 24).
2. **The early runs are over-driven too**: 2.6–11× the test at runs 1–6, falling to ~1.2× by runs 16–20.
3. **US-1 and US-2 differ by 4× at run 24** (29.5 vs 7.6) on nominally identical walls. "Validation against US-1" is a weaker statement than it looks; the record-on-model line is the comparison that isolates the input.
4. **Ihsan and Eleni's reading** (truncate before run 24 and Route 2 looks 2–5× too high; the loops look soft at run 22) is correct as a description. §3 shows it is the input, not the model.

---

## 3. The loops, and the like-for-like base shear

![Loops with cstav](img/02_loops_runs21-25_cstav.png)

*Base shear vs displacement, runs 21–25, model force = `cstav`.* The model loops (red) look flat and "soft" next to the steep test loops at run 22.

Two things made that comparison invalid:

- **The model force is the wrong quantity.** `cstav` is the base + joist support reaction. The test's base shear (Moshfeghi et al. 2024, Eq. 2) is F = Σ mⱼ aⱼ: absolute acceleration at three wall accelerometers (0.60, 1.26, 2.06 m) × tributary wall mass (408.5 / 320.7 / 404.1 kg), × a processing factor **2.317348** whose origin is unknown. The model's own pushover puts ~27 kN at 8 mm; the cstav loop shows 8 kN there — two-thirds of the inertia force is not in the channel.
- **The experimental force carries the ×2.317.** 2.317 × wall mass × g over the cross-section is 0.095 MPa, so it is probably "the 0.1 MPa of vertical load turned into mass" — physically doubtful (10.4 kN of it is spring force), but it scales both sides equally if applied to both. Open question for Moshfeghi: what is 2.317, and did the published loops use accelerometers 15/16/17 only? Also: the paper gives wall density 1884 kg/m³, the masses use 1619.

**Fix.** `scripts/instrument_baseshear_exp.dat` records the absolute OOP velocities at the three accelerometer points; `scripts/baseshear_exp_from_model.py` differentiates them (50 Hz low-pass), applies the same masses and the same 2.317, and plots model and test on the same axes. Runs 21–25 were re-run from `route2_run_20.sav` with the channel added (driver switch `RESUME_AFTER = 20`, `OUT_DIR = route2_full25_revised_bs`).

![Like-for-like loops](img/03_loops_runs21-24_like_for_like.png)

*Same formula on both sides. Grey = the old cstav curve for reference.*

Reading:

- **Stiffness:** at run 22 the model loop climbs at ~4–5 kN/mm through the origin against 7–8 for the test — a factor ~1.5, within what the mass bookkeeping can move. Not a factor 8. Joint normal stiffness kn is not the problem; the run-1 free vibration of the model gives 10.15 Hz, inside the 9.6–12 Hz band the test transfer functions show.
- **Capacity:** the model's envelope flattens at 25–35 kN, where its pushover plateau (29 kN) and the US-1 loops (27–30 kN) sit.
- **What differs is excursion, not slope:** model 8 / 20 / 30 mm, test 3 / 10 / 25 at runs 22–24. The "soft" loop at run 22 is a wall that has uplifted and is rocking (flat post-uplift branch); the test loop is a wall that has not (4 mm < dy). Same wall, different regime, because the pulse pushed the model past uplift and the record did not.
- Two genuine model questions the loops raise, unrelated to kn: the model loops are fatter (more dissipation per cycle — joint sliding? check the joist-slip channels) and less asymmetric than the test's (the slab-joist inertia the paper describes).
- The −65 kN spike at run 24 is a contact impact in an undamped model; discount it.

---

## 4. Diagnosis: the target, the leak, the pulse

Three defects in Route 2 as run, in order of importance.

**(a) The target.** Wave B asks for the record's Sd at T_eff = 0.241 s on every run. T_eff is the secant period the wall has near du = 21 mm. At runs 22–24 the record on the model and the test put the wall at 3.7 / 5.9 / 8.3 mm — on the steep branch, with a secant period around 0.11–0.14 s, where the record's Sd is 2–6 mm, not 13–24. Route 2 hands the wall the demand of a state it has not entered, and the wall obliges by entering it. This is the whole of the 2–5× over-driving.

**(b) The T1 leak.** The two waves are calibrated separately, but wave B also excites a T1 oscillator. The driver logs it: `SdA_total_mm` vs `SdA_mm` in `route2_full_log.csv` — run 22: 3.48 vs 1.75 (2.0×); run 1: 0.65 vs 0.48. The elastic runs get up to twice their T1 target. Fix: solve the two amplitudes together on the combined signal.

**(c) Pulse severity.** To deliver Sd(T_eff) in one cycle the FR76 wave needs 0.44 g where the record peaks at 0.35 g (×1.75: 0.78 vs 0.62 g), in a single coherent swing. A cracked wall rocks on pulses. Two cycles on wave B at the same Sd halve the PGA (0.22 g).

---

## 5. Ihsan's proposal and what it turned into

Ihsan's comment (three versions, same content): keep the two waves and the Sd amplitude rule ("the amplitude pushes the structure to the displacement limit representative of the SDOF demand"), but take the wave periods from the record's predominant period instead of the wall's, so the input carries the record's frequency content on both HU12 and FR76; "FR76 ×1.00 triggers something in the structure that the record does not."

What was found:

- **"Predominant period" is not one number.** For FR76: Sa-peak 0.26 s (Rathje et al. 1998's T_p), mean period T_m 0.37 s (Rathje's recommended parameter), Fourier-amplitude peak 0.50 s (Kramer 1996's definition), velocity-pulse fit 1.28 s (Mavroeidis & Papageorgiou 2003 Eq. 1; R² 0.43 — FR76 is not pulse-like, and Baker's wavelet PI would say the same). For HU12 the Sa-peak is 0.12 s, T_m 0.206 s.
- **A record gives one period; two waves at one period cannot hit two targets.** Literally applied, the proposal collapses wave A on FR76 (the 0.26 s wave already covers T1) and is infeasible on HU12 (a 0.12 s wave B cannot deliver the T_eff target without over-driving T1 by 50%).
- **Solving an amplitude from Sd at a period the wave does not contain is unstable.** A 1.28 s pulse tuned to Sd(T_eff) needs ×9.7 — 0.86 g. Carriers more than ~2× off the structure's period are not usable.
- **Target vs carrier.** His comment changes the *carrier* (the wave's period). §4(a) says the *target* is what is wrong. The two are separable: the target can be corrected (§6) and the carrier can still be the record's period (§8), on FR76 at a *lower* pulse than the wall-period wave.
- **The circularity objection** ("we only take input based on what the structure outputs") is answered by distinguishing structure *properties* (the pushover, measured once, independent of the input — same status as T1 and T_eff, which he accepts) from structure *response to the input being built* (ring-down periods, secants from the last run — Strategy C-SECANT, which did need priming, clamps and growth caps). The capacity-spectrum method (N2, Eurocode 8 Annex B) reads the record's spectrum at the period set by the capacity curve and is not considered circular. Rule: the structure may choose where the record's spectrum is read; it may never modify the spectrum.

---

## 6. Step 1 — capacity-adaptive targets (computed, not run)

Script: `scripts/casm_step1_targets.py`. Output: `data/casm_step1_targets.csv`.

**Method.** For each run, solve the fixed point d\* = scale × Sd_record(T_sec(d\*)) with T_sec(d) = 2π√(m_eff / (F(d)/d)) from the pushover backbone (`pushover_pos.csv`, |F_tot| vs d_ctrl, m_eff 1635 kg). This is the N2 / capacity-spectrum intersection, per run. All crossings of the demand curve with the 1:1 line are reported: the first (from above) is the stable amplitude; a crossing from below is a tipping point; demand still above d at 70 mm is flagged "runaway".

**Damage memory (one input to state).** Once the method's own d\* has passed dy (run 19), the backbone's initial branch is softened by (T1/T_end)² with T_end/T1 = **1.30** — the ring-down elongation measured on the model after cracking in Strategy C. Plateau and descent are geometry + overburden and stay as they are. Without this, runs 22–24 come out at 1.7 / 2.5 / 3.2 mm; with it, 3.4 / 5.1 / 6.0. It is a property measured once, not a per-run feedback, but it is the one number here that came from a model response.

![Step 1 targets](img/04_capacity_adaptive_targets_step1.png)

*Left: what each run is asked to reach — Route 2's targets (Sd at T_eff and at T1), the fixed-point target d\*, and the measured peaks of US-1, US-2, the record on the model and Route 2. Right: the secant period of the wall as a function of amplitude, virgin and cracked backbones; the dots are the 25 runs.*

**Score, runs 7–24, geometric mean of target ÷ measured (×/÷ scatter):**

| | vs record on model | vs US-1 | vs US-2 |
|---|---|---|---|
| capacity-adaptive d\* | 0.71 ×/÷ 1.49 | 0.69 ×/÷ 1.45 | **1.08 ×/÷ 1.15** |
| Route 2 as run | 1.50 ×/÷ 1.62 | 1.47 ×/÷ 1.62 | 2.29 ×/÷ 1.61 |

Runs 22–24: d\* 3.4 / 5.1 / 6.0 vs record-on-model 3.7 / 5.9 / 8.3 vs US-2 2.4 / 5.3 / 7.6. Runs 1–6: d\* 0.4–1.0 vs US-1 0.14–0.96 (Route 2: 0.8–2.8). Runs 1–18 come out elastic, period 0.103 s (the backbone's initial stiffness gives 0.103, closer to the dynamic 0.0985 than the bilinear's 0.112).

**Where it is low, and why that is honest.** Runs 7–11: the record on the model gives 2.3–2.8 mm where the test gives 0.9–1.4; d\* follows the test. Runs 19–21: d\* follows US-2 (0.9–1.0×), not US-1 and the record-on-model (0.55–0.66×) — the HU12 ×5–6 runs are where the two walls diverge.

**The tipping point.** For every FR76 run the demand curve crosses the 1:1 line twice; the upper crossing is at **~21 mm = du**, the bond-failure drop in the pushover. For ×1.75 and ×2.0 the demand stays above capacity beyond it: a wall pushed past ~21 mm cannot find an equilibrium amplitude and runs away. That is US-1 at run 24 (29.5), the record on the model at run 25 (41.6), US-2 at run 25 (31.1) — and not US-2 or the model at run 24. The margin, scale × Sd(T_sec(du)) / du, is 0.65 / 0.97 / 1.14 / 1.30 at runs 22–25. Note that Sd(T_sec(du)) is Sd(T_eff): **Route 2's number was the right number for the wrong job** — a tipping criterion, not an amplitude target.

**Limit.** Equivalent linearisation predicts the stable amplitude. Run 25 (d\* 6.9 vs 41.6 / 31.1) is the runaway branch, which no linearised target follows; the method flags it (margin 1.30) but cannot predict it. At margin ≈ 1.1 (run 24) the real outcome was bistable — US-1 tipped, US-2 did not — so no deterministic surrogate reproduces both walls there.

---

## 7. A record-only reference: the spectrum-compatible synthetic (parked)

Script: `scripts/fig_synthetic_vs_route2.py`; table: `data/synthetic_FR76_x1p00_vel.txt`.

![Synthetic vs Route 2](img/05_synthetic_vs_route2_FR76.png)

14 sinusoids, log-spaced 0.05–1.0 s, random phases (seed 1), envelope from the record's Arias build-up over the 12–16 s window, amplitudes iterated until the sum's Sd matches the record's at the 14 periods. No wall property in it. It reproduces the record's spectrum across the band (all nodes within 4%), PGA 0.35 g (record 0.35), 4 s long (a quarter of the record). It does **not** reproduce the record's dominant −0.25 m/s pulse (synthetic −0.17): a pulse is a phase property and the phases are random. It is standard spectrum-compatible generation, not the research contribution; keep it as an optional reference line that shows what full spectral content costs, or drop it. The extraction of the same 12–16 s window (Route 1 extraction driver, `route1_extract_3dec.py`) has the pulse and the same cost.

---

## 8. Route 2 refined — the specification (not yet implemented)

Keep everything that makes it Route 2: two sinusoidal cycles, amplitudes from the record's Sd, structure periods from the pushover. Change `build_signal` only:

| | Route 2 as run | Route 2 refined |
|---|---|---|
| wave A | 1 cycle at T1, target Sd(T1) | same |
| wave B period | T_eff, every run | T\* = T_sec(d\*) from the per-run fixed point (= T1 while elastic, lengthens as the wall cracks) — or the record's period, see below |
| wave B target | Sd(T_eff), every run | Sd(T\*) = d\* |
| runs that can tip | — | if scale × Sd(T_eff) ≥ du: T_eff and Sd(T_eff) — Route 2 as run, because that is what it was good at |
| amplitudes | each wave alone | solved together on the combined signal |
| cracked-wall memory | none | initial branch softened by (T1/T_end)² once d\* > dy |
| cycles on wave B | 1 | 1 (2 is a separate option that halves the PGA) |

Everything else unchanged: driver, spectra, saves, 0.35 s + tail, cost.

**Carrier option (Ihsan's period).** Wave B can be at the record's period with its amplitude still solved for d\* at the wall's T\*. Single wave, amplitude solved so Sd(T\*) = d\*, FR76 runs 22–24 (T\* = 0.134 s):

| wave B at | PGA / record | PGV / record | Sd at 0.4 s |
|---|---|---|---|
| T\* (wall) | 1.03 | 0.31 | 2.8–4.9 mm |
| T_p = 0.26 s (Sa peak) | 0.93 | 0.53 | 14–25 mm |
| T_m = 0.374 s | 1.37 | 1.12 | 38–67 mm |

The 0.26 s carrier reaches the target below the record's PGA and carries a third to half of the record's long-period demand. On HU12 (run 21, T\* 0.156 s) the Sa-peak (0.12 s) is shorter than T\* and needs 1.47× the record's PGA; T_m (0.206 s) works at 0.85×. Usable rule: carrier period between T\* and ~2 T\*, taken from the record. Implement as `CARRIER = "wall" | "record"`.

**Name.** Capacity-adaptive spectral matching: the record's spectrum is the target; the capacity curve chooses where it is read; nothing from the wall's dynamic response feeds back. "Adaptive sinusoidal spectral matching" remains accurate; "capacity-" says where the adaptation comes from.

---

## 9. Next steps, in order

1. **Step 2 — write the refined `build_signal`** as a drop-in for `route2_full25_revised_3dec.py` with switches `REFINED`, `CARRIER`, `N_CYC_B`, `TIPPING_RULE`. Inputs it needs are all in `data/`.
2. **Step 3 — FR76 22–25 from `route2_run_21.sav`**, three signals: Route 2 as run (known: 8.1 / 28.6 / 32.1 / 38.6), Route 2 refined wall-carrier, Route 2 refined record-carrier. Target line: record on the model 3.7 / 5.9 / 8.3 / 41.6; US-1 and US-2 alongside. Twelve 3 s runs. Plus two cross runs to settle state vs signal: the FR76 record from `route2_run_21.sav`, and the Route 2 run-22 pulse from the control's `stratC_run_21.sav`. Postprocess with the like-for-like base shear (call `instrument_baseshear_exp.dat`; export lines in its header).
3. If step 3 lands near the record-on-model line: full 25-run sequence with the refined signal, SS no-damp first, then Maxwell 1.5% (`route2_full25_cases_3dec.py` has the four damping/geometry cases).
4. Parked: the ×2.317 / density question for Moshfeghi; the tiltmeter comparison (model `tilt_local_incl` against the tiltmeter, with the segment extent confirmed — the tiltmeter reads the rotation of the top courses above the slab-2 joist, not the wall's chord; the deck still says "tilt agrees" and must be corrected); joist slip and loop asymmetry; EUCENTRE kn = 0.36 build; the deck update.

---

## 10. Files in this folder

**scripts/**
- `route2_full25_revised_3dec.py` — the Route 2 driver as run, plus `RESUME_AFTER` / `RESUME_SAVE` for partial re-runs into a new folder and the base-shear instrumentation called at setup and after every `history delete`.
- `instrument_baseshear_exp.dat` — FISH: velocities at the three accelerometer points; export lines in the header.
- `baseshear_exp_from_model.py` — like-for-like loops and envelope slopes; `EXP_DIR` env var points at `TestXXRunYY_processed_globalzero.xlsx`.
- `casm_step1_targets.py` — the fixed-point targets, tipping points, score table and figure 04.
- `fig_synthetic_vs_route2.py` — the spectrum-compatible synthetic and figure 05.

**data/**
- `route2_full_log.csv` — per-run signal parameters and peaks from the driver.
- `route2_full_metrics.csv` — per-run peaks and residuals for US-1, US-2, record-on-model, Route 2.
- `casm_step1_targets.csv` — d\*, T\*, tipping, runaway flag, Route 2 targets, per run.
- `bilinear_periods.json`, `pushover_pos.csv`, `spectrum_*.csv` — the wall and record inputs.
- `synthetic_FR76_x1p00_vel.txt` — 3DEC velocity table of the synthetic (§7).

**Experimental processing conventions** (reverse-engineered exactly from `Test9Run1_processed_globalzero.xlsx`): U_avg = ½(ch3+ch4) − Z0, ch5 not subtracted, Z0 = value at t0 of run 1; base_shear_N = 9.81·[(m15+m16+m17)·ch12 + m15·ch15 + m16·ch16 + m17·ch17] with m = 946.679 / 743.092 / 936.500 kg (= tributary wall mass × 2.317348); `make_processed_globalzero.py` rebuilds these files and was used for Test 12 (US-2). Sign convention in the postprocess: `EXP_SHEAR_SIGN = MODEL_SHEAR_SIGN = −1`.

**On the user's machine, not here:** the save files (`route2_full25_revised/route2_run_NN.sav`, `stratF_full_results_US1/stratC_run_NN.sav`), the channel exports, `EXPRAW/`, `instrument_history_new.dat`, `instrument_history_export_v2.dat`, `instrument_tilt_v2.dat`, `route2_full25_cases_3dec.py`, `route1_extract_3dec.py`, `make_processed_globalzero.py`.
