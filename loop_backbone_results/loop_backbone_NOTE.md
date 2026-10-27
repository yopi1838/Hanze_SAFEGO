# Loop backbone, mass check, test-on-test fixed point — addendum to the 20 Oct handover (§11 next steps)

Scripts: `loop_backbone.py`, `fixed_point_exp.py`. Results: `loop_backbone_results/` (v2 = three-term masses, v2s = five-term).

## Data findings
- **Run 18 DAQ slip.** In `Test9Run18.xlsx` channels 1–7 (wall pots, table pot, overload cells) lead channels 8–17 (joist pots, all accelerometers incl. ch12) by 266–275 samples (≈1.36 s). Run 19 is in sync. Every other run has U vs Σma within ±1 sample (run 24: −3). `loop_backbone.py` detects the slip and shifts U by 272 samples; after the shift run 18 sits between runs 17 and 19. Peaks (Table 4, Fig. 13) are unaffected; anything comparing run-18 time histories against model or table time is 1.36 s off.
- **Potentiometer stiction below ~0.1 mm.** In runs 1–6 the force rings at ~10 Hz after the pulse while U stays flat or steps. The EC40-type runs (3, 6, 8, 11) are mostly in this range. Elastic stiffness is taken only from runs whose smaller peak exceeds 0.1 mm (1, 2, 6).
- **Three- vs five-term base shear** differ by 1–8 % (largest in run 24). Joist motion relative to the wall doesn't drive the force numbers.

## Experimental backbone (simultaneous pairs)
- The experimental backbone has no linear range. The secant at the zero-velocity points falls from ~7.5 kN/mm at 0.3 mm to ~4 at 1 mm, ~2.3 at 3 mm and ~0.3 at 16 mm.
- Plateau ≈ 8.5–9 kN reached at ~3 mm in both directions (runs 16–22). The FR76 runs reach 10.5–11.5 kN at 10–25 mm (peak 11.5 kN, run 23, +11 mm).
- **Strong direction asymmetry.** Up to run 22 the + side stays at ≤ 3.6 mm while F reaches 8–9 kN. The − side ratchets from −4.8 (run 16) to −16.4 mm (run 21), and F at the zero-velocity point drops to −5 kN (rocking softening past the peak). Elastic secants: + 9–14, − 5.8–6.7 kN/mm. The negative-branch pushover in the queue is needed.

## Mass check
| source | m_eff |
|---|---|
| Fig. 13 envelope slope 3.2 kN/mm × (0.091/2π)² | 672 kg (envelope artefact) |
| loop pp-secant, runs 1, 2, 6 (7.5, 8.4, 6.8 kN/mm) × (T_Table5/2π)² | 1580, 1800, 1530 → **≈1.6 t** |
| + side / − side secants only | ~2.1 t / ~1.35 t |
| model K_el·(T1/2π)² | 1406 kg |
| v2 Σ m_j / wall | 1732 / 1340 kg |

With the 672 kg value gone, the test's mass agrees with the model's and the masses. The 13 % period gap (0.091 vs 0.103 s) is therefore a stiffness gap: the loop's 7.5–8.4 kN/mm vs the model's 5.1, ≈ 1.5×. m_eff between 1406 and 1732 kg moves d* by < 10 % on runs 10–21. The mass question no longer decides the backbone.

## Fixed point on the test's own backbone and the achieved table motion (ch12, 5 %)
- **HU12 runs 1–14:** d* over-predicts by 1.4–2.4× (amplitude-dependent damping above 5 %; hysteretic loops visible from run 4).
- **HU12 runs 15–17:** d*/measured 1.27 → 1.12. **Runs 18–21:** 0.89 → 0.64 on the max. d* follows the half peak-to-peak (18: 7.6 vs 5.6; 21: 10.5 vs 9.8); what it misses is the one-sided negative drift, i.e. cycle accumulation.
- **FR76 runs 22–24:** d* = 29 / 57 / 61 mm vs measured 3.6 / 13 / 29. The achieved FR76 table motion has Sa ≈ 0.7–1.25 g flat over 0.15–0.5 s, above the capacity acceleration F/m ≈ 0.55–0.64 g. On a flat backbone in a flat-Sa band, G(d)/d = Sa/a_cap > 1 at every d, so no crossing exists until the spectrum drops (~0.55 s, ~55–60 mm). The fixed point degenerates into a strength check. This is the run-22 ill-conditioning of §7 in general form: **d* is defined only where the record's local slope dlnSd/dlnT < 2 at T*, or where Sa < a_cap.** That's the regime N2 handles with an R–μ–T rule rather than equal displacement.

## Consequence for the surrogate
The capacity-adaptive target works, on the test itself, through the HU12 build-up to the first strong rocking runs (+27 / −36 % on runs 15–21). It cannot carry FR76 on a rocking backbone, independently of the DEM model. For the FR76 runs, use the extraction route, or an R–μ–T/strength-ratio target with an explicit tipping check.

## Still open
v5 pushover CSV not found on disk (no `pushover_oldnew_summary.csv` / `m2fp_targets_preflight_v5backbone.csv` in the folder), so the model-backbone overlay (`--pushover`) and the two-backbone d* spread are pending. No output folders for `F_LS_MAXWELL15_v12` or `route2_m2_fixedpoint_v13_from21_nodamp` in the folder yet. `vel_EC.txt` *is* in the folder (handover §9 said it wasn't available).
