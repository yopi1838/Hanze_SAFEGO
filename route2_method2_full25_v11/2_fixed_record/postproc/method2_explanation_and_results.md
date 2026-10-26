# Route 2 Method 2: explanation and current results

## Purpose

This note explains Method 2 and compares the earlier four-run FR76 test, which started from the Strategy F run-21 state, with the later uninterrupted Method 2 sequence from runs 1–25.

## What Method 2 does

Method 2 builds a short, two-part table-velocity input for each IDA run:

- **Wave A** is one sinusoidal cycle carried at the fixed wall period, `T1`. Its amplitude is fitted so the input has the run-specific target spectral displacement `s × Sd(T1)`.
- **Wave B** targets `s × Sd(Teff)`, where `Teff` is fixed from the pushover idealisation. Its carrier period comes from the earthquake record's predominant-period rule, `Tp`. For FR76, `Tp = 0.260 s`, while the spectral target is evaluated at `Teff = 0.24138 s`.
- The amplitudes of A and B are fitted **jointly on the combined input**. The scale `s` is the prescribed IDA intensity for that run.
- In the full-25 version, the B-wave cycle count and relative polarity are selected from a fixed candidate order using the record and the spectral targets. The selection does not use the wall's response from the previous run. Therefore Method 2 is record-informed, but it is **not response-adaptive**.
- The full-25 sequence starts from the initial small-strain model save and carries the model's damage from one run to the next.

The spectral fit is a check on two 5%-damped linear SDOF ordinates. It is not a fit to the nonlinear wall displacement, cycle-by-cycle damage, energy input, or shaking duration.

## Stage 1: four FR76 runs starting from the Strategy F run-21 state

The earlier four-method comparison restored a common Strategy F run-21 state, then carried each method's own state through FR76 runs 22–25. For Method 2, the results were:

| Run | Method 2 peak OOP (mm) | US-1 test (mm) | US-2 test (mm) |
|---:|---:|---:|---:|
| 22 | 2.47 | 3.65 | 2.38 |
| 23 | 12.10 | 13.15 | 5.26 |
| 24 | 32.75 | 29.48 | 7.58 |
| 25 | 41.01 | — | 31.14 |

This run-21-started sequence produced strong growth through the later FR76 runs. The run-22 input log, however, shows that this was an earlier implementation: one B cycle, 0.37 s input duration, about 0.331 g combined PGA, and a non-converged joint fit (+16.70% at `T1`, −0.71% at `Teff`). It is useful evidence that the starting model state matters, but it is **not a controlled comparison** with the later full-25 Method 2 inputs.

## Stage 2: uninterrupted Method 2 sequence from run 1 through run 25

The v11 sequence began from the initial model state, ran the prescribed HU12/EC40 sequence through run 21, then carried its own resulting state into FR76. The checkpoint and run log agree that all 25 runs completed. For FR76, v11 selected two same-phase B cycles at 0.260 s. The amplitude fit converged, with negligible error at both target periods. Active input duration was 0.63 s.

| Run | Method 2 peak OOP (mm) | Strategy F control (mm) | US-1 test (mm) | US-2 test (mm) |
|---:|---:|---:|---:|---:|
| 22 | 2.54 | 3.72 | 3.65 | 2.38 |
| 23 | 3.00 | 5.89 | 13.15 | 5.26 |
| 24 | 3.16 | 8.29 | 29.48 | 7.58 |
| 25 | 3.74 | 41.63 | — | 31.14 |

The full-25 model stays at roughly 2.5–3.7 mm through FR76, well below the Strategy F control and the later test peaks. Run 22 is close to US-2; the large mismatch develops over runs 23–25.

Crucially, the divergence is already present **before FR76**:

| Run | Record | Method 2 peak OOP (mm) | Strategy F control (mm) | US-1 test (mm) |
|---:|---|---:|---:|---:|
| 18 | HU12 × 4.5 | 3.92 | 3.97 | 8.48 |
| 19 | HU12 × 5.0 | 5.50 | 9.29 | 10.74 |
| 20 | HU12 × 5.5 | 5.91 | 13.73 | 13.52 |
| 21 | HU12 × 6.0 | 6.00 | 16.97 | 16.37 |

The same-model control and Method 2 are close at run 18, then separate at run 19 and diverge further at runs 20–21. Method 2 therefore reaches FR76 from a less displaced, likely less damaged state than the Strategy F control. This is a plausible major reason that the later FR76 response remains low.

## What the comparison means

The two sequences changed **both** the starting state and the waveform implementation:

1. The four-run test starts at the Strategy F run-21 state; the full-25 run builds its own history from the initial model save. By run 21, that Method 2 history already underpredicts the control response.
2. The four-run Method 2 input used one B cycle and did not meet the `T1` target within the solver tolerance. The v11 full-sequence input uses two B cycles for FR76 and achieves both targets essentially exactly.

Consequently, the high FR76 response in the four-run test does not prove that the v11 full sequence should produce the same response, and the difference cannot be attributed to starting state alone. The clean conclusion is that matching `Sd` at `T1` and `Teff` is not sufficient to reproduce the response of the wall over this sequence. The mismatch begins during high-intensity HU12 runs and becomes much larger during FR76.

The short surrogate waveforms may not reproduce important features of the original records, such as the sequence and number of strong cycles that accumulate damage. That is a working explanation, not yet a demonstrated cause. A controlled same-state comparison is needed before changing the waveform rule or making a predictive-performance claim.

## Recommended diagnostic

1. Start two branches from the **same saved state immediately before run 19**. Apply the original HU12 run-19 input in one branch and the exact v11 Method 2 run-19 input in the other. Compare peak and residual OOP displacement and available damage/energy indicators.
2. To isolate the FR76 effect, apply the exact v11 FR76 inputs for runs 22–25 to a common Strategy F run-21 state, then compare with the v11 continuation from its own run-21 state.
3. Preserve the model saves and exact velocity tables for both branches. The JSON checkpoint is a protocol/results summary, not a 3DEC restart state; the Method 2 save-file paths recorded in the log are not present in this exported output folder.

The existing input-fit plot and sequence plot are included for reference: [two-period input fit](fig_method2_full25_input_fit.png) and [full-sequence comparison](fig_method2_full25_sequence.png). Base-shear comparison is not used here because the mass conversion convention is still being investigated.

## Data references

- [Full-25 Method 2 run log](../route2_method_log.csv)
- [Full-25 Method 2 checkpoint summary](../route2_method2_checkpoint.json)
- [Full-25 postprocessed metrics](route2_method2_full25_metrics.csv)
- [Earlier four-method postprocessed metrics](../../../route2_four_methods/postproc/route2_four_method_metrics.csv)
- [Earlier four-run Method 2 run log](../../../route2_four_methods/2_fixed_record/route2_method_log.csv)
