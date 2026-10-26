import fs from "node:fs/promises";
import path from "node:path";
import { execFileSync } from "node:child_process";
import { pathToFileURL } from "node:url";
import { Presentation, PresentationFile } from "@oai/artifact-tool";

const workspaceDir = "C:/Users/yopi1/Documents/09.Research_Dataset/Hanze_SAFEGO";
const buildDir = path.join(workspaceDir, ".presentation_build_us2_postprocess_v3");
const skillDir = "C:/Users/yopi1/.codex/plugins/cache/openai-primary-runtime/presentations/26.909.12148/skills/Presentations";
const pythonExecutable = "C:/Users/yopi1/.cache/codex-runtimes/codex-primary-runtime/dependencies/python/python.exe";
const finalPath = path.join(workspaceDir, "00.Meeting", "Progress_meeting_Route2_four_methods_2026-10-05_with_US2_postprocess_v3.pptx");
const sourceDeck = path.join(workspaceDir, "00.Meeting", "Progress_meeting_Hanze_SAFEGO.pptx");
const figureDir = path.join(workspaceDir, "route2_four_methods");
const C = {
  charcoal: "#2C2927",
  ink: "#272523",
  muted: "#716D69",
  cream: "#F1EFEC",
  terracotta: "#B7684B",
  terracottaLight: "#D9B5A5",
  white: "#FFFFFF",
  blue: "#276B8E",
  orange: "#E07A2D",
  green: "#3E8C61",
  purple: "#A14D8D",
  rule: "#D4D0CC",
};
const FONT_TITLE = "Cambria";
const FONT_BODY = "Calibri";

const { resolvePresentationFont, applyPresentationChartFont, finalizePresentation } = await import(
  pathToFileURL(path.join(skillDir, "container_tools/artifact_tool_utils.mjs")).href,
);
await fs.mkdir(buildDir, { recursive: true });
const familyCheck = resolvePresentationFont({ fontFamily: FONT_BODY });
if (!familyCheck) throw new Error(`Required font unavailable: ${FONT_BODY}`);

const p = Presentation.create({ slideSize: { width: 1280, height: 720 } });

function shape(slide, geometry, name, position, fill = "none", lineFill = "none", lineWidth = 0) {
  return slide.shapes.add({
    geometry,
    name,
    position,
    fill,
    line: { style: "solid", fill: lineFill, width: lineWidth },
  });
}

function text(slide, name, value, position, opts = {}) {
  const s = shape(slide, "textbox", name, position);
  s.text = value;
  s.text.style = {
    typeface: opts.font ?? FONT_BODY,
    fontSize: opts.size ?? 20,
    color: opts.color ?? C.ink,
    bold: opts.bold ?? false,
    italic: opts.italic ?? false,
    alignment: opts.align ?? "left",
    verticalAlignment: opts.valign ?? "top",
    autoFit: "shrinkText",
    wrap: "square",
    insets: { top: 0, right: 0, bottom: 0, left: 0 },
    ...(opts.lineSpacing ? { lineSpacing: opts.lineSpacing } : {}),
  };
  return s;
}

function addBodyHeader(slide, slideNumber, titleText, subtitle = "") {
  shape(slide, "ellipse", `slide-number-circle-${slideNumber}`, { left: 52, top: 37, width: 42, height: 42 }, C.terracotta);
  text(slide, `slide-number-${slideNumber}`, String(slideNumber), { left: 52, top: 37, width: 42, height: 42 }, {
    size: 21, bold: true, color: C.white, align: "center", valign: "middle",
  });
  text(slide, `slide-title-${slideNumber}`, titleText, { left: 112, top: 30, width: 1110, height: 56 }, {
    font: FONT_TITLE, size: 36, bold: true, color: C.ink, valign: "middle",
  });
  if (subtitle) {
    text(slide, `slide-subtitle-${slideNumber}`, subtitle, { left: 112, top: 88, width: 1100, height: 30 }, {
      size: 18, italic: true, color: C.muted, valign: "middle",
    });
  }
}

function addWave(slide, name, left, top, width, height, cycles, color) {
  const commands = [];
  const n = 60;
  for (let i = 0; i <= n; i++) {
    const x = width * i / n;
    const y = height / 2 - Math.sin(2 * Math.PI * cycles * i / n) * height * 0.38;
    commands.push(i === 0 ? { moveTo: { x, y } } : { lineTo: { x, y } });
  }
  return slide.shapes.add({
    geometry: "custom",
    name,
    position: { left, top, width, height },
    fill: "none",
    line: { style: "solid", fill: color, width: 2.4 },
    customPaths: [{ width, height, commands }],
  });
}

function addRule(slide, name, x1, y1, x2, y2, color = C.rule, width = 1) {
  return slide.shapes.add({
    geometry: "line",
    name,
    position: { left: x1, top: y1, width: x2 - x1, height: y2 - y1 },
    fill: "none",
    line: { style: "solid", fill: color, width },
  });
}

async function addPng(slide, name, filePath, position, alt) {
  const bytes = new Uint8Array(await fs.readFile(filePath));
  return slide.images.add({
    blob: bytes,
    contentType: "image/png",
    alt,
    fit: "contain",
    position,
  });
}

function setNotes(slide, note) {
  slide.speakerNotes.textFrame.setText(note);
}

function parseCsvRows(csvText) {
  const lines = csvText.trim().split(/\r?\n/);
  const headers = lines[0].split(",");
  return lines.slice(1).filter(Boolean).map((line) => {
    const cells = line.split(",");
    return Object.fromEntries(headers.map((header, i) => [header, cells[i]]));
  });
}

function chartPrecision(value) {
  return Number(Number(value).toPrecision(12));
}

// 1. Cover, carrying forward the previous deck's charcoal field and terracotta accent.
{
  const slide = p.slides.add();
  slide.background.fill = C.charcoal;
  text(slide, "cover-title", "Adaptive sinusoidal spectral matching\nfor out-of-plane masonry walls", {
    left: 84, top: 177, width: 1030, height: 174,
  }, { font: FONT_TITLE, size: 50, bold: true, color: C.white, lineSpacing: 0.95 });
  text(slide, "cover-subtitle", "Four Route 2 methods tested from the common Strategy F Run 21 state", {
    left: 88, top: 366, width: 1000, height: 60,
  }, { size: 25, color: "#D7D2CE", valign: "middle" });
  text(slide, "cover-meta", "Progress meeting, Hanze SAFEGo US-1 + US-2, 5 October 2026", {
    left: 88, top: 538, width: 900, height: 42,
  }, { size: 18, color: "#C4BFBB", valign: "middle" });
  shape(slide, "ellipse", "cover-terracotta-circle", { left: 1075, top: 454, width: 205, height: 205 }, C.terracotta);
  setNotes(slide, "This update focuses on the four Route 2 carrier and target combinations. Each sequence starts from the same Strategy F Run 21 save and carries the damage state forward through FR76 Runs 22–25.");
}

// 2. One editable figure slide to explain all four method definitions.
{
  const slide = p.slides.add();
  slide.background.fill = C.white;
  addBodyHeader(slide, 1, "Four Route 2 methods", "Two target rules crossed with two Wave B carrier periods");

  const methods = [
    { n: 1, color: C.blue, x: 62, y: 137, target: "B target: scale × Sd(Teff)", carrier: "B carrier: Teff = 0.241 s", carrierPeriod: 0.241 },
    { n: 2, color: C.orange, x: 656, y: 137, target: "B target: scale × Sd(Teff)", carrier: "B carrier: record Tp = 0.260 s", carrierPeriod: 0.260 },
    { n: 3, color: C.green, x: 62, y: 369, target: "B target: (T*, d*) by intensity", carrier: "B carrier: T* = 0.134 s", carrierPeriod: 0.134 },
    { n: 4, color: C.purple, x: 656, y: 369, target: "B target: (T*, d*) by intensity", carrier: "B carrier: record Tp = 0.260 s", carrierPeriod: 0.260 },
  ];
  for (const m of methods) {
    const w = 562, h = 210;
    shape(slide, "rect", `method-${m.n}-surface`, { left: m.x, top: m.y, width: w, height: h }, C.cream);
    shape(slide, "rect", `method-${m.n}-accent`, { left: m.x, top: m.y, width: 7, height: h }, m.color);
    text(slide, `method-${m.n}-title`, `Method ${m.n}`, { left: m.x + 24, top: m.y + 12, width: 165, height: 34 }, {
      font: FONT_TITLE, size: 26, bold: true, color: m.color, valign: "middle",
    });
    text(slide, `method-${m.n}-target`, m.target, { left: m.x + 202, top: m.y + 13, width: 336, height: 30 }, {
      size: 17, bold: true, color: C.ink, valign: "middle",
    });
    text(slide, `method-${m.n}-carrier`, m.carrier, { left: m.x + 202, top: m.y + 44, width: 336, height: 28 }, {
      size: 16, color: C.muted, valign: "middle",
    });
    const baseY = m.y + 140;
    const waveW = 190, waveH = 34;
    addRule(slide, `method-${m.n}-axis-a`, m.x + 24, baseY, m.x + 214, baseY, "#AAA49E", 1);
    addRule(slide, `method-${m.n}-axis-b`, m.x + 334, baseY, m.x + 524, baseY, "#AAA49E", 1);
    addWave(slide, `method-${m.n}-wave-a`, m.x + 24, baseY - waveH / 2, waveW, waveH, 3, C.terracotta);
    addWave(slide, `method-${m.n}-wave-b`, m.x + 334, baseY - waveH / 2, waveW, waveH, 0.336 / m.carrierPeriod, m.color);
    text(slide, `method-${m.n}-wave-a-label`, "Wave A · T1 = 0.112 s", { left: m.x + 24, top: m.y + 158, width: 210, height: 26 }, {
      size: 15, color: C.ink, valign: "middle",
    });
    text(slide, `method-${m.n}-wave-b-label`, "Wave B · carrier shown above", { left: m.x + 334, top: m.y + 158, width: 214, height: 26 }, {
      size: 15, color: C.ink, valign: "middle",
    });
  }
  addRule(slide, "methods-footer-rule", 62, 600, 1218, 600, C.rule, 1);
  text(slide, "methods-footer-note", "Target = desired spectral displacement at an oscillator period. Carrier = input sine period. The fit adjusts wave velocities with carrier periods fixed; Methods 3–4 use precomputed targets.", {
    left: 64, top: 612, width: 1150, height: 58,
  }, { size: 15, color: C.muted, italic: true, valign: "middle" });
  setNotes(slide, "Wave A is common to all four methods and targets scale × Sd(T1), where T1 is about 0.1119 s. Methods 1 and 2 use fixed target points at T1 and Teff, with Teff about 0.2414 s. Methods 3 and 4 use the precomputed capacity-proxy target d* and T* for each intensity. Method 1 carries Wave B at Teff; Method 2 uses the FR76 predominant period Tp = 0.260 s; Method 3 carries Wave B at T*; Method 4 uses Tp. The amplitude fit changes the two input velocity amplitudes while the carrier periods remain prescribed. The adaptive schedule is precomputed and independent of the preceding 3DEC response.");
}

// 3. The 25 run-scaled spectra show how the fixed periods supply Route 2 targets.
{
  const slide = p.slides.add();
  slide.background.fill = C.white;
  addBodyHeader(slide, 2, "Reading targets at fixed periods", "T1 and Teff come from the wall idealisation; each run supplies a scaled record spectrum");

  const periodLimit = 0.35;
  const targetRows = parseCsvRows(await fs.readFile(path.join(workspaceDir, "Claude outputs", "casm_step1_targets.csv"), "utf8"));
  const periodData = JSON.parse(await fs.readFile(path.join(workspaceDir, "bilinear_idealisation", "bilinear_periods.json"), "utf8"));
  const T1 = periodData.model.T1_s;
  const TEFF = periodData.model.Teff_s;
  const groupColors = { HU12: C.blue, EC40: C.green, FR76: C.terracotta };
  const spectrumByRecord = {};
  for (const record of ["HU12", "EC40", "FR76"]) {
    const spectrumText = await fs.readFile(path.join(workspaceDir, `spectrum_${record}.csv`), "utf8");
    spectrumByRecord[record] = parseCsvRows(spectrumText)
      .map((r) => ({ period: Number(r.T_s), sdMm: Number(r.Sd_m) * 1000 }))
      .filter((r) => r.period <= periodLimit);
  }
  const maxDemand = Math.max(...targetRows.flatMap((run) =>
    spectrumByRecord[run.record].map((point) => point.sdMm * Number(run.scale))));
  const spectrumAxisMax = Math.ceil((maxDemand * 1.05) / 10) * 10;
  const spectrumSeries = targetRows.map((run) => ({
    name: `Run ${run.run}`,
    xValues: spectrumByRecord[run.record].map((point) => chartPrecision(point.period)),
    values: spectrumByRecord[run.record].map((point) => chartPrecision(point.sdMm * Number(run.scale))),
    line: { style: "solid", fill: groupColors[run.record], width: 1.15 },
  }));
  spectrumSeries.push(
    { name: "T1 = 0.112 s", xValues: [chartPrecision(T1), chartPrecision(T1)], values: [0, spectrumAxisMax], line: { style: "solid", fill: C.ink, width: 1.7 } },
    { name: "Teff = 0.241 s", xValues: [chartPrecision(TEFF), chartPrecision(TEFF)], values: [0, spectrumAxisMax], line: { style: "solid", fill: C.terracotta, width: 2.2 } },
  );
  const spectrumChart = slide.charts.add("scatter", {
    position: { left: 62, top: 148, width: 1156, height: 420 },
    series: spectrumSeries,
    scatterOptions: { style: "line" },
    hasLegend: false,
    chartFill: C.white,
    chartLine: { style: "solid", fill: "none", width: 0 },
    plotAreaFill: C.white,
    plotAreaLine: { style: "solid", fill: "none", width: 0 },
    xAxis: {
      title: { text: "Oscillator period T (s) — target range", textStyle: { typeface: FONT_BODY, fontSize: 16, fill: C.ink } },
      min: 0, max: periodLimit, majorUnit: 0.05, numberFormatCode: "0.00",
      textStyle: { typeface: FONT_BODY, fontSize: 13, fill: C.muted },
      line: { style: "solid", fill: "#918B85", width: 1 },
      majorGridlines: null,
    },
    yAxis: {
      title: { text: "Scaled spectral displacement Sd (mm)", textStyle: { typeface: FONT_BODY, fontSize: 16, fill: C.ink } },
      min: 0, max: spectrumAxisMax, majorUnit: 10, numberFormatCode: "0",
      textStyle: { typeface: FONT_BODY, fontSize: 13, fill: C.muted },
      line: { style: "solid", fill: "#918B85", width: 1 },
      majorGridlines: { style: "solid", fill: "#DDD9D4", width: 1 },
    },
  });
  applyPresentationChartFont(spectrumChart, { fontFamily: FONT_BODY });
  text(slide, "spectrum-period-key", "Guides: T1 = 0.112 s (dark)  ·  Teff = 0.241 s (terracotta)", {
    left: 80, top: 570, width: 680, height: 28,
  }, { size: 15, color: C.ink, bold: true, valign: "middle" });
  text(slide, "spectrum-record-key", "Curves: HU12 (blue, 21 runs)  ·  EC40 (green, 4 runs)  ·  FR76 (terracotta, 4 runs)", {
    left: 80, top: 599, width: 1120, height: 25,
  }, { size: 14, color: C.muted, valign: "middle" });
  text(slide, "spectrum-target-equations", "Wave A target = scale × Sd(T1)     ·     Method 1 Wave B target = scale × Sd(Teff)", {
    left: 80, top: 638, width: 1120, height: 33,
  }, { size: 17, color: C.terracotta, bold: true, valign: "middle" });
  setNotes(slide, "All 25 lines are the 5%-damped displacement spectra for the run's ground-motion record multiplied by that run's protocol scale. T1 = 0.111887 s and Teff = 0.241382 s are the fixed periods from the wall's bilinear pushover idealisation and are reused for every run. For each run, the ordinate read at T1 sets the Wave A spectral-displacement target, and the ordinate read at Teff sets Method 1's Wave B target. The chart zooms the period axis to 0–0.35 s to make both fixed target periods legible. The lines are grouped by record; repeated curves have different protocol scales.");
}

// 4. Method 3 Run 22's capacity-proxy target is a numerical fixed point.
{
  const slide = p.slides.add();
  slide.background.fill = C.white;
  addBodyHeader(slide, 3, "Method 3 target for Run 22", "A stable fixed point between the scaled FR76 spectrum and the cracked pushover backbone");
  text(slide, "dstar-equation", "d* = 1.00 × Sd,FR76(Tsec(d*))", { left: 80, top: 145, width: 1120, height: 54 }, {
    font: FONT_TITLE, size: 31, bold: true, color: C.ink, align: "center", valign: "middle",
  });
  text(slide, "dstar-period-equation", "Tsec(d) = 2π √[ meff / (Fcracked(d) / d) ]", { left: 80, top: 200, width: 1120, height: 38 }, {
    font: FONT_TITLE, size: 22, color: C.muted, align: "center", valign: "middle",
  });

  const calcColumns = [
    { x: 82, heading: "TRIAL DISPLACEMENT", value: "3.423 mm", detail: "Run 22 · FR76 × scale 1.00" },
    { x: 478, heading: "PUSHOVER SECANT", value: "12.360 kN", detail: "Fcracked at d* · ksec = F/d = 3.611 kN/mm" },
    { x: 874, heading: "READ THE SPECTRUM", value: "0.133704 s", detail: "Tsec → Sd,FR76 = 3.423 mm" },
  ];
  for (const [i, col] of calcColumns.entries()) {
    text(slide, `dstar-step-${i}-heading`, col.heading, { left: col.x, top: 302, width: 320, height: 28 }, {
      size: 15, bold: true, color: C.terracotta, valign: "middle",
    });
    text(slide, `dstar-step-${i}-value`, col.value, { left: col.x, top: 343, width: 330, height: 53 }, {
      font: FONT_TITLE, size: 31, bold: true, color: i === 2 ? C.green : C.ink, valign: "middle",
    });
    text(slide, `dstar-step-${i}-detail`, col.detail, { left: col.x, top: 403, width: 330, height: 69 }, {
      size: 16, color: C.muted, valign: "top",
    });
    if (i < 2) addRule(slide, `dstar-step-divider-${i}`, col.x + 354, 303, col.x + 354, 472, C.rule, 1);
  }
  text(slide, "dstar-mass-check", "meff = 1,635 kg  →  Tsec = 2π√(1,635 / 3.611×10⁶) = 0.133704 s", {
    left: 84, top: 490, width: 1110, height: 36,
  }, { size: 17, color: C.ink, valign: "middle" });
  addRule(slide, "dstar-result-rule", 82, 545, 1198, 545, C.rule, 1);
  text(slide, "dstar-result-note", "At that period, FR76 gives Sd = 3.423 mm. Demand equals the trial displacement, so this is d*.", {
    left: 84, top: 562, width: 1110, height: 40,
  }, { font: FONT_TITLE, size: 21, bold: true, color: C.green, valign: "middle" });
  text(slide, "dstar-schedule-caveat", "The cracked backbone follows the precomputed target schedule after an earlier target crossed dy; Run 21's measured response is not used.", {
    left: 84, top: 623, width: 1110, height: 40,
  }, { size: 15, color: C.muted, italic: true, valign: "middle" });
  setNotes(slide, "The value d*=3.423165 mm is read from Claude outputs/casm_step1_targets.csv for Run 22, FR76 at scale 1.00. The calculation uses the actual pushover force-displacement backbone and the 5%-damped FR76 spectrum. At d=3.423165 mm the cracked backbone gives F=12.360003 kN, hence F/d=3.610694 kN/mm. With m_eff=1635 kg, the secant-period calculation gives T*=0.133703694 s. Linear interpolation of the FR76 spectrum at that period gives Sd=3.423165 mm, so scale × Sd(Tsec(d)) equals d. The script numerically scans and interpolates the first stable crossing. The target schedule's cracked-backbone state is property-based damage memory: it switches after an earlier scheduled d* crossed dy; it does not read the measured Run 21 or Run 22 3DEC response. Method 3 uses this T* as its Wave B carrier before snapping to the 0.135 s time-grid period.");
}

// 5. The saved Run 22 velocity tables show the actual four 3DEC inputs.
{
  const slide = p.slides.add();
  slide.background.fill = C.white;
  addBodyHeader(slide, 4, "Run 22 input velocity", "Actual FR76 × 1.00 velocity tables imported into 3DEC · Wave A is common to all methods");
  const inputMethods = [
    { key: "1_fixed_wall", label: "Method 1 · T_B = 0.240 s", color: C.blue },
    { key: "2_fixed_record", label: "Method 2 · T_B = 0.260 s", color: C.orange },
    { key: "3_adaptive_wall", label: "Method 3 · T_B = 0.135 s", color: C.green },
    { key: "4_adaptive_record", label: "Method 4 · T_B = 0.260 s", color: C.purple },
  ];
  const inputSeries = [];
  for (const method of inputMethods) {
    const inputPath = path.join(figureDir, method.key, `Run22_FR76_s1p00_${method.key}_vel.txt`);
    const lines = (await fs.readFile(inputPath, "utf8")).trim().split(/\r?\n/).slice(2);
    const samples = lines.map((line) => line.trim().split(/\s+/).map(Number));
    inputSeries.push({
      name: method.label,
      xValues: samples.map((r) => r[0]),
      values: samples.map((r) => r[1]),
      line: { style: "solid", fill: method.color, width: 2.1 },
    });
  }
  const inputChart = slide.charts.add("scatter", {
    position: { left: 58, top: 144, width: 1162, height: 452 },
    series: inputSeries,
    scatterOptions: { style: "line" },
    hasLegend: true,
    legend: { position: "bottom", overlay: false, textStyle: { typeface: FONT_BODY, fontSize: 14, fill: C.ink } },
    chartFill: C.white,
    chartLine: { style: "solid", fill: "none", width: 0 },
    plotAreaFill: C.white,
    plotAreaLine: { style: "solid", fill: "none", width: 0 },
    xAxis: {
      title: { text: "Time (s)", textStyle: { typeface: FONT_BODY, fontSize: 16, fill: C.ink } },
      min: 0, max: 0.37, majorUnit: 0.05, numberFormatCode: "0.00",
      textStyle: { typeface: FONT_BODY, fontSize: 13, fill: C.muted },
      line: { style: "solid", fill: "#918B85", width: 1 },
      majorGridlines: null,
    },
    yAxis: {
      title: { text: "Applied table velocity (m/s)", textStyle: { typeface: FONT_BODY, fontSize: 16, fill: C.ink } },
      min: -0.16, max: 0.16, majorUnit: 0.08, numberFormatCode: "0.00",
      textStyle: { typeface: FONT_BODY, fontSize: 13, fill: C.muted },
      line: { style: "solid", fill: "#918B85", width: 1 },
      majorGridlines: { style: "solid", fill: "#DDD9D4", width: 1 },
    },
  });
  applyPresentationChartFont(inputChart, { fontFamily: FONT_BODY });
  text(slide, "input-velocity-note", "Wave A carrier = 0.110 s in every method. The fitting routine changes the Wave A/B velocity amplitudes; the selected carrier periods stay fixed.", {
    left: 78, top: 611, width: 1125, height: 55,
  }, { size: 16, color: C.muted, italic: true, valign: "middle" });
  setNotes(slide, "Each series is read from the corresponding Run22_FR76_s1p00_*_vel.txt table saved in route2_four_methods and imported by the 3DEC Itasca Python script. These are the applied boundary/table velocities at the native 0.005 s sample interval, not wall response velocities. The common Wave A carrier is 0.110 s on the input grid. Wave B grid periods are 0.240 s, 0.260 s, 0.135 s, and 0.260 s for Methods 1–4. The carrier periods are set by each method; the joint spectral fit adjusts the component velocities.");
}

// 6. Peak-response comparison, using the requested shortened labels.
{
  const slide = p.slides.add();
  slide.background.fill = C.white;
  addBodyHeader(slide, 5, "Peak response and spectral calibration", "Runs 22–24 compare with US-1 Test 9; Run 25 has no matching specimen test");
  await addPng(slide, "four-method-comparison-figure", path.join(figureDir, "postproc", "fig_route_2_four_method_comparison.png"), {
    left: 49, top: 130, width: 1182, height: 390,
  }, "Four-method comparison of peak relative out-of-plane displacement, peak base shear, and mean spectral target error. Legends read Method 1 through Method 4.");
  text(slide, "best-displacement-value", "17.1%", { left: 95, top: 540, width: 235, height: 48 }, {
    font: FONT_TITLE, size: 40, bold: true, color: C.orange,
  });
  text(slide, "best-displacement-label", "Method 2 mean OOP error", { left: 98, top: 587, width: 310, height: 30 }, { size: 18, color: C.ink });
  addRule(slide, "comparison-mid-rule", 640, 541, 640, 618, C.rule, 1);
  text(slide, "best-baseshear-value", "9.3%", { left: 700, top: 540, width: 200, height: 48 }, {
    font: FONT_TITLE, size: 40, bold: true, color: C.blue,
  });
  text(slide, "best-baseshear-label", "Method 1 mean base-shear peak error", { left: 704, top: 587, width: 485, height: 30 }, { size: 18, color: C.ink });
  text(slide, "fit-provisional-note", "All 16 spectral fits missed the 2% convergence tolerance. Treat these response rankings as provisional.", {
    left: 74, top: 642, width: 1130, height: 35,
  }, { size: 16, color: C.muted, italic: true, valign: "middle" });
  setNotes(slide, "Metrics use the three experimental comparisons available for US-1 Test 9, Runs 22–24. Method 2 has the lowest mean absolute percentage error in peak relative OOP displacement at 17.1%. Method 1 has the lowest mean absolute percentage error in peak base shear at 9.3%. Run 25 is shown for the model sequence only because no US-1 Test 9 experiment exists for that run. The four methods used the common Strategy F Run 21 starting state. All 16 signal fits report fit_converged=False at the 2% tolerance, so the method ranking is preliminary.");
}

// 7. Repeat the comparison with Methods 1, 3 and 4 faded to focus attention on Method 2.
{
  const slide = p.slides.add();
  slide.background.fill = C.white;
  addBodyHeader(slide, 6, "Peak-response comparison: Method 2 highlighted", "Methods 1, 3 and 4 are muted; the experimental series remains visible");
  await addPng(slide, "method2-focused-comparison-figure", path.join(figureDir, "postproc", "fig_route_2_four_method_comparison_method2_focus.png"), {
    left: 49, top: 130, width: 1182, height: 408,
  }, "The same peak displacement, experiment-equivalent base-shear, and spectral calibration data as the prior slide. Method 2 remains orange; Methods 1, 3, and 4 share a pale grey treatment; US-1 experimental points remain black.");
  text(slide, "method2-highlight-takeaway", "Method 2 has the lowest mean peak OOP error (17.1%); Method 1 still has the lowest base-shear error (9.3%).", {
    left: 76, top: 562, width: 1130, height: 42,
  }, { size: 19, color: C.ink, bold: true, valign: "middle" });
  text(slide, "method2-highlight-caveat", "All 16 spectral fits missed the 2% convergence tolerance; treat these rankings as provisional.", {
    left: 76, top: 637, width: 1130, height: 32,
  }, { size: 16, color: C.muted, italic: true, valign: "middle" });
  setNotes(slide, "This is a visual focus transition from the preceding all-method comparison. Methods 1, 3 and 4 use one pale grey with partial transparency; Method 2 retains its orange emphasis, and the US-1 experiment remains black. The underlying data are unchanged. Method 2 has the lowest mean absolute peak relative OOP displacement error at 17.1%, while Method 1 has the lowest mean absolute peak base-shear error at 9.3%. All 16 spectral fits missed the 2% convergence tolerance, so these rankings remain preliminary.");
}

// 8. Explain the support-reaction channel versus the test's inertia-derived definition.
{
  const slide = p.slides.add();
  slide.background.fill = C.white;
  addBodyHeader(slide, 7, "Two base-shear measures", "cstav is a model support reaction; the test reports an inertia-derived estimate");
  addRule(slide, "baseshear-definition-divider", 640, 148, 640, 492, C.rule, 1);

  text(slide, "cstav-heading", "DIRECT 3DEC OUTPUT", { left: 78, top: 151, width: 500, height: 28 }, {
    size: 15, bold: true, color: C.muted,
  });
  text(slide, "cstav-label", "cstav", { left: 80, top: 190, width: 470, height: 57 }, {
    font: FONT_TITLE, size: 42, bold: true, color: "#8A8F97", valign: "middle",
  });
  text(slide, "cstav-equation", "V_cstav(t) = sign × [cstav(t) − cstav(0)]", {
    left: 82, top: 253, width: 490, height: 50,
  }, { font: FONT_TITLE, size: 23, color: C.ink, valign: "middle" });
  addRule(slide, "cstav-divider", 82, 316, 570, 316, C.rule, 1);
  text(slide, "cstav-meaning-label", "What it measures", { left: 82, top: 337, width: 450, height: 28 }, {
    size: 18, bold: true, color: C.ink,
  });
  text(slide, "cstav-meaning", "The net reaction from the model’s base and joist supports, read directly from the 3DEC history.", {
    left: 82, top: 369, width: 485, height: 62,
  }, { size: 18, color: C.ink });
  text(slide, "cstav-use", "Kept as a grey reference trace. It is not the quantity reported by the experimental processing formula.", {
    left: 82, top: 438, width: 485, height: 54,
  }, { size: 16, italic: true, color: C.muted });

  text(slide, "inertia-heading", "TEST-COMPATIBLE INERTIA ESTIMATE", { left: 682, top: 151, width: 535, height: 28 }, {
    size: 15, bold: true, color: C.terracotta,
  });
  text(slide, "inertia-equation", "V_eq = sign × (λ / 1000) Σ m_j a_j,abs", {
    left: 680, top: 188, width: 555, height: 48,
  }, { font: FONT_TITLE, size: 25, color: C.ink, valign: "middle" });
  text(slide, "inertia-unit-note", "kN · sign aligned and each run zeroed at its start", {
    left: 682, top: 235, width: 535, height: 24,
  }, { size: 15, color: C.muted });

  // Sensor elevations and tributary masses used in the test's force reconstruction.
  addRule(slide, "sensor-wall-line", 715, 285, 715, 444, C.rule, 5);
  shape(slide, "rect", "sensor-wall-foot", { left: 685, top: 444, width: 61, height: 10 }, C.cream);
  const sensorRows = [
    { y: 300, label: "Acc 17 · z = 2.06 m", mass: "m17 = 404.13 kg" },
    { y: 362, label: "Acc 16 · z = 1.26 m", mass: "m16 = 320.66 kg" },
    { y: 424, label: "Acc 15 · z = 0.60 m", mass: "m15 = 408.52 kg" },
  ];
  for (let i = 0; i < sensorRows.length; i++) {
    const r = sensorRows[i];
    shape(slide, "ellipse", "sensor-point-" + i, { left: 707, top: r.y - 8, width: 17, height: 17 }, C.terracotta);
    text(slide, "sensor-label-" + i, r.label, { left: 740, top: r.y - 18, width: 205, height: 25 }, {
      size: 15, bold: true, color: C.ink, valign: "middle",
    });
    text(slide, "sensor-mass-" + i, r.mass, { left: 740, top: r.y + 5, width: 205, height: 23 }, {
      size: 14, color: C.muted, valign: "middle",
    });
  }
  addRule(slide, "sensor-process-divider", 955, 282, 955, 464, C.rule, 1);
  text(slide, "model-accel-heading", "3DEC MODEL", { left: 978, top: 286, width: 235, height: 24 }, {
    size: 14, bold: true, color: C.blue,
  });
  text(slide, "model-accel-copy", "Differentiate absolute sensor velocities; resample at 0.5 ms and low-pass at 50 Hz.", {
    left: 978, top: 314, width: 235, height: 62,
  }, { size: 15, color: C.ink });
  text(slide, "test-accel-heading", "EXPERIMENT", { left: 978, top: 385, width: 235, height: 24 }, {
    size: 14, bold: true, color: C.orange,
  });
  text(slide, "test-accel-copy", "Use table acceleration plus the wall accelerometers to reconstruct absolute acceleration.", {
    left: 978, top: 413, width: 235, height: 55,
  }, { size: 15, color: C.ink });

  addRule(slide, "baseshear-comparison-rule", 78, 515, 1204, 515, C.rule, 1);
  text(slide, "baseshear-likeforlike-label", "LIKE-FOR-LIKE COMPARISON", {
    left: 82, top: 532, width: 280, height: 28,
  }, { size: 15, bold: true, color: C.terracotta, valign: "middle" });
  text(slide, "baseshear-likeforlike-copy", "The same Σm·a definition is used for the experimental data and reconstructed from the model’s sensor histories. Base-shear errors on the next slides use this estimate; cstav stays as context.", {
    left: 362, top: 526, width: 830, height: 53,
  }, { size: 17, color: C.ink, valign: "middle" });
  text(slide, "baseshear-mass-note", "Mass basis: ρ = 1,619 kg/m³ gives the tributary masses above. λ = 2.317348 reproduces the factor in the processed tests and is applied to both sides; its origin is not confirmed.", {
    left: 82, top: 596, width: 1115, height: 48,
  }, { size: 15, color: C.muted, italic: true, valign: "middle" });
  text(slide, "baseshear-xaxis-note", "Loops use relative OOP displacement at z = 2.06 m on the horizontal axis.", {
    left: 82, top: 648, width: 1115, height: 25,
  }, { size: 15, color: C.muted });
  setNotes(slide, "cstav is the direct model support reaction from the base and joist supports; the existing source plots it as a grey context trace. The experimental base shear is derived from accelerometer inertia (Moshfeghi et al. 2024, Eq. 2), not cstav. It uses the three sensor-point absolute accelerations at z=0.60, 1.26, and 2.06 m, weighted by tributary masses 408.52, 320.66, and 404.13 kg. The model side reconstructs those absolute accelerations by differentiating absolute out-of-plane velocities, resampling at 0.5 ms, and applying a 50 Hz low-pass filter. The processed test formula combines table acceleration and wall accelerometer channels to get absolute acceleration. The factor lambda=2.317348 is applied to both sides to reproduce the factor embedded in the processed tests; its origin is unknown and remains a question. Curves are start-zeroed, sign-aligned, and plotted against relative OOP displacement at 2.06 m. The reported model-to-test peak base-shear errors use the experimental-derived sum(m*a) quantity, not cstav.");
}

// 9–12. Base-shear versus displacement plots from the existing project script.
const loopMethods = [
  { n: 8, method: "1_fixed_wall", label: "Method 1", color: C.blue, mean: "9.3%" },
  { n: 9, method: "2_fixed_record", label: "Method 2", color: C.orange, mean: "13.3%" },
  { n: 10, method: "3_adaptive_wall", label: "Method 3", color: C.green, mean: "14.7%" },
  { n: 11, method: "4_adaptive_record", label: "Method 4", color: C.purple, mean: "17.1%" },
];
for (const m of loopMethods) {
  const slide = p.slides.add();
  slide.background.fill = C.white;
  addBodyHeader(slide, m.n, `Base shear versus OOP displacement: ${m.label}`, "Like-for-like loops from baseshear_exp_from_model.py");
  await addPng(slide, `base-shear-loop-${m.label.toLowerCase().replaceAll(" ", "-")}`, path.join(figureDir, m.method, "postproc", "fig_baseshear_exp_loops.png"), {
    left: 45, top: 135, width: 1190, height: 316,
  }, `${m.label} base-shear versus OOP displacement loops for FR76 Runs 22–25, including the source plot's US-1, US-2, and support series.`);
  text(slide, `${m.label}-loop-context`, "US-1 Test 9 is the primary experimental comparator. US-2 Test 12 and the support series remain visible in the source figure.", {
    left: 68, top: 472, width: 1140, height: 33,
  }, { size: 16, color: C.muted, italic: true, valign: "middle" });
  text(slide, `${m.label}-mass-formula-label`, "Base-shear mass basis", { left: 70, top: 530, width: 260, height: 30 }, {
    size: 18, bold: true, color: C.terracotta,
  });
  text(slide, `${m.label}-mass-formula`, "Tributary masses use ρ = 1,619 kg/m³; the ×2.317348 test factor is applied to both model and experiment.", {
    left: 70, top: 561, width: 820, height: 52,
  }, { size: 17, color: C.ink });
  addRule(slide, `${m.label}-mean-error-rule`, 927, 524, 927, 614, C.rule, 1);
  text(slide, `${m.label}-mean-error-caption`, "US-1 mean base-shear peak error", { left: 956, top: 532, width: 270, height: 24 }, {
    size: 15, color: C.muted, align: "center",
  });
  text(slide, `${m.label}-mean-error`, m.mean, { left: 956, top: 558, width: 270, height: 54 }, {
    font: FONT_TITLE, size: 38, bold: true, color: m.color, align: "center",
  });
  setNotes(slide, `${m.label} uses the per-method plot generated from the project copy of Claude outputs/baseshear_exp_from_model.py. The script sums tributary masses at sensors 15–17 using ρ = 1,619 kg/m³ and applies MASS_SCALE = 2.317348 to both the model and experimental sides. The source figure retains both US-1 Test 9 and US-2 Test 12, plus the support reaction traces. The displayed mean base-shear peak error is specific to US-1 Test 9, Runs 22–24; a separate slide scores US-2 over Runs 22–25.`);
}

// 12. US-2 peak-response benchmark, kept separate from the US-1 scorecard.
{
  const slide = p.slides.add();
  slide.background.fill = C.white;
  addBodyHeader(slide, 12, "US-2 is a separate benchmark", "The Test 12 displacement path stays low through Run 24, then rises sharply at Run 25");
  const metricRows = parseCsvRows(await fs.readFile(path.join(figureDir, "postproc", "route2_four_method_metrics.csv"), "utf8"));
  const methods = [
    { key: "1_fixed_wall", label: "Method 1", color: C.blue },
    { key: "2_fixed_record", label: "Method 2", color: C.orange },
    { key: "3_adaptive_wall", label: "Method 3", color: C.green },
    { key: "4_adaptive_record", label: "Method 4", color: C.purple },
  ];
  const testRows = metricRows.filter((r) => r.method === methods[0].key).sort((a, b) => Number(a.run) - Number(b.run));
  const modelSeries = (field) => methods.map((m) => {
    const subset = metricRows.filter((r) => r.method === m.key).sort((a, b) => Number(a.run) - Number(b.run));
    return {
      name: m.label,
      xValues: subset.map((r) => Number(r.run)),
      values: subset.map((r) => chartPrecision(r[field])),
      line: { style: "solid", fill: m.color, width: 2.0 },
      marker: { symbol: "circle", size: 6 },
    };
  });
  const expSeries = (field) => [
    {
      name: "US-1 Test 9",
      xValues: testRows.filter((r) => r[field] !== "").map((r) => Number(r.run)),
      values: testRows.filter((r) => r[field] !== "").map((r) => chartPrecision(r[field])),
      line: { style: "dash", fill: C.ink, width: 1.8 },
      marker: { symbol: "diamond", size: 8 },
    },
    {
      name: "US-2 Test 12",
      xValues: testRows.map((r) => Number(r.run)),
      values: testRows.map((r) => chartPrecision(r[field.replace("US1", "US2")])),
      line: { style: "dot", fill: "#8A5A44", width: 2.0 },
      marker: { symbol: "triangle", size: 8 },
    },
  ];
  const baseChartStyle = {
    scatterOptions: { style: "line" },
    hasLegend: true,
    legend: { position: "bottom", overlay: false, textStyle: { typeface: FONT_BODY, fontSize: 10, fill: C.ink } },
    chartFill: C.white,
    chartLine: { style: "solid", fill: "none", width: 0 },
    plotAreaFill: C.white,
    plotAreaLine: { style: "solid", fill: "none", width: 0 },
    xAxis: {
      title: { text: "FR76 sequence run", textStyle: { typeface: FONT_BODY, fontSize: 14, fill: C.ink } },
      min: 22, max: 25, majorUnit: 1, numberFormatCode: "0",
      textStyle: { typeface: FONT_BODY, fontSize: 12, fill: C.muted },
      line: { style: "solid", fill: "#918B85", width: 1 },
      majorGridlines: null,
    },
    yAxis: {
      min: 0, max: 45, majorUnit: 10, numberFormatCode: "0",
      textStyle: { typeface: FONT_BODY, fontSize: 12, fill: C.muted },
      line: { style: "solid", fill: "#918B85", width: 1 },
      majorGridlines: { style: "solid", fill: "#DDD9D4", width: 1 },
    },
  };
  const oopChart = slide.charts.add("scatter", {
    position: { left: 56, top: 143, width: 565, height: 378 },
    series: [...modelSeries("peak_relative_OOP_mm_model"), ...expSeries("peak_relative_OOP_mm_US1")],
    ...baseChartStyle,
    yAxis: { ...baseChartStyle.yAxis, title: { text: "Peak relative OOP displacement (mm)", textStyle: { typeface: FONT_BODY, fontSize: 14, fill: C.ink } } },
  });
  applyPresentationChartFont(oopChart, { fontFamily: FONT_BODY });
  const shearChart = slide.charts.add("scatter", {
    position: { left: 652, top: 143, width: 572, height: 378 },
    series: [...modelSeries("peak_base_shear_kN_model"), ...expSeries("peak_base_shear_kN_US1")],
    ...baseChartStyle,
    yAxis: { ...baseChartStyle.yAxis, title: { text: "Peak experiment-equivalent base shear (kN)", textStyle: { typeface: FONT_BODY, fontSize: 14, fill: C.ink } } },
  });
  applyPresentationChartFont(shearChart, { fontFamily: FONT_BODY });

  addRule(slide, "us2-benchmark-lower-rule", 72, 532, 1205, 532, C.rule, 1);
  const takeaways = [
    { x: 76, w: 330, label: "RUN 24", value: "3.9×", detail: "US-1: 29.48 mm  ·  US-2: 7.58 mm", color: C.terracotta },
    { x: 463, w: 332, label: "US-2 RUN 25", value: "31.14 mm", detail: "Large excursion; US-1 has no Run 25 test.", color: C.blue },
    { x: 847, w: 365, label: "US-2 MEAN ABSOLUTE ERROR", value: "OOP: M3 32.9%  ·  Shear: M1 7.7%", detail: "M3 tracks through Run 24, then misses Run 25.", color: C.green },
  ];
  for (const [i, d] of takeaways.entries()) {
    if (i > 0) addRule(slide, `us2-benchmark-divider-${i}`, d.x - 28, 550, d.x - 28, 637, C.rule, 1);
    text(slide, `us2-benchmark-label-${i}`, d.label, { left: d.x, top: 550, width: d.w, height: 22 }, {
      size: 13, bold: true, color: C.muted, valign: "middle",
    });
    text(slide, `us2-benchmark-value-${i}`, d.value, { left: d.x, top: 578, width: d.w, height: 34 }, {
      font: FONT_TITLE, size: i === 2 ? 23 : 29, bold: true, color: d.color, valign: "middle",
    });
    text(slide, `us2-benchmark-detail-${i}`, d.detail, { left: d.x, top: 615, width: d.w, height: 34 }, {
      size: 14, color: C.ink, valign: "middle",
    });
  }
  text(slide, "us2-benchmark-scope-note", "The model sequences start from the US-1 Strategy F Run 21 state; Test 12 is an independent specimen check. US-1 and US-2 are not pooled.", {
    left: 76, top: 667, width: 1130, height: 24,
  }, { size: 13, color: C.muted, italic: true, valign: "middle" });
  setNotes(slide, "US-2 Test 12 workbooks were post-processed with the same start-zero peak definition used for US-1: U_avg converted from metres to millimetres, and processed base_shear_kN zeroed at each run start. The four method response peaks are compared with the four Test 12 observations, Runs 22–25; these sequences all start from the common US-1 Strategy F Run 21 model state. US-2 peaks are 2.385, 5.264, 7.581, and 31.142 mm displacement, with 21.191, 27.865, 29.990, and 38.072 kN base shear. US-2 Run 24 displacement is 7.58 mm against US-1's 29.48 mm, while Test 12 then reaches 31.14 mm at Run 25. Mean absolute percentage errors over Runs 22–25: Method 3 is lowest for displacement (32.9%) and Method 1 is lowest for base shear (7.7%). The 32.9% displacement error is still substantial: Method 3 matches US-2 Runs 23–24 closely, but its Run 25 model displacement is only 9.92 mm and it misses the large excursion. Methods 1, 2, and 4 instead produce large model excursions by Run 24. These results show different timing of the transition and specimen-specific sensitivity, not a universal winner.");
}

// 13. Settled readings from both actual experimental tiltmeters; the model's local-chord proxy stays out.
{
  const slide = p.slides.add();
  slide.background.fill = C.white;
  addBodyHeader(slide, 13, "Experimental tiltmeter residuals", "Settled sensor readings for US-1 and US-2; only US-2 has a Run 25 observation");
  const us1TiltRows = parseCsvRows(await fs.readFile(path.join(figureDir, "postproc", "exp_US1_tiltmeter.csv"), "utf8"));
  const us2TiltRows = parseCsvRows(await fs.readFile(path.join(figureDir, "postproc", "exp_US2_tiltmeter.csv"), "utf8"));
  const chart = slide.charts.add("scatter", {
    position: { left: 64, top: 150, width: 750, height: 440 },
    series: [{
      name: "US-1 Test 9",
      xValues: us1TiltRows.filter((r) => Number(r.run) >= 21 && Number(r.run) <= 24).map((r) => Number(r.run)),
      values: us1TiltRows.filter((r) => Number(r.run) >= 21 && Number(r.run) <= 24).map((r) => chartPrecision(r.resid_tilt_y_deg)),
      line: { style: "solid", fill: C.terracotta, width: 3 },
      marker: { symbol: "diamond", size: 8 },
    }, {
      name: "US-2 Test 12",
      xValues: us2TiltRows.filter((r) => Number(r.run) >= 21 && Number(r.run) <= 25).map((r) => Number(r.run)),
      values: us2TiltRows.filter((r) => Number(r.run) >= 21 && Number(r.run) <= 25).map((r) => chartPrecision(r.resid_tilt_y_deg)),
      line: { style: "solid", fill: C.blue, width: 3 },
      marker: { symbol: "triangle", size: 8 },
    }],
    scatterOptions: { style: "line" },
    hasLegend: true,
    legend: { position: "bottom", overlay: false, textStyle: { typeface: FONT_BODY, fontSize: 13, fill: C.ink } },
    chartFill: C.white,
    chartLine: { style: "solid", fill: "none", width: 0 },
    plotAreaFill: C.white,
    plotAreaLine: { style: "solid", fill: "none", width: 0 },
    xAxis: {
      title: { text: "FR76 sequence run", textStyle: { typeface: FONT_BODY, fontSize: 17, fill: C.ink } },
      min: 21, max: 25, majorUnit: 1, numberFormatCode: "0",
      textStyle: { typeface: FONT_BODY, fontSize: 16, fill: C.muted },
      line: { style: "solid", fill: "#918B85", width: 1 },
      majorGridlines: null,
    },
    yAxis: {
      title: { text: "Settled residual tilt Y (degrees)", textStyle: { typeface: FONT_BODY, fontSize: 17, fill: C.ink } },
      min: 0, max: 7, majorUnit: 1,
      textStyle: { typeface: FONT_BODY, fontSize: 15, fill: C.muted },
      line: { style: "solid", fill: "#918B85", width: 1 },
      majorGridlines: { style: "solid", fill: "#DDD9D4", width: 1 },
    },
  });
  applyPresentationChartFont(chart, { fontFamily: FONT_BODY });
  addRule(slide, "tilt-note-rule", 850, 168, 850, 616, C.rule, 1);
  text(slide, "tilt-run24-value", "US-1  +3.266°", { left: 886, top: 183, width: 310, height: 55 }, {
    font: FONT_TITLE, size: 33, bold: true, color: C.terracotta,
  });
  text(slide, "tilt-run24-caption", "Run 24 change from its Run 21 reading", { left: 889, top: 242, width: 300, height: 44 }, {
    size: 17, color: C.ink,
  });
  text(slide, "tilt-us2-value", "US-2  +1.065°", { left: 886, top: 315, width: 310, height: 55 }, {
    font: FONT_TITLE, size: 33, bold: true, color: C.blue,
  });
  text(slide, "tilt-us2-caption", "Run 25 change; settled reading 2.288°", { left: 889, top: 374, width: 300, height: 44 }, {
    size: 17, color: C.ink,
  });
  text(slide, "tilt-clip-note", "Within-shaking peaks clip near ±25.75°; this plot uses settled residual readings only.", {
    left: 888, top: 447, width: 302, height: 65,
  }, { size: 17, color: C.muted });
  text(slide, "tilt-model-note", "These are instrument readings. The model's local-chord angle is excluded from this comparison.", {
    left: 888, top: 538, width: 302, height: 64,
  }, { size: 17, bold: true, color: C.ink });
  setNotes(slide, "US-1 and US-2 tiltmeter residuals are recomputed with the same rule from EXP_DATA/US1_Tilt_values.csv and EXP_DATA/US2_Tilt_values.csv: use accelerometer swing to identify shaking events, then take the settled median after each event relative to the median of the first five samples. US-1 Runs 21–24 = 2.768, 2.533, 2.417, and 6.034 degrees; Run 24 is +3.266 degrees relative to Run 21. US-2 Runs 21–25 = 1.223, 1.0225, 0.385, 0.944, and 2.288 degrees; Run 25 is +1.065 degrees relative to Run 21. The within-shaking inclinometer peaks clip near ±25.75 degrees, so settled readings are shown. This slide compares the two experimental instruments only; it does not plot the model's local-chord tilt proxy.");
}

// 14. Figure-led close: specimen and response metric change the ranking.
{
  const slide = p.slides.add();
  slide.background.fill = C.white;
  addBodyHeader(slide, 14, "What these runs support", "The method ranking depends on response metric and specimen");
  const cols = [
    { x: 94, value: "17.1%", color: C.orange, label: "Lowest US-1 mean OOP\ndisplacement error", method: "Method 2" },
    { x: 477, value: "9.3%", color: C.blue, label: "Lowest US-1 mean base-shear\npeak error", method: "Method 1" },
    { x: 860, value: "16 / 16", color: C.terracotta, label: "Spectral fits outside\n2% tolerance", method: "All methods" },
  ];
  for (const [i, c] of cols.entries()) {
    text(slide, `summary-value-${i}`, c.value, { left: c.x, top: 176, width: 300, height: 76 }, {
      font: FONT_TITLE, size: 54, bold: true, color: c.color, align: "center", valign: "middle",
    });
    text(slide, `summary-method-${i}`, c.method, { left: c.x, top: 264, width: 300, height: 31 }, {
      size: 20, bold: true, color: C.ink, align: "center", valign: "middle",
    });
    text(slide, `summary-label-${i}`, c.label, { left: c.x, top: 303, width: 300, height: 67 }, {
      size: 18, color: C.muted, align: "center", valign: "middle",
    });
    if (i < 2) addRule(slide, `summary-divider-${i}`, c.x + 337, 183, c.x + 337, 370, C.rule, 1);
  }
  addRule(slide, "summary-lower-rule", 87, 414, 1192, 414, C.rule, 1);
  text(slide, "summary-reading", "US-2 adds a different result: Method 3 has the lowest mean displacement error, while Method 1 leads on base shear.", {
    left: 94, top: 438, width: 1090, height: 62,
  }, { font: FONT_TITLE, size: 25, bold: true, color: C.ink, valign: "middle" });
  text(slide, "summary-next-step", "Method 3 still misses the sharp US-2 Run 25 displacement jump; refine the fit and investigate this transition.", {
    left: 94, top: 541, width: 1088, height: 50,
  }, { size: 21, color: C.terracotta, bold: true, valign: "middle" });
  text(slide, "summary-scope-note", "US-1 is scored over Runs 22–24; US-2 over Runs 22–25. Tests stay separate; models start from US-1 Strategy F Run 21.", {
    left: 94, top: 633, width: 1090, height: 34,
  }, { size: 16, italic: true, color: C.muted, valign: "middle" });
  setNotes(slide, "The headline US-1 metrics compare Runs 22–24: Method 2 has the lowest mean absolute OOP error at 17.1%, and Method 1 the lowest mean absolute base-shear peak error at 9.3%. The separate US-2 score over Runs 22–25 gives Method 3 the lowest OOP MAPE at 32.9% and Method 1 the lowest base-shear MAPE at 7.7%. Method 3 predicts only 9.92 mm at US-2 Run 25 versus the measured 31.14 mm. This indicates no universal winner and a specimen-sensitive transition. All 16 spectral fits missed the 2% convergence tolerance. The next step is to improve the joint two-wave fit and investigate the Run 25 transition.");
}

if (p.slides.items.length !== 15) throw new Error("Unexpected slide count: " + p.slides.items.length);

// Export a draft, slide previews, montage, and layout data for visual QA.
const draftPath = path.join(buildDir, "candidate.pptx");
await (await PresentationFile.exportPptx(p)).save(draftPath);
for (let i = 0; i < p.slides.items.length; i++) {
  const slide = p.slides.items[i];
  const preview = await p.export({ slide, format: "png", scale: 1.4 });
  await fs.writeFile(path.join(buildDir, `slide-${String(i + 1).padStart(2, "0")}.png`), new Uint8Array(await preview.arrayBuffer()));
  const layout = await slide.export({ format: "layout" });
  await fs.writeFile(path.join(buildDir, `slide-${String(i + 1).padStart(2, "0")}.layout.json`), await layout.text());
}
const montage = await p.export({ format: "webp", montage: true, scale: 1 });
await fs.writeFile(path.join(buildDir, "montage.webp"), new Uint8Array(await montage.arrayBuffer()));
const inspection = await p.inspect({ kind: "slide,textbox,shape,image,chart,notes", maxChars: 30000 });
await fs.writeFile(path.join(buildDir, "deck-inspect.ndjson"), inspection.ndjson);

// Give editable native charts embedded workbooks before the finalizer validates them.
const materializedCandidatePath = path.join(buildDir, "candidate-with-chart-workbooks.pptx");
const chartSnapshotReceiptPath = path.join(buildDir, "chart-data-snapshot.json");
const chartMaterializerPath = path.join(skillDir, "container_tools/materialize_literal_chart_workbooks.py");
const snapshotReport = execFileSync(pythonExecutable, [
  "-B",
  chartMaterializerPath,
  draftPath,
  materializedCandidatePath,
  "--workspace",
  workspaceDir,
  "--receipt",
  chartSnapshotReceiptPath,
], { encoding: "utf8", maxBuffer: 2 * 1024 * 1024, windowsHide: true });
console.log(snapshotReport.trim());

// Finalize to a new file in 00.Meeting and retain the private validation receipt.
const validationPath = path.join(buildDir, "Progress_meeting_Route2_four_methods_2026-10-05-US2-postprocess-v3.validation.json");
const result = await finalizePresentation({
  workspaceDir,
  candidatePath: materializedCandidatePath,
  finalPath,
  pythonExecutable,
  integrityValidatorPath: path.join(skillDir, "container_tools/inspect_presentation_package_integrity.py"),
  layoutValidatorPath: path.join(skillDir, "container_tools/inspect_presentation_layout_geometry.py"),
  layoutArgs: ["--expected-slide-size-emu", "12192000,6858000", "--validate-heading-fit"],
  explicitTotalSlideCount: 15,
  requiredNativeChartOwnerSlides: [3, 5, 13, 14],
  materializeLiteralChartWorkbooks: false,
  fontPolicy: {
    basis: "reference",
    families: [FONT_TITLE, FONT_BODY],
    referencePath: sourceDeck,
    referenceSha256: "cc000df3cafaa143c4e72e09b9b4f16de2ddec1d95fdb0bf1c723252f70fb715",
  },
  verifyArtifactToolImport: true,
  receiptPath: validationPath,
});
console.log(JSON.stringify({ draftPath, finalPath, validationPath, slideCount: p.slides.items.length, finalizer: result }, null, 2));
