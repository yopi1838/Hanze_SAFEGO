# -*- coding: ascii -*-
"""
capacity_law.py -- bilinear/secant period law from a capacity curve.

Implements the period law Ihsan proposed (meeting 2026-08-31): the initial
stiffness comes from the linear branch of the force-displacement curve, the
secant stiffness from F(d)/d at the largest displacement of interest, and

    T1    = 2 pi sqrt(m_eff / K_init)
    T_eff = 2 pi sqrt(m_eff / K_sec(d))

The curve can come from
  (a) the MODEL's own pushover (pushover_oop_3dec.py output:
      lam_g, d_ctrl_mm, F_base_kN) -- the transferable, non-circular source,
      required for any PREDICTION claim; or
  (b) the digitised EXPERIMENTAL Fig 13 envelope (drift_pct, F_kN, d_mm) --
      legitimate for VALIDATION against US-1 only, because the largest
      observed displacement is itself the quantity being predicted.
The loader detects which format it was given from the header.

Direction: the OOP response is asymmetric (Table 4: +29.5 kN at 25 mm,
-29.5 kN at -30 mm). DIRECTION picks the branch. "pos"/"neg" use one branch;
"envelope" interpolates F on |d| from whichever branch is weaker at that
|d| (conservative).

Beyond the digitised/pushed range the force is held at the last value and
the note says so -- the secant keeps softening with d, which is the correct
qualitative behaviour, but numbers past the curve end are extrapolation and
are flagged.

Self-tests run on import (SELFTEST = True): a synthetic bilinear curve must
return its own K_init and K_sec exactly, and the units round-trip must give
T = 2 pi sqrt(m/K) for a hand-checked case. A failed self-test raises --
the same discipline as sdof.py, after the Nigam-Jennings episode.
"""
import csv as _csv
import math

SELFTEST = True
_G = 9.81


class CapacityLaw(object):
    def __init__(self, csv_path, m_eff, d_linear_mm=2.0, t_eff_max=1.0,
                 direction="pos", label=""):
        self.path = csv_path
        self.m_eff = float(m_eff)
        self.d_linear_mm = float(d_linear_mm)
        self.t_eff_max = float(t_eff_max)
        self.direction = direction
        self.label = label or csv_path
        self._load(csv_path)
        self._fit_kinit()

    # ------------------------------------------------------------------
    def _load(self, path):
        with open(path, "r") as f:
            rows = list(_csv.reader(f))
        header = [h.strip().lower() for h in rows[0]]
        if "d_ctrl_mm" in header:            # pushover driver output
            i_d = header.index("d_ctrl_mm")
            # v2.1 drivers carry F_tot (base+top reaction, the Eq.-2
            # comparable quantity); older CSVs only F_base. Prefer total.
            i_f = (header.index("f_tot_kn") if "f_tot_kn" in header
                   else header.index("f_base_kn"))
            self.source_kind = "model_pushover"
        elif "d_mm" in header:               # digitised Fig 13 envelope
            i_d = header.index("d_mm")
            i_f = header.index("f_kn")
            self.source_kind = "experimental_envelope"
        else:
            raise RuntimeError("capacity_law: unrecognised header in {}: {}"
                               .format(path, header))
        pts = []
        for r in rows[1:]:
            try:
                pts.append((float(r[i_d]), float(r[i_f])))
            except (ValueError, IndexError):
                continue
        if len(pts) < 4:
            raise RuntimeError("capacity_law: too few points in " + path)
        pos = sorted([(d, f) for d, f in pts if d > 0])
        neg = sorted([(-d, -f) for d, f in pts if d < 0])   # mirrored to +
        if self.direction == "pos":
            branch = pos or neg
        elif self.direction == "neg":
            branch = neg or pos
        else:                                # "envelope": weaker of the two
            branch = self._envelope(pos, neg)
        if not branch:
            raise RuntimeError("capacity_law: no points on requested branch")
        self.d = [p[0] for p in branch]      # mm, ascending, > 0
        self.F = [p[1] for p in branch]      # kN, sign folded to +
        self.d_max = self.d[-1]

    @staticmethod
    def _envelope(pos, neg):
        if not pos:
            return neg
        if not neg:
            return pos
        out = []
        for d, f in pos:
            fn = CapacityLaw._interp_on(neg, d)
            out.append((d, min(f, fn) if fn is not None else f))
        return out

    @staticmethod
    def _interp_on(branch, d):
        if not branch or d < branch[0][0] or d > branch[-1][0]:
            return None
        for i in range(1, len(branch)):
            if branch[i][0] >= d:
                d0, f0 = branch[i - 1]
                d1, f1 = branch[i]
                if d1 == d0:
                    return f0
                return f0 + (f1 - f0) * (d - d0) / (d1 - d0)
        return branch[-1][1]

    # ------------------------------------------------------------------
    def _fit_kinit(self):
        """Least-squares slope through the origin over the linear branch."""
        pts = [(d, f) for d, f in zip(self.d, self.F)
               if 0.0 < d <= self.d_linear_mm]
        if len(pts) < 2:                     # curve starts coarse: first 3 pts
            pts = list(zip(self.d[:3], self.F[:3]))
        num = sum(d * f for d, f in pts)
        den = sum(d * d for d, f in pts)
        if den <= 0:
            raise RuntimeError("capacity_law: cannot fit K_init")
        self.K_init_kNmm = num / den                       # kN/mm
        self.n_linear_pts = len(pts)

    # ------------------------------------------------------------------
    def F_of_d(self, d_mm):
        """kN at d_mm (>0). Flat extrapolation past the curve end."""
        if d_mm <= self.d[0]:
            return self.F[0] * (d_mm / self.d[0]) if self.d[0] > 0 else 0.0
        v = self._interp_on(list(zip(self.d, self.F)), d_mm)
        return v if v is not None else self.F[-1]

    def K_sec_kNmm(self, d_mm):
        if d_mm <= 0:
            return self.K_init_kNmm
        d_eval = min(d_mm, self.d_max)
        return self.F_of_d(d_eval) / d_eval

    @staticmethod
    def _period(m_eff, K_kNmm):
        K_SI = K_kNmm * 1e6                                # kN/mm -> N/m
        return 2.0 * math.pi * math.sqrt(m_eff / K_SI)

    def T1(self):
        return self._period(self.m_eff, self.K_init_kNmm)

    def t_eff(self, d_mm):
        """(T_eff, note) at displacement d_mm. Capped at t_eff_max."""
        note = ""
        if d_mm > self.d_max:
            note = "beyond_curve_end({:.1f}>{:.1f}mm;F_held)".format(
                d_mm, self.d_max)
        K = self.K_sec_kNmm(d_mm)
        T = self._period(self.m_eff, K)
        if T > self.t_eff_max:
            note = (note + "|" if note else "") + "T_capped_{:.2f}s".format(
                self.t_eff_max)
            T = self.t_eff_max
        return T, note

    def report(self):
        lines = [
            "capacity law [{}]  source={}  direction={}".format(
                self.label, self.source_kind, self.direction),
            "  points: {}  d range: {:.2f}-{:.2f} mm  F max: {:.2f} kN".format(
                len(self.d), self.d[0], self.d_max, max(self.F)),
            "  K_init = {:.3f} kN/mm  (fit over d <= {:.1f} mm, {} pts)".format(
                self.K_init_kNmm, self.d_linear_mm, self.n_linear_pts),
            "  T1(m={:.0f} kg) = {:.4f} s".format(self.m_eff, self.T1()),
        ]
        for dq in (5.0, 10.0, 20.0, min(29.5, self.d_max)):
            T, note = self.t_eff(dq)
            lines.append("  d = {:5.1f} mm : K_sec = {:.3f} kN/mm  "
                         "T_eff = {:.4f} s  {}".format(
                             dq, self.K_sec_kNmm(dq), T, note))
        return "\n".join(lines)


# =====================================================================
# SELF-TESTS -- exact bilinear synthetic; fail loudly on import
# =====================================================================
def _selftest():
    import tempfile, os
    # bilinear: K1 = 10 kN/mm to 2 mm (F=20), then K2 = 1 kN/mm
    rows = ["d_mm,F_kN"]
    d = 0.0
    while d <= 30.0001:
        F = 10.0 * d if d <= 2.0 else 20.0 + 1.0 * (d - 2.0)
        rows.append("{:.4f},{:.6f}".format(d, F))
        d += 0.25
    fd, path = tempfile.mkstemp(suffix=".csv")
    with os.fdopen(fd, "w") as f:
        f.write("\n".join(rows))
    try:
        law = CapacityLaw(path, m_eff=1000.0, d_linear_mm=2.0, label="selftest")
        # 1. K_init recovered exactly
        assert abs(law.K_init_kNmm - 10.0) < 1e-9, law.K_init_kNmm
        # 2. K_sec at 20 mm: F = 20 + 18 = 38 kN -> 1.9 kN/mm
        assert abs(law.K_sec_kNmm(20.0) - 1.9) < 1e-9, law.K_sec_kNmm(20.0)
        # 3. period units: m=1000 kg, K=10 kN/mm = 1e7 N/m -> T = 0.06283 s
        assert abs(law.T1() - 2.0 * math.pi * math.sqrt(1000.0 / 1e7)) < 1e-12
        # 4. monotone softening: T_eff grows with d
        T5 = law.t_eff(5.0)[0]
        T25 = law.t_eff(25.0)[0]
        assert T25 > T5 > law.T1()
        # 5. beyond curve end: flat F, note set
        Tb, nb = law.t_eff(100.0)
        assert "beyond_curve_end" in nb
    finally:
        os.remove(path)


if SELFTEST:
    _selftest()