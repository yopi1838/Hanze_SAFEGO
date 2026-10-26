"""Transcription of Onur Arslan's SafeGo_Mechanic_Model_v1.xlsx (sheet SafeGo_Model_v1), cell by cell.
Units as in the sheet: lengths mm, line loads N/mm (per mm of wall width), forces N/mm -> x b(mm)/1000 = kN."""
import math

def graziotti(t, h, W, O, fw, ab=1, at=1, e=0.0):
    """Graziotti et al. (2016) / Tomassetti et al. (2018) one-way strip, sheet rows 24-47.
    t, h in mm; W (self-weight) and O (overburden) in N/mm; fw in MPa. Returns per-unit-width quantities."""
    Mb = (fw + (W + O) / t) * t ** 2 / 6                         # B24
    Mt = (fw + O / t) * t ** 2 / 6                               # B25
    zeta = fw + 1 / t * (W / 2 + O + 3 * (ab * Mb + at * Mt) / t)   # B28
    dM = ab * Mb - at * Mt
    Fcr = (zeta + math.sqrt(zeta ** 2 - W ** 2 / (4 * t ** 2) - dM / t ** 3 - dM ** 2 / t ** 4)) / (3 * h / 2 / t ** 2)   # B29
    h1 = (0.5 + 1 / h / Fcr * (W * t / 6 + dM)) * h              # B30, B31
    F0 = 2 / h1 * (W + O) * t + O / (h - h1) * (t + 2 * e)       # B34 / B36
    num = 2 / h1 * (W + O) * t + O * (t + 2 * e) / (h - h1)
    den = 2 / h1 * (W + O) + 2 * O / (h - h1)
    uins = num / den                                             # B46 / B47
    return dict(Mb=Mb, Mt=Mt, zeta=zeta, Fcr=Fcr, h1=h1, h2=h - h1, F0=F0, uins=uins)

def cantilever_a0(Wb, O, t, hb, FQ=0.0, edge=False):
    """Casapulla et al. (2017) rigid block rocking on its base, sheet P3/P4 (upper) and P6/P7 (lower, with joist
    friction FQ). Load multiplier a0; the block's capacity is a0 x Wb. Wb, O, FQ in kN; t, hb in m."""
    return (Wb * t / 2 + O * (t if edge else t / 2) + FQ * hb) / (O * hb + Wb * hb / 2)
