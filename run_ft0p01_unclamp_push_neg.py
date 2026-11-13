# -*- coding: ascii -*-
# run_ft0p01_unclamp_push_neg.py -- pk2 build with joint tension f_t = 0.01 MPa (Ihsan, 9 Nov 2026), rigid blocks, large
# strain, joist pockets UNCLAMPED: JOB "push-". The first session builds the f_t 0.01 save and its unclamped variant
# (ANALYSIS_PART_I_MASON_zup_pk2_ft0p01.dat, then ..._ft0p01_unclamp.dat); the other waits (see pocket_jobs.py).
# Output: pushover_results_v5zup_pk2_ft0p01_unclamp_TopMason/pushover_neg.csv
# 3DEC 9, Hanze_SAFEGO folder:  program call 'run_ft0p01_unclamp_push_neg.py'
BASE = "ft0p01"
POCKET = "unclamp"
JOB = "push-"
exec(compile(open("pocket_jobs.py").read(), "pocket_jobs.py", "exec"))
