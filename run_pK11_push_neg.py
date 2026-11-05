# -*- coding: ascii -*-
# run_pK11_push_neg.py -- pocket variant pK11 (ks 1.1e8 Pa/m, cohesion 0.1 MPa, friction 31 deg, tension 0): JOB "push-".
# One 3DEC session per job; the three can run at the same time (see pocket_jobs.py).
# 3DEC 9, Hanze_SAFEGO folder:  program call 'run_pK11_push_neg.py'
POCKET = "pK11"
JOB = "push-"
DAMP_PCT = 1.5
exec(compile(open("pocket_jobs.py").read(), "pocket_jobs.py", "exec"))
