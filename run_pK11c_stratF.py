# -*- coding: ascii -*-
# run_pK11c_stratF.py -- pocket variant pK11c (pK11 with cohesion-residual 0.1 MPa): JOB "stratF".
# One 3DEC session per job; they can run at the same time (see pocket_jobs.py).
# 3DEC 9, Hanze_SAFEGO folder:  program call 'run_pK11c_stratF.py'
POCKET = "pK11c"
JOB = "stratF"
DAMP_PCT = 1.5
exec(compile(open("pocket_jobs.py").read(), "pocket_jobs.py", "exec"))
