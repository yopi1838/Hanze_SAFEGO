# -*- coding: ascii -*-
# run_unclamp_push_pos.py -- pk2 f_t 0.05 rigid build, large strain, joist pockets UNCLAMPED (top-face contacts soft): JOB "push+".
# One 3DEC session per direction; both can run at the same time (see pocket_jobs.py: the first builds the save).
# 3DEC 9, Hanze_SAFEGO folder:  program call 'run_unclamp_push_pos.py'
POCKET = "unclamp"
JOB = "push+"
exec(compile(open("pocket_jobs.py").read(), "pocket_jobs.py", "exec"))
