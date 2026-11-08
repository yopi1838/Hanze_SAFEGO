# -*- coding: ascii -*-
# run_gluedDefLS_push_pos.py -- DEFORMABLE blocks, LARGE strain, GLUED pockets (joint set, f_t 0.05 MPa everywhere): JOB "push+".
# 3DEC 9, Hanze_SAFEGO folder:  program call 'run_gluedDefLS_push_pos.py'   (+ first: it builds the save)
JOB = "push+"
STRAIN = "LS"
POCKET = "glued"
exec(compile(open("pocket_def_jobs.py").read(), "pocket_def_jobs.py", "exec"))
