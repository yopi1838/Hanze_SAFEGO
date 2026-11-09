# -*- coding: ascii -*-
# run_unclamp_stratF.py -- pk2 f_t 0.05 rigid build, large strain, joist pockets UNCLAMPED (top-face contacts soft):
# JOB "stratF", Strategy F runs 1-25 at 1.5 % Maxwell. Uses the save the unclamped pushovers built
# (Part_I_MASON_v8_LargeStrain_springforce_zup_pk2_ft0p05_unclamp.sav); builds it first if it is missing.
# Output: stratF_full_LS_MAXWELL15_zup_pk2_ft0p05_unclamp_tr/ (resumes from its own checkpoint if stopped).
# 3DEC 9, Hanze_SAFEGO folder:  program call 'run_unclamp_stratF.py'
POCKET = "unclamp"
JOB = "stratF"
DAMP_PCT = 1.5
exec(compile(open("pocket_jobs.py").read(), "pocket_jobs.py", "exec"))
