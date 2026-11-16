# -*- coding: ascii -*-
# run_deformed_shape_unclamp.py -- deformed shape and joint damage of the 1.5 % Maxwell Strategy F run with the
# UNCLAMPED pocket (pk2 f_t 0.05, soft tensionless joist top faces). Reference = the unclamped build save, so the
# pocket contacts are compared with their own build values. 16 Nov 2026.
# 3DEC 9, Hanze_SAFEGO folder, no other model open:  program call 'run_deformed_shape_unclamp.py'
DFM_RUN_DIR = "stratF_full_LS_MAXWELL15_zup_pk2_ft0p05_unclamp_tr"
DFM_REF_SAV = "Part_I_MASON_v8_LargeStrain_springforce_zup_pk2_ft0p05_unclamp"
exec(compile(open("deformed_shape_stratF_zup.py").read(), "deformed_shape_stratF_zup.py", "exec"))
