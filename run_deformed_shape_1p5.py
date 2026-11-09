# -*- coding: ascii -*-
# run_deformed_shape_1p5.py -- deformed shape and joint damage of the 1.5 % Maxwell Strategy F run (glued pocket,
# pk2 f_t 0.05). Same build save as the 5 % run (REF_SAV in deformed_shape_stratF_zup.py).
# 3DEC 9, Hanze_SAFEGO folder, no other model open:  program call 'run_deformed_shape_1p5.py'
DFM_RUN_DIR = "stratF_full_LS_MAXWELL15_zup_pk2_ft0p05_tr"
exec(compile(open("deformed_shape_stratF_zup.py").read(), "deformed_shape_stratF_zup.py", "exec"))
