# -*- coding: ascii -*-
# run_pushover_zup_pk2_char_neg.py -- the - push of the characterisation-strength pk2 build, on its own (16 Nov 2026).
# Run it and run_pushover_zup_pk2_char_pos.py in two 3DEC sessions at the same time; the first to start builds the save.
# 3DEC 9, Hanze_SAFEGO folder:  program call 'run_pushover_zup_pk2_char_neg.py'
#   -> pushover_results_v5zup_pk2_char_TopMason/pushover_neg.csv
ONLY = -1
exec(compile(open("run_pushover_zup_pk2_char_queue.py").read(), "run_pushover_zup_pk2_char_queue.py", "exec"))
