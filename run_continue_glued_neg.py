# -*- coding: ascii -*-
# run_continue_glued_neg.py -- continue the glued-pocket neg pushover from its 80 mm end save to F = 0 (pushover_continue_zup.py).
# 3DEC 9, Hanze_SAFEGO folder:  program call 'run_continue_glued_neg.py'
CASE = "glued"
PUSH_DIR = -1
exec(compile(open("pushover_continue_zup.py").read(), "pushover_continue_zup.py", "exec"))
