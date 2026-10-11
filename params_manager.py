from params_const import myinf
from params_config import gsetnum, use_multi_proc, ncycles, MinPathThreshhold, alg_num
from params_config import lp_max_iters, lp_max_stall
import numpy as np

gset_toroidal = [11, 12, 13, 32, 33, 34, 48, 49, 50, 57, 62, 65, 66, 67, 72, 77, 81]
for gsetnum_current in gsetnum:
    assert gsetnum_current in gset_toroidal
gsetnum = np.array(gsetnum)
assert 1 <= alg_num <= 7
