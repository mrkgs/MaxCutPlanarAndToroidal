#  Main script that calculates upper bound of a given GSet toroidal Max-Cut problem
import numpy as np
import time
from datetime import datetime
from toroid_maxcut_lib import get_toroid_max_cut_bound
from params_manager import gsetnum
from file_lib import load_gset_problem

t0 = time.time()
test_date = datetime.now().date()
test_time = datetime.now().time()
formatted_datetime = datetime.now().strftime("%Y %B %d, %H:%M:%S")
print(formatted_datetime)
print('=====================')


if __name__ == '__main__':
    for GNum in gsetnum:
        H, edges = load_gset_problem(GNum)
        Habs = np.abs(H).astype(bool)
        Hsign = (H<0).astype(bool)
        mat_sz = Habs.shape[0]
        cycle_size = []
        for cycle_start in range(0,mat_sz):
            ii_connected = np.argwhere(H[cycle_start, :] != 0)
            ii_connected = ii_connected[ii_connected > cycle_start]
            if len(ii_connected) > 1:
                if ii_connected[1] < mat_sz/2:
                    cycle_size.append(int(ii_connected[1]) - cycle_start)
        D = int(np.round(np.mean(cycle_size)))
        print('GSet #'+str(GNum)+', #Nodes='+str(mat_sz)+', D='+str(D))
        total_t0 = time.time()
        MCLimit, t101 = get_toroid_max_cut_bound(H, Habs, Hsign, D, GNum)
        formatted_datetime = datetime.now().strftime("%Y %B %d, %H:%M:%S")
        with open("log.txt", "a") as flog:
            flog.write(f"{formatted_datetime} GSet={GNum} MaxCutBound={MCLimit} time={t101}sec\n")



