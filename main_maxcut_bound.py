#  Main script that calculates upper bound of a given GSet toroidal Max-Cut problem
import os
import numpy as np
import time
from datetime import datetime
from toroid_maxcut_lib import get_toroid_max_cut_bound
from toroid_faces_lib import find_edges_and_faces_of_toroid
from toroid_lp_lib import get_toroid_max_cut_bound_lp
from params_manager import gsetnum, alg_num
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
        summary = []
        if alg_num & 3:  # Bound of Algorithms #1, #2 (shortest paths between odd faces)
            total_t0 = time.time()
            MCLimit, t101 = get_toroid_max_cut_bound(H, Habs, Hsign, D, GNum)
            total_t = time.time() - total_t0
            formatted_datetime = datetime.now().strftime("%Y %B %d, %H:%M:%S")
            with open("log.txt", "a") as flog:
                flog.write(f"{formatted_datetime} GSet={GNum} MaxCutBound={MCLimit} time={t101}sec\n")
            summary.append(f'Bound={MCLimit}, time={total_t:.1f} sec.')

        if alg_num & 4:  # Liers-Pardella bound. Total times include finding the faces in both cases
            lp_t0 = time.time()
            edges_array, faces_array, _ = find_edges_and_faces_of_toroid(Habs, Hsign, D)
            lp_match_t0 = time.time()
            lp_bound, lp_is_cut = get_toroid_max_cut_bound_lp(H, edges_array, faces_array)
            lp_match_t = int(round(time.time() - lp_match_t0))
            lp_t = time.time() - lp_t0
            print(f'answer: LP Max-Cut bound = {lp_bound}{" (it is a cut, so it is the Max-Cut)" if lp_is_cut else ""}, time = {lp_match_t} sec.')
            print('=====================')
            if not os.path.exists('./bounds'):
                os.makedirs('./bounds')
            np.savetxt(f'./bounds/gset{GNum}_max_cut_bound_lp.txt', np.atleast_1d(lp_bound))
            formatted_datetime = datetime.now().strftime("%Y %B %d, %H:%M:%S")
            with open("log.txt", "a") as flog:
                flog.write(f"{formatted_datetime} GSet={GNum} Alg=LP MaxCutBound={lp_bound} time={lp_match_t}sec\n")
            summary.append(f'LP: Bound={lp_bound}, time={lp_t:.1f} sec.')

        print(f'GSet #{GNum} bounds | ' + ' | '.join(summary))
