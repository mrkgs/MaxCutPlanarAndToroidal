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
    runtimes = []  # per GSet: (GNum, #nodes, edges+faces time, Alg#1/#2 time, LP time)
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
        # Edges and faces are shared by all bounds and are excluded from their run times
        faces_t0 = time.time()
        edges_array, faces_array, neighbor_faces_array = find_edges_and_faces_of_toroid(Habs, Hsign, D)
        faces_t = time.time() - faces_t0
        total_t, lp_t = None, None
        if alg_num & 3:  # Bound of Algorithms #1, #2 (shortest paths between odd faces)
            for cnt in range(len(edges_array)):
                assert neighbor_faces_array[cnt, 0] != neighbor_faces_array[cnt, 1], 'Currently, working only with graphs consisting of cycles - each edge must have 2 faces'
            if not os.path.exists('./info'):
                os.makedirs('./info')
            np.savetxt('info/edges_'+str(GNum)+'.txt', edges_array)
            np.savetxt('info/faces_' + str(GNum) + '.txt', faces_array)
            np.savetxt('info/neighbors_' + str(GNum) + '.txt', neighbor_faces_array)
            total_t0 = time.time()
            MCLimit, t101 = get_toroid_max_cut_bound(H, Habs, Hsign, edges_array, faces_array, neighbor_faces_array, GNum)
            total_t = time.time() - total_t0
            formatted_datetime = datetime.now().strftime("%Y %B %d, %H:%M:%S")
            with open("log.txt", "a") as flog:
                flog.write(f"{formatted_datetime} GSet={GNum} MaxCutBound={MCLimit} time={t101}sec\n")
            summary.append(f'Bound={MCLimit}, time={total_t:.1f} sec.')

        if alg_num & 4:  # Liers-Pardella bound
            lp_t0 = time.time()
            lp_bound, lp_is_cut = get_toroid_max_cut_bound_lp(H, edges_array, faces_array)
            lp_match_t = int(round(time.time() - lp_t0))
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
        runtimes.append((GNum, mat_sz, faces_t, total_t, lp_t))

    def fmt_t(t):
        return '-' if t is None else f'{t:.1f}'
    print('=====================')
    print('Runtimes summary [sec] (Alg#1/#2 and LP times exclude edges and faces)')
    print(f"{'GSet':>6} {'#Nodes':>7} {'Faces':>8} {'Alg#1/#2':>9} {'LP':>8}")
    for GNum, mat_sz, faces_t, total_t, lp_t in runtimes:
        print(f'{GNum:>6} {mat_sz:>7} {fmt_t(faces_t):>8} {fmt_t(total_t):>9} {fmt_t(lp_t):>8}')
    sum_col = lambda k: None if all(r[k] is None for r in runtimes) else sum(r[k] for r in runtimes if r[k] is not None)
    print(f"{'Total':>6} {'':>7} {fmt_t(sum_col(2)):>8} {fmt_t(sum_col(3)):>9} {fmt_t(sum_col(4)):>8}")
    print(f'Script total time: {time.time() - t0:.1f} sec.')
