# Main script that gives an approximated solution for toroidal Max-Cut problem
import random
import os
import atexit
from time import sleep

import numpy as np
import matplotlib.pyplot as plt
import time
from datetime import datetime
from toroid_maxcut_lib import get_toroid_max_cut_topological
from file_lib import write_bin_to_hex_file, load_gset_problem
from toroid_faces_lib import find_edges_and_faces_of_toroid
from toroid_lp_lib import get_toroid_max_cut_lp
from params_manager import gsetnum, alg_num, lp_max_iters, lp_max_stall

t0 = time.time()
test_date = datetime.now().date()
test_time = datetime.now().time()
plot_en = 0
analyse_problems_en = 1
solve_enable = 1
# Format as "Day, Month DD, YYYY HH:MM:SS"
formatted_datetime = datetime.now().strftime("%Y %B %d, %H:%M:%S")
# with open("log.txt", "a") as flog:
#     flog.write(f"{formatted_datetime}\n")


def print_summary_table(rows):
    # Summary of all the problems that were finished (also if the run was stopped in the middle)
    if len(rows) == 0:
        return
    cols = ['GSet', '#Nodes', '#Edges']
    if alg_num & 3:
        cols += [alg_name + ' MC', alg_name + ' t[s]']
    if alg_num & 4:
        cols += ['LP MC', 'LP bound', 'LP iters', 'LP t[s]']
    lines = ['  '.join(f'{c:>10}' for c in cols)]
    lines.append('-' * len(lines[0]))
    for row in rows:
        lines.append('  '.join(f'{str(row.get(c, "")):>10}' for c in cols))
    if alg_num & 3 and alg_num & 4:
        both = [r for r in rows if alg_name + ' MC' in r and 'LP MC' in r]
        lp_better = sum(r['LP MC'] > r[alg_name + ' MC'] for r in both)
        lp_worse = sum(r['LP MC'] < r[alg_name + ' MC'] for r in both)
        lines.append(f'LP vs {alg_name}: better in {lp_better}, equal in {len(both) - lp_better - lp_worse}, worse in {lp_worse} of {len(both)} problems. '
                     f'Total time: {alg_name} {sum(r[alg_name + " t[s]"] for r in both):.1f} sec., LP {sum(r["LP t[s]"] for r in both):.1f} sec.')
    if alg_num & 4:
        lines.append("'*' - LP cut is equal to the LP bound (optimal)")
    print('\n'.join(['', 'Summary:'] + lines))
    with open("log.txt", "a") as flog:
        flog.write('\n'.join(lines) + '\n')


alg_name = 'Alg' + ['', '1', '2', '1+2'][alg_num & 3]

if __name__ == '__main__':
    summary_rows = []
    atexit.register(print_summary_table, summary_rows)
    for GNum in gsetnum:
        H, edges = load_gset_problem(GNum)
        Habs = np.abs(H).astype(bool)
        Hsign = (H<0).astype(bool)
        mat_sz = H.shape[0]
        disp_sz = 800
        if plot_en:
            # plt.imshow(H)
            plt.imshow(abs(H[:disp_sz,:disp_sz]),cmap = 'gray')
            plt.title('GSet #'+str(GNum)+', |V|='+str(mat_sz))
            plt.colorbar()
            plt.show()
            print(gsetnum)
            time.sleep(1)
        # planar_matrix, non_planar_matrix =  plot_graph_from_adjacency_matrix(H)

        if analyse_problems_en:
            cycle_size = []
            for cycle_start in range(0,mat_sz):
                ii_connected = np.argwhere(H[cycle_start, :] != 0)
                ii_connected = ii_connected[ii_connected > cycle_start]
                if len(ii_connected) > 1:
                    if ii_connected[1] < mat_sz/2:
                        cycle_size.append(int(ii_connected[1]) - cycle_start)
            D = int(np.round(np.mean(cycle_size)))
            print('GSet #'+str(GNum)+', #Nodes='+str(mat_sz)+', D='+str(D))

        if solve_enable:
                total_t0 = time.time()
                # cycles_start_step = D - 1
                # cycles_vec = np.arange(0, mat_sz, cycles_start_step)
                cycles_vec = np.arange(mat_sz)
                random.shuffle(cycles_vec)
                force_rerun_arg = 1
                solve_en = 1
                # alg_num = 2
                NEdges = int(np.sum(abs(np.sign(H)))/2)
                summary = []
                summary_row = {'GSet': int(GNum), '#Nodes': mat_sz, '#Edges': NEdges}
                if solve_en and alg_num & 3:  # Algorithms #1, #2
                    mc_val_best, num_ov, num_act_short_paths, mc_groups_best, mc_masks_best = get_toroid_max_cut_topological(H, Habs, Hsign, D, 0, GNum, cycles_vec, force_rerun_arg, alg_num & 3)
                    print('max cut:', mc_val_best)
                    mc_masks_best = np.array(mc_masks_best)
                    print('masks array shape:', mc_masks_best.shape)
                    mc_groups_best = np.array(mc_groups_best)
                    print('groups array shape:', mc_groups_best.shape)
                    if len(mc_masks_best.shape) > 1:
                        fname = 'sols/gset_{:d}_maxcut_masks.txt'.format(GNum)
                        write_bin_to_hex_file(mc_masks_best, fname, "w")
                    if len(mc_groups_best) > 0:
                        fname = 'sols/gset_{:d}_maxcut_groups_hex.txt'.format(GNum)
                        write_bin_to_hex_file(mc_groups_best, fname, "w")

                    vbc_time = time.time()-total_t0
                    delta_t = 'time='+str(int(vbc_time))+' sec.'
                    # res_str = 'GSet #'+str(GNum)+', #Nodes='+str(mat_sz)+', #Edges='+str(NEdges)+', D='+str(D)+'+'+', NActShortPaths='+str(num_act_short_paths)+', MC='+str(int(mc_val_best))+', '+delta_t
                    res_str = f'GSetNum={GNum}, #Nodes={mat_sz}, #Edges={NEdges}, D={D}, #OddVertices={num_ov}, #ActShortPaths={num_act_short_paths}, MC={int(mc_val_best)}, {delta_t}'
                    print(res_str)
                    formatted_datetime = datetime.now().strftime("%Y %B %d, %H:%M:%S")
                    with open("log.txt", "a") as flog:
                        flog.write(f"{formatted_datetime}: "+res_str+'\n')
                    summary.append(f'{alg_name}: MC={int(mc_val_best)}, time={vbc_time:.1f} sec.')
                    summary_row[alg_name + ' MC'] = int(mc_val_best)
                    summary_row[alg_name + ' t[s]'] = round(vbc_time, 1)

                if solve_en and alg_num & 4:  # Liers-Pardella. Time includes finding the faces, as for Algorithms #1, #2
                    print(f'>>>>>>>>>>>>>>>>>>>>>>>>#{GNum}: Running LP...>>>')
                    lp_t0 = time.time()
                    edges_array, faces_array, _ = find_edges_and_faces_of_toroid(Habs, Hsign, D)
                    lp_val, lp_bound, lp_groups, lp_mask, lp_info = get_toroid_max_cut_lp(H, edges_array, faces_array, lp_max_iters, lp_max_stall)
                    if not os.path.exists('./sols'):
                        os.makedirs('./sols')
                    write_bin_to_hex_file(lp_mask[None, :], 'sols/gset_{:d}_maxcut_masks_lp.txt'.format(GNum), "w")
                    write_bin_to_hex_file(lp_groups[None, :], 'sols/gset_{:d}_maxcut_groups_hex_lp.txt'.format(GNum), "w")
                    lp_time = time.time()-lp_t0
                    res_str = f'GSetNum={GNum}, #Nodes={mat_sz}, #Edges={NEdges}, D={D}, Alg=LP, MCBound={lp_bound}, #Iters={lp_info["iters"]}, MC={lp_val}, time={int(lp_time)} sec.'
                    print(res_str)
                    formatted_datetime = datetime.now().strftime("%Y %B %d, %H:%M:%S")
                    with open("log.txt", "a") as flog:
                        flog.write(f"{formatted_datetime}: "+res_str+'\n')
                    summary.append(f'LP: MC={lp_val}{" (optimal)" if lp_info["optimal"] else ""}, MCBound={lp_bound}, time={lp_time:.1f} sec.')
                    summary_row['LP MC'] = lp_val
                    summary_row['LP bound'] = str(lp_bound) + ('*' if lp_info["optimal"] else '')
                    summary_row['LP iters'] = lp_info["iters"]
                    summary_row['LP t[s]'] = round(lp_time, 1)

                print(f'GSet #{GNum} run times | ' + ' | '.join(summary))
                summary_rows.append(summary_row)
                print('Total Time =',int(time.time()-total_t0),'sec.')
        # sleep(15)