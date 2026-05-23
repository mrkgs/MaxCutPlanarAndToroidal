# Main script that gives an approximated solution for toroidal Max-Cut problem
import random
import numpy as np
import matplotlib.pyplot as plt
import time
import os
from datetime import datetime
from toroid_maxcut_lib import get_toroid_max_cut_topological, get_toroid_max_cut_direct_painting
from file_lib import write_bin_to_hex_file, load_gset_problem
from params_config import gsetnum

t0 = time.time()
test_date = datetime.now().date()
test_time = datetime.now().time()
assert gsetnum in [11, 12, 13, 32, 33, 34, 48, 49, 50, 57, 62, 65, 66, 67, 72, 77, 81]
nums = np.array([gsetnum])
plot_en = 0
analyse_problems_en = 1
solve_enable = 1
# Format as "Day, Month DD, YYYY HH:MM:SS"
formatted_datetime = datetime.now().strftime("%Y %B %d, %H:%M:%S")
# with open("log.txt", "a") as flog:
#     flog.write(f"{formatted_datetime}\n")

if __name__ == '__main__':
    for GNum in nums:
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
            print(nums+1)
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
                painting_en = 1
                alg_num = 2
                if solve_en:
                    inner_cycles_starts = get_toroid_max_cut_topological(H, Habs, Hsign, D, 0, GNum, cycles_vec, force_rerun_arg, alg_num)
                if painting_en:
                    if not os.path.exists('./sols'):
                        os.makedirs('./sols')
                    max_cut_val_best, max_cut_log, max_cut_masks_best, max_cut_grps_best = get_toroid_max_cut_direct_painting(H, Habs, Hsign, D, 0, GNum, inner_cycles_starts, alg_num)
                    max_cut_masks_best = np.array(max_cut_masks_best)
                    print('masks array shape:', max_cut_masks_best.shape)
                    max_cut_grps_best = np.array(max_cut_grps_best)
                    print('grps array shape:', max_cut_grps_best.shape)
                    if len(max_cut_masks_best.shape) > 1:
                        fname = 'sols/gset_{:d}_maxcut_masks.txt'.format(GNum)
                        write_bin_to_hex_file(max_cut_masks_best, fname, "w")
                    if len(max_cut_grps_best)>0:
                        fname = 'sols/gset_{:d}_maxcut_groups_hex.txt'.format(GNum)
                        write_bin_to_hex_file(max_cut_grps_best, fname, "w")

                    delta_t = 'time: '+str(int(time.time()-t0))+' sec.'
                    NEdges = int(np.sum(abs(np.sign(H)))/2)
                    res_str = 'GSet #'+str(GNum)+', #Nodes='+str(mat_sz)+', #Edges='+str(NEdges)+', D='+str(D)+', MC='+str(max_cut_val_best)+','+delta_t
                    print(res_str)
                    with open("log.txt", "a") as flog:
                        flog.write(res_str+'\n')
                print('Total Time =',int(time.time()-total_t0),'sec.')

