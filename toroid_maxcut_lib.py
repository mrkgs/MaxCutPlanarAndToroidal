from planar_maxcut_lib import get_min_odd_paths, get_full_ov_min_cover, get_weighted_signed_max_cut_incremental, paint_nodes_for_maxcut_using_pol
from planar_maxcut_lib import get_painting_list
from toroid_faces_lib import find_edges_and_faces_of_toroid
from toroid_cycles_lib import find_inner_cycle_path_of_toroid
import time, os
import numpy as np
import pickle
import networkx as nx
from multiprocessing import Process
from file_lib import read_from_hex_file_to_bin
from graph_lib import get_clock_wise_dir_edges
# from params_config import use_multi_proc, ncycles
from params_manager import myinf, MinPathThreshhold
from params_manager import use_multi_proc as use_multi_proc_prm
from params_manager import ncycles as number_of_cycles

def get_toroid_max_cut_bound(H, Habs, Hsign, D, gsetnum):
    # MinPathThreshhold = 5
    # MinPathThreshhold = 30000
    rerun_get_min_odd_paths = 1
    force_find_faces = 1
    save_faces = 1
    if not os.path.exists('./info'):
        os.makedirs('./info')
    if not os.path.exists('./bounds'):
        os.makedirs('./bounds')

    if os.path.exists('info/edges_'+str(gsetnum)+'.txt') and not force_find_faces:
        edges_array = np.loadtxt('info/edges_'+str(gsetnum)+'.txt').astype('int64')
        faces_array = np.loadtxt('info/faces_' + str(gsetnum) + '.txt').astype('int64')
        neighbor_faces_array = np.loadtxt('info/neighbors_' + str(gsetnum) + '.txt').astype('int64')
    else:
        # edges_array, faces_array, neighbor_faces_array = find_edges_and_faces_of_toroid(H, D)
        edges_array, faces_array, neighbor_faces_array = find_edges_and_faces_of_toroid(Habs, Hsign, D)
        for cnt in range(len(edges_array)):
            assert neighbor_faces_array[cnt, 0] != neighbor_faces_array[cnt, 1], 'Currently, working only with graphs consisting of cycles - each edge must have 2 faces'
        if save_faces:
            np.savetxt('info/edges_'+str(gsetnum)+'.txt', edges_array)
            np.savetxt('info/faces_' + str(gsetnum) + '.txt', faces_array)
            np.savetxt('info/neighbors_' + str(gsetnum) + '.txt', neighbor_faces_array)

    min_paths_t0 = time.time()
    short_data_fname = f'./info/short_paths_data_{gsetnum}.npz'
    short_base_fname = f'./info/short_paths_base_{gsetnum}.pkl'
    if os.path.exists(short_data_fname) and not rerun_get_min_odd_paths:
        short_dataz = np.load(short_data_fname)
        H_dual, H_dual_cnt, H_edges, odd_vertices, P = short_dataz['arr_0'], short_dataz['arr_1'], short_dataz['arr_2'], short_dataz['arr_3'], short_dataz['arr_4']
        with open(short_base_fname,"rb") as fpkl:
            short_paths_base = pickle.load(fpkl)
    else:
        short_paths_base, H_dual, H_dual_cnt, H_edges, odd_vertices, P = get_min_odd_paths(Hsign, edges_array, faces_array, neighbor_faces_array, MinPathThreshhold)
        np.savez(short_data_fname, H_dual, H_dual_cnt, H_edges, odd_vertices, P)
        with open(short_base_fname,"wb") as fpkl:
            pickle.dump(short_paths_base, fpkl)
    print('finished min_odd_paths(), time:', int(time.time()-min_paths_t0),'sec.')
    print(f'number of odd vertices = {len(odd_vertices)}.')
    print(f'shortest paths threshold = {MinPathThreshhold},')
    print('started searching odd vertices min cover...')
    t100 = time.time()
    full_min_cover, MCLimit  = get_full_ov_min_cover(P.astype(np.uint8), odd_vertices, short_paths_base, H_edges, Habs, Hsign, 0)
    t101 = int(round(time.time() - t100))
    print(f'answer: Max-Cut bound = {MCLimit}, time = {t101} sec.')
    print('=====================')
    # np.savetxt(f'./bounds/gset{gsetnum}_min_cover_time.txt', np.atleast_1d(t101))
    np.savetxt(f'./bounds/gset{gsetnum}_max_cut_bound.txt', np.atleast_1d(MCLimit))
    return MCLimit, t101

def get_toroid_max_cut_topological(H, Habs, Hsign, D, print_en, gsetnum, cycles_starts_to_check, force_rerun, alg_num):
    if gsetnum in [77, 81]:
        use_multi_proc = 0
        print(f'Setting use_multi_proc=0 for gsetnum={gsetnum}')
    else:
        use_multi_proc = use_multi_proc_prm
    if gsetnum in [48, 49, 50]:
        ncycles = 1
    else:
        ncycles = number_of_cycles
    if not os.path.exists('./info'):
        os.makedirs('./info')
    if not os.path.exists('./bounds'):
        os.makedirs('./bounds')
    if not os.path.exists('./tmp'):
        os.makedirs('./tmp')
    if not os.path.exists('./sols'):
        os.makedirs('./sols')
    print_str = 'GSet'+str(gsetnum)
    t0 = time.time()
    N = len(H)
    # MinPathThreshhold = int(0.09*N/D)
    # MinPathThreshhold = 5
    # MinPathThreshhold = 30000
    # MinPathThreshhold = 100
    mc_val_best = -myinf
    force_find_faces = 1
    save_faces = 1
    if os.path.exists('info/edges_'+str(gsetnum)+'.txt') and not force_find_faces:
        edges_array = np.loadtxt('info/edges_'+str(gsetnum)+'.txt').astype('int64')
        faces_array = np.loadtxt('info/faces_' + str(gsetnum) + '.txt').astype('int64')
        neighbor_faces_array = np.loadtxt('info/neighbors_' + str(gsetnum) + '.txt').astype('int64')
    else:
        # edges_array, faces_array, neighbor_faces_array = find_edges_and_faces_of_toroid(H, D)
        edges_array, faces_array, neighbor_faces_array = find_edges_and_faces_of_toroid(Habs, Hsign, D)
        for cnt in range(len(edges_array)):
            assert neighbor_faces_array[cnt, 0] != neighbor_faces_array[cnt, 1], 'Currently, working only with graphs consisting of cycles - each edge must have 2 faces'
        if save_faces:
            np.savetxt('info/edges_'+str(gsetnum)+'.txt', edges_array)
            np.savetxt('info/faces_' + str(gsetnum) + '.txt', faces_array)
            np.savetxt('info/neighbors_' + str(gsetnum) + '.txt', neighbor_faces_array)

    min_paths_t0 = time.time()
    short_data_fname = f'info/short_paths_data_{gsetnum}.npz'
    short_base_fname = f'info/short_paths_base_{gsetnum}.pkl'
    if os.path.exists(short_data_fname):
        short_dataz = np.load(short_data_fname)
        H_dual, H_dual_cnt, H_edges, odd_vertices, P = short_dataz['arr_0'], short_dataz['arr_1'], short_dataz['arr_2'], short_dataz['arr_3'], short_dataz['arr_4']
        with open(short_base_fname,"rb") as fpkl:
            short_paths_base = pickle.load(fpkl)
    else:
        short_paths_base, H_dual, H_dual_cnt, H_edges, odd_vertices, P = get_min_odd_paths(Hsign, edges_array, faces_array, neighbor_faces_array, MinPathThreshhold)
        np.savez(short_data_fname, H_dual, H_dual_cnt, H_edges, odd_vertices, P)
        with open(short_base_fname,"wb") as fpkl:
            pickle.dump(short_paths_base, fpkl)
        # exit(0)
    print('Finished min_odd_paths(), time:', int(time.time()-min_paths_t0),'sec.')
    check_mc_time = 0
    if check_mc_time:
        print(f'#OV={len(odd_vertices)}')
        print('Starting searching ov min cover...')
        t100 = time.time()
        full_min_cover = get_full_ov_min_cover(P.astype(np.uint8), odd_vertices, short_paths_base, H_edges, Habs, Hsign, gsetnum)
        t101 = int(round(time.time() - t100))
        print(t101,'seconds')
        np.savetxt(f'gset{gsetnum}_min_cover_time.txt', np.atleast_1d(t101))
        exit(0)

    inner_cycles = []
    inner_cycles_e = []
    inner_cycles_dual = []
    start_idxs = []
    vbc_pos_odd = []
    force_check_cycles = 1
    cycles_starts_to_check=cycles_starts_to_check[cycles_starts_to_check<N]
    for in_cycle_0 in cycles_starts_to_check:
        cycle_file1 = 'info/' + str(gsetnum) + '_cycle_' + str(in_cycle_0) + '_nodes.txt'
        cycle_file2 = 'info/' + str(gsetnum) + '_cycle_' + str(in_cycle_0) + '_duals.txt'
        cycle_file3 = 'info/' + str(gsetnum) + '_cycle_' + str(in_cycle_0) + '_edges.txt'
        if not os.path.exists(cycle_file3) or force_check_cycles:
            inner_cycle_path, inner_cycle_dual_edges, inner_cycle_edges = find_inner_cycle_path_of_toroid(Habs, edges_array, neighbor_faces_array, in_cycle_0)
            if len(inner_cycle_path)>0 and len(inner_cycle_dual_edges)>0 and len(inner_cycle_edges)>0:
                np.savetxt(cycle_file1, inner_cycle_path)
                np.savetxt(cycle_file2, inner_cycle_dual_edges)
                np.savetxt(cycle_file3, inner_cycle_edges)
        else:
            inner_cycle_path = [int(x) for x in np.loadtxt(cycle_file1)]
            inner_cycle_dual_edges = [[int(x[0]), int(x[1])]  for x in np.loadtxt(cycle_file2)]
            inner_cycle_edges = [[int(x[0]), int(x[1])]  for x in np.loadtxt(cycle_file3)]

        pos_edges_cnt = 0
        if len(inner_cycle_path) == 0:
            valid_inner_cycle = False
        else:
            path_len = len(inner_cycle_path)
            for vcnt in range(path_len):
                if Hsign[inner_cycle_path[vcnt], inner_cycle_path[(vcnt+1)%path_len]] == 0: # if positive
                    pos_edges_cnt+= 1
            # valid_inner_cycle = (pos_edges_cnt % 2 == 0)  # Take only even cycles
            valid_inner_cycle = 1                         # Take all
        evenness_str = ['EVEN', 'ODD ']
        if valid_inner_cycle:
            inner_cycles.append(inner_cycle_path)
            inner_cycles_dual.append(inner_cycle_dual_edges)
            inner_cycles_e.append(inner_cycle_edges)
            print(evenness_str[pos_edges_cnt % 2], inner_cycle_path)
            start_idxs.append(in_cycle_0)
            vbc_pos_odd.append(pos_edges_cnt % 2)
        num_inner_cycles = len(inner_cycles)
        if num_inner_cycles == ncycles:
            break
    cycle_cross_numbers = np.zeros(num_inner_cycles)
    print(f"GSet {gsetnum}, #Nodes={N}, #Edges={len(edges_array)}: found {num_inner_cycles} VBC")
    n_proc = 6
    # use_multi_proc = 0
    num_of_runs = int(np.ceil(num_inner_cycles/n_proc))
    cycles_left_to_run = num_inner_cycles
    cycles_list = np.arange(num_inner_cycles)
    # cycles_list = [0, 1, 2, 3, 4, 5]
    skip_all_runs = 0
    # force_rerun = 0
    # np.save('Puint16.npy',P.astype(np.uint16))
    NOddVert = len(odd_vertices)
    print(f'#OV={NOddVert}')
    NActShortPaths = (P!=0).sum()
    print(f'#ActiveShortPaths = {NActShortPaths} / {NOddVert**2}')
    print('Starting searching ov min cover...')
    t100 = time.time()
    full_min_cover, MCLimit = get_full_ov_min_cover(P.astype(np.uint16), odd_vertices, short_paths_base, H_edges, Habs, Hsign, gsetnum)
    t101 = time.time()
    print(round(t101-t100),'seconds')
    # return
    POL = get_painting_list(Habs)
    reached_the_bound = False
    for run_cnt in range(num_of_runs):
        for alg_cnt in [1, 2]:
            if alg_num == alg_cnt or alg_num == 3:
                for vbc_search_depth in [1, 2]:
                    print(f'>>>>>>>>>>>>>>>>>>>>>>>>#{gsetnum}: Running Alg{alg_cnt} VBC_Depth={vbc_search_depth}...>>>')
                    # if alg_cnt == 1:
                    #     alg_str = ''
                    # else:
                    #     alg_str = '_alg' + str(alg_cnt)
                    alg_str = '_alg' + str(alg_cnt)
                    test_func_called = False
                    num_proc_to_run = min(n_proc, cycles_left_to_run)
                    my_processes = []
                    for pid in range(num_proc_to_run-1):
                        cycle_id = n_proc * run_cnt + pid
                        current_inner_cycle = inner_cycles_dual[cycle_id]
                        fname_startidx = 'tmp/gset_{:d}_max_cut_mask_hex_start{:d}'.format(gsetnum, start_idxs[cycle_id])+alg_str+'.txt'
                        if cycle_id in cycles_list and not skip_all_runs and not (os.path.exists(fname_startidx) and not force_rerun):
                            test_func_called = True
                            if use_multi_proc:
                                my_processes.append(Process(target=get_weighted_signed_max_cut_incremental,
                                                            args=(gsetnum, P, start_idxs[cycle_id], current_inner_cycle, Habs, Hsign,
                                                                  short_paths_base, H_edges, odd_vertices, edges_array, full_min_cover,
                                                                  alg_cnt, vbc_search_depth, vbc_pos_odd[cycle_id])))
                            else:
                                get_weighted_signed_max_cut_incremental(gsetnum, P, start_idxs[cycle_id], current_inner_cycle, Habs, Hsign,
                                                                        short_paths_base, H_edges, odd_vertices, edges_array, full_min_cover,
                                                                        alg_cnt, vbc_search_depth, vbc_pos_odd[cycle_id])
                    if test_func_called:
                        if use_multi_proc:
                            print('running {:d} extra processes'.format(num_proc_to_run-1))
                        if use_multi_proc:
                            for proc in my_processes:
                                proc.start()

                    pid = num_proc_to_run-1  # Run one cycle the main process
                    cycle_id = n_proc * run_cnt + pid
                    current_inner_cycle = inner_cycles_dual[cycle_id]
                    fname_startidx = 'tmp/gset_{:d}_max_cut_mask_hex_start{:d}'.format(gsetnum, start_idxs[cycle_id])+alg_str+'.txt'

                    if cycle_id in cycles_list and not skip_all_runs and not (os.path.exists(fname_startidx) and not force_rerun):
                        test_func_called = True
                        get_weighted_signed_max_cut_incremental(gsetnum, P, start_idxs[cycle_id], current_inner_cycle, Habs, Hsign,
                                                            short_paths_base, H_edges, odd_vertices, edges_array, full_min_cover,
                                                                alg_cnt, vbc_search_depth, vbc_pos_odd[cycle_id])

                    if test_func_called:
                        if use_multi_proc:
                            for proc in my_processes:
                                proc.join()
                            print('multi proc finished')

                    if test_func_called:
                        MC = []
                        CN = []
                        for pid in range(num_proc_to_run):
                            cycle_id = n_proc * run_cnt + pid
                            if cycle_id in cycles_list and not skip_all_runs:
                                fname = 'tmp/gset_{:d}_max_cut_mask_hex_start{:d}'.format(gsetnum, start_idxs[cycle_id])+alg_str+'.txt'
                                if os.path.exists(fname):
                                    maxcut_mask_array = read_from_hex_file_to_bin(fname)
                                    H_maxcut = np.zeros_like(H)
                                    for ecnt in range(len(edges_array)):
                                        edge = edges_array[ecnt]
                                        if maxcut_mask_array[0, ecnt]:
                                            H_maxcut[edge[0],edge[1]] = H[edge[0],edge[1]]
                                            H_maxcut[edge[1],edge[0]] = H[edge[1],edge[0]]

                                    mc_mask_mat = (H_maxcut != 0)
                                    mc_grps_vec = paint_nodes_for_maxcut_using_pol(POL, mc_mask_mat)

                                    maxcut_mask_array_updated = np.zeros_like(maxcut_mask_array)
                                    mc_valid = 0
                                    edg_cnt = -1
                                    for edg in edges_array:
                                        edg_cnt += 1
                                        if mc_grps_vec[edg[0]] != mc_grps_vec[edg[1]]:
                                            maxcut_mask_array_updated[0, edg_cnt] = 1
                                            mc_valid += H[edg[0], edg[1]]
                                    mc_valid = int(mc_valid)
                                    if mc_valid >= mc_val_best:
                                        reset_answers = False
                                        if mc_valid > mc_val_best:
                                            reset_answers = True
                                            mc_val_best = mc_valid
                                            mc_mask_best = []
                                            mc_grps_best = []
                                        new_sol = True
                                        for x in mc_grps_best:
                                            if (all([x[nn] == mc_grps_vec[nn] for nn in range(N)]) or
                                                    all([x[nn] == 1 - mc_grps_vec[nn] for nn in range(N)])):
                                                new_sol = False
                                                break
                                        if new_sol:
                                            mc_grps_best.append(mc_grps_vec)
                                            mc_mask_best.append(maxcut_mask_array_updated[0])
                                            if reset_answers:
                                                sol_str = 'better'
                                            else:
                                                sol_str = 'new'
                                            print(f'Found {sol_str} solution: {mc_valid}. Total #solutions =',len(mc_grps_best))
                                    real_cut = int(mc_valid)
                                    if real_cut == MCLimit:
                                        reached_the_bound = True
                                else:
                                    real_cut = -2
                                MC.append(real_cut)
                                cross_num = int(np.loadtxt('info/gset_{:d}_cyclestart{:d}_crossnum.txt'.format(gsetnum, start_idxs[cycle_id])))
                                CN.append(cross_num)
                                cycle_cross_numbers[cycle_id] = cross_num
                        print(f"CutValues={MC},CrossNums={CN},[{run_cnt + 1}/{num_of_runs}])")
                    if reached_the_bound:
                        break
            if reached_the_bound:
                break
        cycles_left_to_run -= num_proc_to_run
        if reached_the_bound:
            print('>>>>> Reached the MC bound !!!!!')
            break

    print('======')

    return mc_val_best, NOddVert, NActShortPaths, mc_grps_best, mc_mask_best
