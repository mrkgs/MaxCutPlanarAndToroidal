from planar_maxcut_lib import get_min_odd_paths, get_full_ov_min_cover, get_weighted_signed_max_cut_incremental, paint_nodes_for_maxcut_using_pol
from planar_maxcut_lib import get_painting_list, paint_nodes_for_maxcut_simple, verify_nodes_for_maxcut_with_edges
from toroid_faces_lib import find_edges_and_faces_of_toroid
from toroid_cycles_lib import find_inner_cycle_path_of_toroid
import time, os
import numpy as np
import pickle
import networkx as nx
from multiprocessing import Process
from file_lib import read_from_hex_file_to_bin
from graph_lib import get_clock_wise_dir_edges
from params_config import use_multi_proc, ncycles

def get_toroid_max_cut_bound(H, Habs, Hsign, D, gsetnum):
    MinPathThreshhold = 5
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
    np.savetxt(f'./bounds/gset{gsetnum}_min_cover_time.txt', np.atleast_1d(t101))
    np.savetxt(f'./bounds/gset{gsetnum}_max_cut_bound.txt', np.atleast_1d(MCLimit))
    exit(0)

def get_toroid_max_cut_topological(H, Habs, Hsign, D, print_en, gsetnum, cycles_starts_to_check, force_rerun, alg_num):
    if not os.path.exists('./info'):
        os.makedirs('./info')
    if not os.path.exists('./bounds'):
        os.makedirs('./bounds')
    if not os.path.exists('./tmp'):
        os.makedirs('./tmp')
    if alg_num == 1:
        alg_str = ''
    else:
        alg_str = '_alg'+str(alg_num)
    print_str = 'GSet'+str(gsetnum)
    t0 = time.time()
    N = len(H)
    MinPathThreshhold = int(0.09*N/D)
    MinPathThreshhold = 5
    # MinPathThreshhold = 30000
    # MinPathThreshhold = 100
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

    # print('# unconnected odd nodes before threshold=',sum(np.sum(P+P.T, axis=0)==0))
    # P = P * (P < MinPathThreshhold)
    # print('# unconnected odd nodes after threshold=',sum(np.sum(P+P.T, axis=0)==0))
    # np.save('P_full.npy',P)
    # assert sum(np.sum(P+P.T, axis=0)==0) == 0
    # assert 0
    max_cut_log = []
    max_cut_val_J_best = 0
    inner_cycles = []
    inner_cycles_e = []
    inner_cycles_dual = []
    start_idxs = []
    force_check_cycles = 1
    # for in_cycle_0 in range(0,NCycles2Try,1):
    # for in_cycle_0 in range(0,N,1):
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

        if len(inner_cycle_path) == 0:
            valid_inner_cycle = False
        else:
            pos_edges_cnt = 0
            path_len = len(inner_cycle_path)
            for vcnt in range(path_len):
                if Hsign[inner_cycle_path[vcnt], inner_cycle_path[(vcnt+1)%path_len]] == 0: # if positive
                    pos_edges_cnt+= 1
            valid_inner_cycle = (pos_edges_cnt % 2 == 0)  # Take only even cycles

        if valid_inner_cycle:
            inner_cycles.append(inner_cycle_path)
            inner_cycles_dual.append(inner_cycle_dual_edges)
            inner_cycles_e.append(inner_cycle_edges)
            print(inner_cycle_path)
            start_idxs.append(in_cycle_0)
        num_inner_cycles = len(inner_cycles)
        if num_inner_cycles == ncycles:
            break
    cycle_cross_numbers = np.zeros(num_inner_cycles)
    print(f"GSet {gsetnum}, #Nodes={N}, #Edges={len(edges_array)}: found {num_inner_cycles} inner even cycles")
    n_proc = 6
    # use_multi_proc = 0
    num_of_runs = int(np.ceil(num_inner_cycles/n_proc))
    cycles_left_to_run = num_inner_cycles
    cycles_list = np.arange(num_inner_cycles)
    # cycles_list = [0, 1, 2, 3, 4, 5]
    skip_all_runs = 0
    # force_rerun = 0
    # np.save('Puint16.npy',P.astype(np.uint16))
    print(f'#OV={len(odd_vertices)}')
    print('Starting searching ov min cover...')
    t100 = time.time()
    full_min_cover, MCLimit = get_full_ov_min_cover(P.astype(np.uint16), odd_vertices, short_paths_base, H_edges, Habs, Hsign, gsetnum)
    t101 = time.time()
    print(round(t101-t100),'seconds')
    # return
    if 1:
        for run_cnt in range(num_of_runs):
            test_func_called = False
            num_proc_to_run = min(n_proc, cycles_left_to_run)
            cycles_left_to_run -= num_proc_to_run
            my_processes = []
            for pid in range(num_proc_to_run-1):
                cycle_id = n_proc * run_cnt + pid
                # with open("log.txt", "a") as flog:
                #     flog.write(f"cycle {cycle_id}/{len(inner_cycles)}\n")

                current_inner_cycle = inner_cycles_dual[cycle_id]
                # current_inner_cycle_2 = inner_cycles_dual_2[cycle_id]
                # fname = 'tmp/gset_{:d}_max_cut_mask_{:d}.txt'.format(gsetnum, cycle_id)
                # if not use_alg2:
                #     fname_startidx = 'tmp/gset_{:d}_max_cut_mask_hex_start{:d}.txt'.format(gsetnum, start_idxs[cycle_id])
                # else:
                #     fname_startidx = 'tmp/gset_{:d}_max_cut_mask_hex_start{:d}_alg2.txt'.format(gsetnum, start_idxs[cycle_id])
                fname_startidx = 'tmp/gset_{:d}_max_cut_mask_hex_start{:d}'.format(gsetnum, start_idxs[cycle_id])+alg_str+'.txt'
                # if os.path.exists(fname):
                #     shutil.move(fname, fname_startidx)
                if cycle_id in cycles_list and not skip_all_runs and not (os.path.exists(fname_startidx) and not force_rerun):
                    test_func_called = True
                    if use_multi_proc:
                        my_processes.append(Process(target=get_weighted_signed_max_cut_incremental,
                                                    args=(gsetnum, P, start_idxs[cycle_id], current_inner_cycle, Habs, Hsign,
                                                          short_paths_base, H_edges, odd_vertices, edges_array, full_min_cover, alg_num)))
                    else:
                        get_weighted_signed_max_cut_incremental(gsetnum, P, start_idxs[cycle_id], current_inner_cycle, Habs, Hsign,
                                                                short_paths_base, H_edges, odd_vertices, edges_array, full_min_cover, alg_num)
            if test_func_called:
                if use_multi_proc:
                    print('running {:d} extra processes'.format(num_proc_to_run-1))
                if use_multi_proc:
                    for proc in my_processes:
                        proc.start()

            pid = num_proc_to_run-1  # Run one cycle the main process
            cycle_id = n_proc * run_cnt + pid
            # print(f"cycle {cycle_id}/{len(inner_cycles)}")
            # with open("log.txt", "a") as flog:
            #     flog.write(f"cycle {cycle_id}/{len(inner_cycles)}\n")

            current_inner_cycle = inner_cycles_dual[cycle_id]
            # fname = 'tmp/gset_{:d}_max_cut_mask_{:d}.txt'.format(gsetnum, cycle_id)
            # if not use_alg2:
            #     fname_startidx = 'tmp/gset_{:d}_max_cut_mask_hex_start{:d}.txt'.format(gsetnum, start_idxs[cycle_id])
            # else:
            #     fname_startidx = 'tmp/gset_{:d}_max_cut_mask_hex_start{:d}_alg2.txt'.format(gsetnum, start_idxs[cycle_id])
            fname_startidx = 'tmp/gset_{:d}_max_cut_mask_hex_start{:d}'.format(gsetnum, start_idxs[cycle_id])+alg_str+'.txt'
            # if os.path.exists(fname):
            #     shutil.move(fname, fname_startidx)

            if cycle_id in cycles_list and not skip_all_runs and not (os.path.exists(fname_startidx) and not force_rerun):
                test_func_called = True
                get_weighted_signed_max_cut_incremental(gsetnum, P, start_idxs[cycle_id], current_inner_cycle, Habs, Hsign,
                                                    short_paths_base, H_edges, odd_vertices, edges_array, full_min_cover, alg_num)

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
                        # max_cut_grps_J = np.loadtxt('tmp/groups_{:d}.txt'.format(cycle_id))
                        # J_maxcut = np.loadtxt('tmp/gset_{:d}_max_cut_matrix_{:d}.txt'.format(gsetnum, cycle_id))
                        # J_maxcut_mask = np.loadtxt('tmp/gset_{:d}_max_cut_mask_{:d}.txt'.format(gsetnum, cycle_id))
                        # J_maxcut_mask_old = np.loadtxt('tmp/gset_{:d}_max_cut_mask_start{:d}.txt'.format(gsetnum, start_idxs[cycle_id]))
                        # J_maxcut_mask = read_from_hex_file_to_bin('tmp/gset_{:d}_max_cut_mask_start{:d}_hex.txt'.format(gsetnum, start_idxs[cycle_id]))
                        fname = 'tmp/gset_{:d}_max_cut_mask_hex_start{:d}'.format(gsetnum, start_idxs[cycle_id])+alg_str+'.txt'
                        if os.path.exists(fname):
                            maxcut_mask_array = read_from_hex_file_to_bin(fname)
                            H_maxcut = np.zeros_like(H)
                            for ecnt in range(len(edges_array)):
                                edge = edges_array[ecnt]
                                if maxcut_mask_array[0, ecnt] == True:
                                    H_maxcut[edge[0],edge[1]] = H[edge[0],edge[1]]
                                    H_maxcut[edge[1],edge[0]] = H[edge[1],edge[0]]

                            # H_maxcut = np.zeros_like(H)
                            # for i1 in range(N):
                            #     for i2 in range(N):
                            #         if J_maxcut[i1, i2] != 0:
                            #             H_maxcut[i1, i2] = H[i1, i2]

                            # real_cut = get_cut(edges_array, max_cut_grps_J, H_maxcut)
                            real_cut = int(H_maxcut.sum() // 2)
                        else:
                            real_cut = -2
                        MC.append(real_cut)
                        cross_num = int(np.loadtxt('info/gset_{:d}_cyclestart{:d}_crossnum.txt'.format(gsetnum, start_idxs[cycle_id])))
                        CN.append(cross_num)
                        cycle_cross_numbers[cycle_id] = cross_num

                print(f"PreliminaryValue={MC},CrossNum={CN},[{run_cnt+1}/{num_of_runs}])")

    # max_cut_val_best, max_cut_log, max_cut_masks_best, max_cut_grps_best = get_toroid_max_cut_direct_painting(H, Habs, Hsign, D,0, gsetnum, cycles_starts_to_check)
    # return max_cut_val_J_best, max_cut_log, max_cut_masks_best, max_cut_grps_best
    print('======')
    # print(cycle_cross_numbers, 'min:', min(cycle_cross_numbers), 'max:', max(cycle_cross_numbers))

    return start_idxs


def get_toroid_max_cut_direct_painting(H, Habs, Hsign, D, print_en, gsetnum, n_cycles_2_try, alg_num):
    if alg_num == 1:
        alg_str = ''
    else:
        alg_str = '_alg'+str(alg_num)
    print_str = 'GSet'+str(gsetnum)
    t0 = time.time()
    N = len(H)
    MinPathThreshhold = int(0.09*N/D)
    MinPathThreshhold = 6
    force_find_faces = 1
    save_faces = 0
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
            np.savetxt('info/neighbors' + str(gsetnum) + '.txt', neighbor_faces_array)

    max_cut_log = []
    max_cut_val_J_best = 0
    inner_cycles = []
    inner_cycles_e = []
    inner_cycles_dual = []
    start_idxs = []
    for in_cycle_0 in n_cycles_2_try: #range(0,N,1):
        cycle_file1 = 'info/'+str(gsetnum)+'_cycle_'+str(in_cycle_0)+'_nodes.txt'
        cycle_file2 = 'info/'+str(gsetnum)+'_cycle_'+str(in_cycle_0)+'_duals.txt'
        cycle_file3 = 'info/'+str(gsetnum)+'_cycle_'+str(in_cycle_0)+'_edges.txt'
        if not os.path.exists(cycle_file3):
            inner_cycle_path, inner_cycle_dual_edges, inner_cycle_edges = find_inner_cycle_path_of_toroid(Habs, edges_array, neighbor_faces_array, in_cycle_0)
            if len(inner_cycle_path)>0 and len(inner_cycle_dual_edges)>0 and len(inner_cycle_edges)>0:
                np.savetxt(cycle_file1, inner_cycle_path)
                np.savetxt(cycle_file2, inner_cycle_dual_edges)
                np.savetxt(cycle_file3, inner_cycle_edges)
        else:
            inner_cycle_path = [int(x) for x in np.loadtxt(cycle_file1)]
            inner_cycle_dual_edges = [[int(x[0]), int(x[1])]  for x in np.loadtxt(cycle_file2)]
            inner_cycle_edges = [[int(x[0]), int(x[1])]  for x in np.loadtxt(cycle_file3)]

        if len(inner_cycle_path) == 0:
            valid_inner_cycle = False
        else:
            pos_edges_cnt = 0
            path_len = len(inner_cycle_path)
            for vcnt in range(path_len):
                if Hsign[inner_cycle_path[vcnt], inner_cycle_path[(vcnt+1)%path_len]] == 0: # if positive
                    pos_edges_cnt+= 1
            valid_inner_cycle = (pos_edges_cnt % 2 == 0)  # Take only even cycles

        if valid_inner_cycle:
            inner_cycles.append(inner_cycle_path)
            inner_cycles_dual.append(inner_cycle_dual_edges)
            inner_cycles_e.append(inner_cycle_edges)
            print(inner_cycle_path)
            start_idxs.append(in_cycle_0)

    num_inner_cycles = len(inner_cycles)
    print(f"GSet {gsetnum}, #Nodes={N}, #Edges={len(edges_array)}: found {num_inner_cycles} inner even cycles")

    max_cut_masks_best = []
    max_cut_grps_best = []
    paint_en = 1

    failed_check_cnt = 0
    POL = get_painting_list(Habs)
    for cycle_id in range(num_inner_cycles):
        # if not use_alg2:
        #     fname = 'tmp/gset_{:d}_max_cut_mask_hex_start{:d}.txt'.format(gsetnum, start_idxs[cycle_id])
        # else:
        fname = 'tmp/gset_{:d}_max_cut_mask_hex_start{:d}'.format(gsetnum, start_idxs[cycle_id])+alg_str+'.txt'
        if os.path.exists(fname):
            if cycle_id==13:
                tmp=0
            maxcut_mask_array = read_from_hex_file_to_bin(fname)
            H_maxcut = np.zeros_like(H)
            for ecnt in range(len(edges_array)):
                edge = edges_array[ecnt]
                if maxcut_mask_array[0, ecnt] == True:
                    H_maxcut[edge[0], edge[1]] = H[edge[0], edge[1]]
                    H_maxcut[edge[1], edge[0]] = H[edge[1], edge[0]]
            maxcut_mask_matrix = H_maxcut!=0
            rm_tmp_files = 0
            if rm_tmp_files:
                os.remove('tmp/max_cut_matrix_{:d}.txt'.format(cycle_id))

            edges_array_dir = get_clock_wise_dir_edges(edges_array)

            assert  ((H_maxcut!=H) * (H_maxcut!=0)).sum() == 0
            if paint_en:

                curr_cycle = inner_cycles[cycle_id]

                max_cut_grps_J = paint_nodes_for_maxcut_using_pol(POL, maxcut_mask_matrix)
                maxcut_mask_array_updated = np.zeros_like(maxcut_mask_array)
                mc_valid = 0
                edg_cnt = -1
                for edg in edges_array:
                    edg_cnt += 1
                    if max_cut_grps_J[edg[0]] != max_cut_grps_J[edg[1]]:
                        maxcut_mask_array_updated[0, edg_cnt] = 1
                        mc_valid += H[edg[0], edg[1]] 

                verify_max_cut_validity = 1
                if verify_max_cut_validity:
                    Check_Ok = verify_nodes_for_maxcut_with_edges(max_cut_grps_J, edges_array, maxcut_mask_array_updated)
                else:
                    Check_Ok = 1

                if Check_Ok:
                    # real_cut = get_cut_simple(max_cut_grps_J, H_maxcut)
                    # assert real_cut == int(H_maxcut.sum()//2)
                    real_cut = int(H_maxcut.sum() // 2)
                    real_cut = mc_valid
                else:
                    print('Check_Ok = 0')
                    exit(0)
                    failed_check_cnt += 1
                    real_cut = -1
            else:
                real_cut = int(H_maxcut.sum()//2)

            if real_cut >= max_cut_val_J_best:
                if real_cut > max_cut_val_J_best:
                    max_cut_val_J_best = real_cut
                    max_cut_masks_best = []
                    max_cut_grps_best = []
                if paint_en:
                    if Check_Ok:
                        new_sol = True
                        for x in max_cut_grps_best:
                            if all([x[nn]==max_cut_grps_J[nn] for nn in range(N)]) or all([x[nn]==1-max_cut_grps_J[nn] for nn in range(N)]):
                                new_sol = False
                        if new_sol:
                            max_cut_grps_best.append(max_cut_grps_J)
                            max_cut_masks_best.append(maxcut_mask_array[0])
                            print('Found new solution =', new_sol, '. Total #solutions =', len(max_cut_grps_best))
                else:
                    max_cut_masks_best.append(maxcut_mask_array)
            max_cut_log.append(int(real_cut))
            nmax = sum([max_cut_log[n] == max_cut_val_J_best for n in range(len(max_cut_log))])
            print(f"{print_str} cycle#{cycle_id} start={curr_cycle[0]} #max={nmax}/{len(max_cut_log)-failed_check_cnt}[{num_inner_cycles}], Cut={real_cut}, Best={max_cut_val_J_best}")
    if paint_en:
        print(f"failed_check_cnt={failed_check_cnt}/{num_inner_cycles}")
    print('Painting time=',int(time.time()-t0),'sec., MaxCut=', max_cut_val_J_best)
    print('max_cut_log=', max_cut_log)
    if paint_en:
        assert len(max_cut_masks_best)  == len(max_cut_grps_best)
    return max_cut_val_J_best, max_cut_log, max_cut_masks_best, max_cut_grps_best