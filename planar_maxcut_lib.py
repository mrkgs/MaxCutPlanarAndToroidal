# Planar functions
import numpy as np
import networkx as nx
from file_lib import write_bin_to_hex_file


def get_path_edges(curr_path):
    ee_sorted = []
    for ii in range(len(curr_path) - 1):
        curr_edg = [curr_path[ii], curr_path[ii + 1]]
        ee_sorted.append([min(curr_edg), max(curr_edg)])
    return ee_sorted


def get_cross_num(in_cycle, max_match, min_paths, odd_vertices, Pmat):
    cross_n = 0
    for edge in in_cycle:
        for check_pair in max_match:
            ov0, ov1 = min(check_pair), max(check_pair)
            if Pmat[ov0, ov1] != 0:
                current_path = min_paths[odd_vertices[ov0], odd_vertices[ov1]]
                # current_path_edges_sorted = []
                # for ii in range(len(current_path) - 1):
                #     curr_edge = [current_path[ii], current_path[ii + 1]]
                #     current_path_edges_sorted.append([min(curr_edge), max(curr_edge)])
                current_path_edges_sorted = get_path_edges(current_path)
                if edge in current_path_edges_sorted:
                    cross_n += 1
                    break
    return cross_n


def get_min_odd_paths(Wsign, edg, fac, nei, MinPathTh):
    n_fac = len(fac)
    #Generate Dual Graph
    H_dual = np.zeros((n_fac, n_fac), dtype=bool)
    H_dual_cnt = np.zeros((n_fac, n_fac), dtype=bool)  # Number of dual edges between any two faces %2
    H_edges = np.zeros((n_fac, n_fac, 2), dtype=np.uint16)
    cnt = -1
    for dual_edge in nei:
        cnt += 1
        primal_edge_sign = Wsign[edg[cnt, 0], edg[cnt, 1]]
        if primal_edge_sign == 0:                                                                 # (positive)
            H_dual_cnt[dual_edge[0], dual_edge[1]] = 1 - H_dual_cnt[dual_edge[0], dual_edge[1]]   # +1 % 2
        w = H_dual[dual_edge[0], dual_edge[1]]
        if w == 0:
            H_dual[dual_edge[0], dual_edge[1]] = 1                                                # =|weight|
            H_edges[dual_edge[0], dual_edge[1], :] = edg[cnt, :]
            H_edges[dual_edge[1], dual_edge[0], :] = edg[cnt, :]

    H_dual += H_dual.T
    H_dual_cnt += H_dual_cnt.T
    # Find odd faces
    neigbors_num_odd = H_dual_cnt.sum(axis=0) % 2
    odd_vertices = np.arange(n_fac)[neigbors_num_odd == 1]
    n_odd_dual_vertices = len(odd_vertices)

    P = np.zeros((n_odd_dual_vertices, n_odd_dual_vertices), dtype='uint8')
    G = nx.from_numpy_array(H_dual, create_using=nx.DiGraph)
    ov_short_paths = dict()
    for cnt1 in range(n_odd_dual_vertices):
        ov1 = odd_vertices[cnt1]
        sh_paths = nx.shortest_path(G, source=ov1, method='dijkstra')
        for cnt2 in range(cnt1+1,n_odd_dual_vertices):
            ov2 = odd_vertices[cnt2]
            if len(sh_paths[ov2])<=MinPathTh:
                ov_short_paths[int(ov1),int(ov2)] = [int(vx) for vx in sh_paths[ov2]]
                tmp_len = len(ov_short_paths[ov1,ov2]) - 1
                assert tmp_len < 255
                P[cnt1, cnt2] = tmp_len

    # print('# unconnected odd nodes=',sum(np.sum(P, axis=0)+np.sum(P, axis=1)==0))

    return ov_short_paths, H_dual, H_dual_cnt, H_edges, odd_vertices, P


def get_full_ov_min_cover(P0, odd_vertices, min_paths, H_edges, Wabs, Wsign, gset_num):
    n_odd_dual_vertices = len(odd_vertices)
    maxcut_mask = (Wabs * (Wsign == 0)).astype(bool)
    if n_odd_dual_vertices>0:
        P_full = P0 + P0.T
        max_min_path_full = int(P_full.max())
        H_odd_full = (((1 + max_min_path_full) - P_full)*(P_full > 0)).astype(np.uint16)
        Gc_full = nx.from_numpy_array(H_odd_full, create_using=nx.Graph)
        max_match_full = nx.algorithms.max_weight_matching(Gc_full, maxcardinality=True, weight='weight')
        min_cover = []
        selected_lens = np.zeros((len(max_match_full)))
        match_cnt = -1
        for odd_v_pair in max_match_full:
            match_cnt += 1
            selected_lens[match_cnt] = P_full[odd_v_pair[0],odd_v_pair[1]]
            ov1, ov2 = int(odd_vertices[odd_v_pair[0]]), int(odd_vertices[odd_v_pair[1]])
            if ov1>ov2:
                ov1, ov2 = ov2, ov1
            p0 = min_paths[ov1,ov2]
            for idx in range(len(p0) - 1):
                min_cover.append(p0[idx:idx + 2])

        for dual_edge in min_cover:
            my_edge = H_edges[dual_edge[0], dual_edge[1], :]
            maxcut_mask[my_edge[0], my_edge[1]] = 1 - maxcut_mask[my_edge[0], my_edge[1]]
            maxcut_mask[my_edge[1], my_edge[0]] = 1 - maxcut_mask[my_edge[1], my_edge[0]]
        if len(selected_lens)>0:
            max_match = int(selected_lens.max())
            sum_match = int(selected_lens.sum())
        H0_maxcut = Wabs.astype(np.int8) * (1 - 2 * (Wsign.astype(np.int8))) * maxcut_mask
        max_in = int(P_full.max())
    else:
        max_match = -1
        sum_match = -1
        max_in = -1
        H0_maxcut = Wabs.astype(np.int8) * (1 - 2 * (Wsign.astype(np.int8)))
        max_match_full = []
    assert H0_maxcut.sum() % 2 == 0
    MCLimit = int(H0_maxcut.sum() // 2)
    print(f" MCBound={MCLimit}. max(|min-path|): MaxMatch={max_match},  SumMatch={sum_match} , MaxInput={max_in}.")
    return max_match_full, MCLimit


def cut_all_crossing_paths(P0, in_cycle, n_odd_dual_vertices, min_paths, odd_vertices, max_match_full):
    P1 = P0.copy()
    ecnt = -1
    crossed_ecnt = -1
    for in_cycle_edge in in_cycle:  # Cut all crossing paths
        ecnt += 1
        crossed_ecnt += 1
        for ov1_cnt in range(n_odd_dual_vertices):
            for ov2_cnt in range(ov1_cnt + 1, n_odd_dual_vertices):
                if P1[ov1_cnt, ov2_cnt] != 0:
                    current_path = min_paths[odd_vertices[ov1_cnt], odd_vertices[ov2_cnt]]
                    current_path_edges_sorted = get_path_edges(current_path)
                    if in_cycle_edge in current_path_edges_sorted:
                        P1[ov1_cnt, ov2_cnt] = 0
    check_cross_num = get_cross_num(in_cycle, max_match_full, min_paths, odd_vertices, P1)
    assert check_cross_num == 0
    return P1

def alg2(n_odd_dual_vertices, P1_in, min_paths, odd_vertices, in_cycle):
    P1 = P1_in.copy()
    for ov1_cnt in range(n_odd_dual_vertices):
        for ov2_cnt in range(ov1_cnt + 1, n_odd_dual_vertices):
            if P1[ov1_cnt, ov2_cnt] != 0:
                current_path = min_paths[odd_vertices[ov1_cnt], odd_vertices[ov2_cnt]]
                current_path_edges = [[current_path[ii], current_path[ii + 1]] for ii in range(len(current_path) - 1)]
                current_path_edges_sorted = [[min(x), max(x)] for x in current_path_edges]
                for in_cycle_edge in in_cycle:
                    if in_cycle_edge in current_path_edges_sorted:
                        P1[ov1_cnt, ov2_cnt] = 0
                        break
    P1 = P1 + P1.T
    max_match, max_min_path = [], -1
    if n_odd_dual_vertices > 0:
        max_min_path = int(P1.max())
        H_odd = (((1 + max_min_path) - P1) * (P1 > 0)).astype(np.uint8)
        Gc = nx.from_numpy_array(H_odd, create_using=nx.Graph)
        max_match = nx.algorithms.max_weight_matching(Gc, maxcardinality=True, weight='weight')
    return max_match, max_min_path

def get_weighted_signed_max_cut_incremental(gsetnum, P0, cycle_start, in_cycle, Wabs, Wsign, min_paths, H_edges, odd_vertices, edges_array, max_match_full, alg_num, check_vbc_mode, cycle_pos_odd):
    # if alg_num == 1:
    #     alg_str = ''
    # else:
    #     alg_str = '_alg'+str(alg_num)
    valid_cycle = False
    alg_str = '_alg'+str(alg_num)
    count_cross_num = True
    n_odd_dual_vertices = len(odd_vertices)
    # P1 = P0.copy()
    if count_cross_num:
        cross_num = 0
        if n_odd_dual_vertices>0:
            if len(max_match_full)*2 == n_odd_dual_vertices:
                cycle_cross_array = np.zeros(len(in_cycle))
                cross_pairs = []
                ecnt = -1
                for in_cycle_edge in in_cycle:
                    ecnt += 1
                    for pair in max_match_full:
                        ov0, ov1 = min(pair), max(pair)
                        current_path = min_paths[odd_vertices[ov0],odd_vertices[ov1]]
                        current_path_edges_sorted = []
                        for ii in range(len(current_path) - 1):
                            curr_edge = [current_path[ii], current_path[ii+1]]
                            current_path_edges_sorted.append([min(curr_edge), max(curr_edge)])
                        if in_cycle_edge in current_path_edges_sorted:
                            cross_num += 1
                            cycle_cross_array[ecnt] = len(current_path)-1
                            cross_pairs.append([ov0, ov1])
                            break
        np.savetxt('info/gset_{:d}_cyclestart{:d}_crossnum.txt'.format(gsetnum, cycle_start), np.atleast_1d(cross_num))
        # return
    max_idx = -1
    max_idx_2 = -1
    if n_odd_dual_vertices > 0:
        if any(cycle_cross_array>0):
            cors_len = cycle_cross_array[np.where(cycle_cross_array > 0)]
            max_idx = np.argmax(cors_len).astype(int)
            if len(cors_len)>1:
                cors_len[max_idx] = 0
                max_idx_2 = np.argmax(cors_len).astype(int)
    maxcut_mask = (Wabs*(Wsign==0)).astype(bool)
    maxcut_mask_change_flag = np.zeros_like(maxcut_mask, dtype=bool)
    # cross_num = 0
    cross_odd = cross_num % 2
    cross_even = not cross_odd
    cycle_pos_even = not cycle_pos_odd

    use_alg2_as_backup = 0
    if alg_num==1:
        if (cross_even and cycle_pos_even) or (cross_odd and cycle_pos_odd):  # rm two crossing paths
            if check_vbc_mode==1 or cross_num < 2:
                max_match = max_match_full
            else:  #check_vbc_mode==2 and cross_num>0
                P1 = cut_all_crossing_paths(P0, in_cycle, n_odd_dual_vertices, min_paths, odd_vertices, max_match_full)
                pcnt = -1
                for cross_p in cross_pairs:
                    pcnt += 1
                    assert (max_idx>=0 and max_idx_2>=0)
                    if pcnt not in [max_idx, max_idx_2]:  # fix all crossing shortest paths except two longest ones
                        P1[cross_p[0],:] = 0
                        P1[cross_p[1],:] = 0
                        P1[:,cross_p[0]] = 0
                        P1[:,cross_p[1]] = 0
                        P1[cross_p[0], cross_p[1]] = P0[cross_p[0], cross_p[1]]
                P1 = P1 + P1.T
                if n_odd_dual_vertices>0:
                    max_min_path = int(P1.max())
                    H_odd = (((1 + max_min_path) - P1)*(P1 > 0)).astype(np.uint8)
                    Gc = nx.from_numpy_array(H_odd, create_using=nx.Graph)
                    max_match = nx.algorithms.max_weight_matching(Gc, maxcardinality=True, weight='weight')

                check_cross_num = get_cross_num(in_cycle, max_match, min_paths, odd_vertices, P1)
                if cross_even and (check_cross_num % 2 == 1) or cross_odd and (check_cross_num % 2 == 0):
                    print('Cycle not found, using backup procedure (alg2)', check_cross_num)
                    use_alg2_as_backup = 1
            if n_odd_dual_vertices > 0:
                max_min_path = int(P0.max())
            else:
                max_min_path = -1
        else:  #(cross_even and cycle_pos_odd) or (cross_odd and cycle_pos_even)
            if cross_num==0:  #cycle_pos_odd, but there are no crossing paths to remove, leaving it as is
                max_match = max_match_full
            else:
                P1 = cut_all_crossing_paths(P0, in_cycle, n_odd_dual_vertices, min_paths, odd_vertices, max_match_full)

                if check_vbc_mode == 1 or cross_num < 2:
                    assert max_idx >= 0
                    cut_idx = max_idx
                else:  # check_vbc_mode == 2 and cross_num > 1
                    assert max_idx_2 >= 0
                    cut_idx = max_idx_2

                pcnt = -1
                for cross_p in cross_pairs:
                    pcnt += 1
                    if pcnt != cut_idx:  # fix all crossing paths except the longest one
                        P1[cross_p[0],:] = 0
                        P1[cross_p[1],:] = 0
                        P1[:,cross_p[0]] = 0
                        P1[:,cross_p[1]] = 0
                        P1[cross_p[0], cross_p[1]] = P0[cross_p[0], cross_p[1]]

                P1 = P1 + P1.T
                if n_odd_dual_vertices>0:
                    max_min_path = int(P1.max())
                    H_odd = (((1 + max_min_path) - P1)*(P1 > 0)).astype(np.uint8)
                    Gc = nx.from_numpy_array(H_odd, create_using=nx.Graph)
                    max_match = nx.algorithms.max_weight_matching(Gc, maxcardinality=True, weight='weight')

                check_cross_num = get_cross_num(in_cycle, max_match, min_paths, odd_vertices, P1)
                if cross_even and (check_cross_num % 2 == 0) or cross_odd and (check_cross_num % 2 == 1):
                    print('Cycle not found, using backup procedure (alg2)', check_cross_num)
                    use_alg2_as_backup = 1
    if n_odd_dual_vertices > 0:
        max_min_path = int(P0.max())
    else:
        max_min_path = -1

    if alg_num==2 or use_alg2_as_backup:
        # as a backup of alg1 it continues from the alg1 matrix, as alg2 itself it starts from the full one
        max_match, max_min_path = alg2(n_odd_dual_vertices, P1 if use_alg2_as_backup else P0, min_paths, odd_vertices, in_cycle)

    if n_odd_dual_vertices>0:
        if len(max_match)*2 == n_odd_dual_vertices:
            valid_cycle = True
        else:
            print(f"cycle start={cycle_start} len(max_match)x2={len(max_match)*2} #ovpairs={n_odd_dual_vertices}")
    else:
        valid_cycle = True
        max_match = []
        max_min_path = -1
    if valid_cycle:
        min_cover = []
        selected_lens = np.zeros((len(max_match)))
        match_cnt = -1
        for odd_v_pair in max_match:
            match_cnt += 1
            selected_lens[match_cnt] = P0[odd_v_pair[0],odd_v_pair[1]]
            ov1, ov2 = int(odd_vertices[odd_v_pair[0]]), int(odd_vertices[odd_v_pair[1]])
            if ov1>ov2:
                ov1, ov2 = ov2, ov1
            p0 = min_paths[ov1,ov2]
            for idx in range(len(p0) - 1):
                min_cover.append(p0[idx:idx + 2])

        min_cover_sorted = [[min(x), max(x)] for x in min_cover]
        if alg_num==2:
            for in_cycle_edge in in_cycle:
                assert not (in_cycle_edge in min_cover_sorted)

        for dual_edge in min_cover:

            my_edge = H_edges[dual_edge[0], dual_edge[1], :]
            maxcut_mask[my_edge[0], my_edge[1]] = 1 - maxcut_mask[my_edge[0], my_edge[1]]
            maxcut_mask[my_edge[1], my_edge[0]] = 1 - maxcut_mask[my_edge[1], my_edge[0]]
            if maxcut_mask_change_flag[my_edge[0], my_edge[1]]==1:
                tmp=1
            maxcut_mask_change_flag[my_edge[0], my_edge[1]] = 1
            maxcut_mask_change_flag[my_edge[1], my_edge[0]] = 1
        if len(selected_lens)>0:
            max_match = int(selected_lens.max())
            sum_match = int(selected_lens.sum())
        else:
            sum_match = -1
    else:
        max_match = -1
        sum_match = -1
        max_min_path = -1
    # Obtain Max Cut Nodes Groups From Set of Max Cut Edges using DFS
    H0_maxcut = Wabs.astype(np.int8) * (1 - 2*(Wsign.astype(np.int8))) * maxcut_mask

    print(f"cycle start {cycle_start}. {['EVEN','ODD '][cycle_pos_odd]}. PreliminaryValue= {int(H0_maxcut.sum()//2)}. max(|min-path|): max match={max_match},  sum match={sum_match} , max input={max_min_path}.",in_cycle)
    maxcut_mask_array = np.zeros((1,len(edges_array)))
    for ecnt in range(len(edges_array)):
        maxcut_mask_array[0,ecnt] = maxcut_mask[edges_array[ecnt,0], edges_array[ecnt,1]]
    fname = 'tmp/gset_{:d}_max_cut_mask_hex_start{:d}'.format(gsetnum, cycle_start)+alg_str+'.txt'
    write_bin_to_hex_file(maxcut_mask_array, fname, "w")
    return


def get_painting_list(A):
    PaintOrderList = []
    N = len(A)
    painted = np.zeros(N, dtype=int)
    visited = np.zeros(N, dtype=int)
    painted[0] = 1
    unpainted_nodes_num = N - 1
    for line_cnt in range(N):
        if unpainted_nodes_num==0:
            break
        for trycnt in range(N):
            if unpainted_nodes_num==0:
                break
            PaintLine = []
            painted_more_nodes = False
            for cnt0 in range(N):
                if painted[cnt0] and not visited[cnt0]:  # Find unvisited node
                    visited[cnt0] = 1
                    PaintLine.append(cnt0)
                    connected_ind = [int(x) for x in np.argwhere(A[cnt0, :] != 0)]
                    for cnt1 in connected_ind:
                        if not painted[cnt1]:            # Find unpainted node
                            PaintLine.append(cnt1)
                            painted[cnt1] = 1
                            unpainted_nodes_num -= 1
                            painted_more_nodes = True
                if painted_more_nodes:
                    break
            assert painted_more_nodes
            PaintOrderList.append(PaintLine)
            if unpainted_nodes_num % 100 == 0:
                print('Building Painting Order:',unpainted_nodes_num)
        assert unpainted_nodes_num==0
    return PaintOrderList


def verify_nodes_for_maxcut_with_edges(node_grp, edges_arr, mc_mask):
    check_ok = True
    ecnt = -1
    for myedge in edges_arr:
        ecnt += 1
        if mc_mask[0,ecnt]==1:
            if node_grp[myedge[0]] == node_grp[myedge[1]]:
                check_ok = False
        else:
            if node_grp[myedge[0]] != node_grp[myedge[1]]:
                check_ok = False
        check_ok = (node_grp[myedge[1]] == node_grp[myedge[0]] ^ mc_mask[0,ecnt])
        if not check_ok:
            # print(ecnt, myedge, mc_mask[0,ecnt], node_grp[myedge])
            break
    return check_ok


def paint_nodes_for_maxcut_using_pol(pol, maxcut_mask_matrix):
    N = len(maxcut_mask_matrix)
    gg = np.zeros(N, dtype=int)
    for pol_line in pol:
        v0 = pol_line[0]
        for v1 in pol_line[1:]:
            gg[v1] = gg[v0] ^ maxcut_mask_matrix[v0,v1]
    return gg
