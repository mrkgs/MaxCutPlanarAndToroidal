# Library of functions for finding faces and edges
import numpy as np
import networkx as nx
import time
from toroid_cycles_lib import cut_toroid

def add_face(edge, emb, faces_lst, neighbors_mat, skip_faces):
    face0 = [np.uint16(node) for node in np.array(emb.traverse_face(edge[0], edge[1]))]
    face0_sorted = np.sort(face0)
    skip_en = False
    for skip_face in skip_faces:
        if np.array_equal(face0_sorted, skip_face):
            skip_en = True
            break
    if not skip_en:
        faces_lst.append(face0)
        fcnt = len(faces_lst)
        flen = len(face0)
        for vcnt in range(flen):
            v0, v1 = face0[vcnt], face0[(vcnt + 1) % flen]
            neighbors_mat[v0, v1] = fcnt

def get_neighbor_faces_array(edges_array, neighbor_faces_mat):
    neighbor_faces_arr = -np.ones((len(edges_array), 2), dtype=np.uint16)
    ecnt = -1
    for edge in edges_array:
        ecnt += 1
        neighbor_faces_arr[ecnt, :] = np.sort([int(neighbor_faces_mat[edge[0], edge[1]]), int(neighbor_faces_mat[edge[1], edge[0]])])
    return neighbor_faces_arr

def find_planar_faces(emb, neighbor_faces_mat, faces_lst, skip_faces):
    emb_edges = emb.edges()
    n_edges = len(emb_edges)//2
    edges = -np.ones((n_edges,2),dtype=int)
    ecnt = 0
    # find_edges_t0 = time.time()
    for edge in emb_edges:
        if edge[0]<edge[1]:
            edges[ecnt, :] = edge
            ecnt += 1
    # print('unique edges time:', int(time.time() - find_edges_t0), 'sec.')

    # find_faces_t0 = time.time()
    for edge in edges:
        if neighbor_faces_mat[edge[0], edge[1]]==0:
            add_face([edge[0], edge[1]], emb, faces_lst, neighbor_faces_mat, skip_faces)
        if neighbor_faces_mat[edge[1], edge[0]]==0:
            add_face([edge[1], edge[0]], emb, faces_lst, neighbor_faces_mat, skip_faces)

    # print('unique faces time:', int(time.time() - find_faces_t0), 'sec.')

    return


def find_edges_and_faces(H, D=0):    # supports planar and toroidal graphs
    L = len(H)
    G = nx.from_numpy_array(H, create_using=nx.DiGraph)
    is_planar, planar_embedding = nx.check_planarity(G)
    if not is_planar:  # toroid
        n_edges_total = int(sum(sum(H!=0))//2)
        H_1 = H.copy()
        i0,i1 = 0,D-1
        # i0,i1 = 1,1+D-1
        for idx in np.arange(i0,i1+1):
            # assert H_1[idx, idx + D] != 0
            H_1[idx, idx + D] = 0
            H_1[idx + D, idx] = 0
            H_1[idx, idx + D - 1] = 0
            H_1[idx + D - 1, idx] = 0

        # assert H_1[i1, i1+1] != 0
        H_1[i1, i1+1] = 0
        H_1[i1+1, i1] = 0
        G_1 = nx.from_numpy_array(H_1, create_using=nx.DiGraph)
        is_planar, planar_embedding = nx.check_planarity(G_1)
        assert is_planar
        edges_1, faces_1, neighbor_faces_1 = find_planar_faces(planar_embedding, 0)
        cut1 = np.sort((i0 - 1 + np.arange(D + 1)) % L)
        cut2 = np.sort((i1 + 1 + np.arange(D + 1)) % L)
        faces = []
        for f1 in faces_1:
            # When toroid is transformed to sphere by a cut, there are two additional faces, which are removed here
            is_cut_1 = np.array_equal(np.sort(f1), cut1)
            is_cut_2 = np.array_equal(np.sort(f1), cut2)
            if (not is_cut_1) and (not is_cut_2) and len(f1)<=7:
                faces.append(f1)

        H_1 = H.copy()
        i0,i1 = L//2,L//2+D-1
        for idx in np.arange(i0,i1+1):
            H_1[idx, idx + D] = 0
            H_1[idx + D, idx] = 0
            H_1[idx, idx + D - 1] = 0
            H_1[idx + D - 1, idx] = 0
        H_1[i1, i1+1] = 0
        H_1[i1+1, i1] = 0
        G_1 = nx.from_numpy_array(H_1, create_using=nx.DiGraph)
        is_planar, planar_embedding = nx.check_planarity(G_1)
        assert is_planar
        edges_2, faces_2, neighbor_faces_2 = find_planar_faces(planar_embedding, 0)

        # edges = np.concatenate((edges_1, -np.ones((D+1,2),dtype=int)), axis=0)
        edges = np.concatenate((edges_1, -np.ones((n_edges_total-len(edges_1),2),dtype=int)), axis=0)
        edges_total_cnt = len(edges_1)
        for ecnt_new in np.arange(len(edges_2)):
            edge_new = edges_2[ecnt_new,:]
            edge_found = False
            for ecnt in np.arange(len(edges_1)):
                if np.array_equal(edge_new, edges_1[ecnt,:]):
                    edge_found = True
                    break
            if not edge_found:
                edges[edges_total_cnt,:] = edge_new
                edges_total_cnt += 1
        # assert edges_total_cnt == len(edges_1) + D + 1, 'wrong number of edges for toroid embedding'
        assert edges_total_cnt == n_edges_total, 'wrong number of edges for toroid embedding'

        cut1 = np.sort((i0 - 1 + np.arange(D + 1)) % L)
        cut2 = np.sort((i1 + 1 + np.arange(D + 1)) % L)
        for f2 in faces_2:
            is_cut_1 = np.array_equal(np.sort(f2), cut1)
            is_cut_2 = np.array_equal(np.sort(f2), cut2)
            face_found = is_cut_1 or is_cut_2 or len(f2)>7
            if not face_found:
                for f1 in faces_1:
                    if np.array_equal(np.sort(f1), np.sort(f2)):
                        face_found = True
                        break
            if not face_found:
                faces.append(f2)
        if len(faces)<100:
            assert len(faces) == len(faces_1) - 2 + D + 1, 'wrong number of faces for toroid embedding'

        neighbor_faces = -np.ones((len(edges), 2), dtype=int)
        ecnt = -1
        for edge in edges:
            ecnt += 1
            faces_found = 0
            fcnt = -1
            for face in faces:
                fcnt += 1
                for face_ecnt in np.arange(len(face)):
                    face_edge = [face[face_ecnt],face[(face_ecnt+1)%len(face)]]
                    if np.array_equal(np.sort(face_edge), edge):
                        neighbor_faces[ecnt, faces_found] = fcnt
                        faces_found += 1
                        break
                if faces_found == 2:
                    break
            assert faces_found == 2, 'neighbor faces not found'

    else:
        # -----------------------------
        # Find all Faces (here face is uniquely defined by all its vertices). And Append them to Edges.
        edges, faces, neighbor_faces = find_planar_faces(planar_embedding, 0)
    return edges, faces, neighbor_faces



def find_edges_and_faces_of_toroid(Habs, Hsign, D):    # supports planar and toroidal graphs
    find_edges_and_faces_of_toroid_t0 = time.time()
    L = len(Habs)
    n_edges_total = int(Habs.sum()//2)
    neighbor_faces_mat = np.zeros((L, L), dtype=np.uint16)  # id numbers (1,2,...) of two neighbor faces adjacent to this edge
    faces_toroid = []
    edges = np.zeros((n_edges_total,2), dtype=np.uint16)
    ecnt = 0
    for v0 in range(L):
        for v1 in range(v0+1,L):
            if Habs[v0, v1]!=0:
                edges[ecnt, :] = [v0, v1]
                ecnt += 1
    Habs1 = Habs.copy()
    cut1, cut2 = cut_toroid(Habs1, 0, D)
    H1 = Habs1.astype(np.int8) * (1-2*(Hsign.astype(np.int8)))
    G1 = nx.from_numpy_array(H1, create_using=nx.Graph)
    check_planarity_t0 = time.time()
    is_planar, planar_embedding = nx.check_planarity(G1)
    print('nx.check_planarity time:', int(time.time() - check_planarity_t0), 'sec.')
    assert is_planar

    skip_faces = [cut1, cut2]
    find_planar_faces(planar_embedding, neighbor_faces_mat, faces_toroid, skip_faces)

    Habs1 = Habs.copy()
    cut3, cut4 = cut_toroid(Habs1, L//2, D)
    H1 = Habs1.astype(np.int8) * (1-2*(Hsign.astype(np.int8)))
    G1 = nx.from_numpy_array(H1, create_using=nx.DiGraph)
    is_planar, planar_embedding = nx.check_planarity(G1)
    assert is_planar

    skip_faces = [cut3, cut4]
    find_planar_faces(planar_embedding, neighbor_faces_mat, faces_toroid, skip_faces)

    neighbor_faces_mat -= 1  # (1,2,3,...) -> (0,1,2,...)
    neighbor_faces = get_neighbor_faces_array(edges, neighbor_faces_mat)

    print('find_edges_and_faces_of_toroid() finished, total time:',  int(time.time() - find_edges_and_faces_of_toroid_t0), 'sec.')

    return edges, np.array(faces_toroid), neighbor_faces


def count_pos_edges_in_faces(faces, H):
    pos_edge_cnt = np.zeros((1, len(faces)), dtype=int)
    fcnt = -1
    for f0 in faces:
        fcnt += 1
        for ecnt0 in range(len(f0)):
            ecnt1 = (ecnt0+1)%len(f0)
            e0 = [f0[ecnt0],f0[ecnt1]]
            if H[e0[0],e0[1]]>0:
                pos_edge_cnt[0, fcnt] += 1
    return pos_edge_cnt


def get_maxcut_change_from_maxcut_grps(H, grps):
    N = len(grps)
    # setting to the original matrix,
    # for each edge, is
    M = np.zeros_like(H)
    for idx1 in range(N):
        for idx2 in range(N):
            if (H[idx1, idx2]>0 and grps[idx1] == grps[idx2]) or (H[idx1, idx2]<0 and grps[idx1] != grps[idx2]):
                M[idx1, idx2] = 1
    return M
