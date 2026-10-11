# Liers-Pardella (LP) Max-Cut functions, adapted from planar to toroidal graphs
#
# F. Liers and G. Pardella, "A simple max-cut algorithm for planar graphs" (2008):
# every face of the embedding becomes a K4 with zero-weight edges (faces with more
# than 4 edges are first split into a path of K4's), every primal edge joins the
# K4 corners of its two faces. A maximum-weight perfect matching of this expanded
# dual selects the heaviest edge set that has an even number of edges on each face.
#
# On a planar graph such an edge set is always a cut, so the matching is the exact
# Max-Cut. On a toroid it also has to cross the two non-contractible directions an
# even number of times, which the matching does not see. Therefore here:
#   1. the matching on the toroidal dual gives an upper bound of Max-Cut;
#   2. if the matched edges already form a cut, it is the optimal one;
#   3. otherwise the cut is repaired and improved iteratively: the edges of two
#      non-contractible cycles (one per direction) are fixed to their state in the
#      current cut and the matching is solved again for all the other edges. With
#      these edges fixed every matching is a cut, and it is the best cut among all
#      cuts that agree with the current one on the two cycles, so the value never
#      decreases. The cycles are moved to a random place at every iteration.
# The result is exact when it reaches the bound, otherwise it is a heuristic.
import time
from collections import deque
import numpy as np
import networkx as nx


def build_expanded_dual(edges_array, faces_array):
    # Returns: number of expanded dual nodes, zero-weight (K4 and split-path) edges,
    #          two K4 corners per primal edge, two faces per primal edge
    edge_id = {(int(e[0]), int(e[1])): cnt for cnt, e in enumerate(edges_array)}
    corners = [[] for _ in edges_array]
    edge_faces = [[] for _ in edges_array]
    zero_edges = []
    n_cities = 0
    for fcnt, face in enumerate(faces_array):
        d = len(face)
        face_edges = []
        for vcnt in range(d):
            v0, v1 = int(face[vcnt]), int(face[(vcnt + 1) % d])
            face_edges.append(edge_id[(min(v0, v1), max(v0, v1))])
        # a face with d > 4 edges is split into floor((d-1)/2) nodes of degree <= 4
        k = (d - 1) // 2 if d > 4 else 1
        slot = 0
        prev_out = None
        for i in range(k):
            base = 4 * n_cities
            n_cities += 1
            for a in range(4):
                for b in range(a + 1, 4):
                    zero_edges.append((base + a, base + b))
            corner = base
            if prev_out is not None:
                zero_edges.append((prev_out, corner))
                corner += 1
            if k == 1:
                take = d
            elif i == 0:
                take = 3
            elif i < k - 1:
                take = 2
            else:
                take = d - slot
            for ecnt in face_edges[slot:slot + take]:
                corners[ecnt].append(corner)
                edge_faces[ecnt].append(fcnt)
                corner += 1
            slot += take
            prev_out = corner if i < k - 1 else None
    for ecnt in range(len(edges_array)):
        assert len(corners[ecnt]) == 2, 'each edge must be on the boundary of 2 faces'
        assert edge_faces[ecnt][0] != edge_faces[ecnt][1], 'each edge must have 2 different faces'
    return 4 * n_cities, zero_edges, np.array(corners), np.array(edge_faces)


def match_expanded_dual(n_exp, zero_edges, corners, weights, fixed_edges=None, fixed_state=None):
    # Maximum-weight perfect matching of the expanded dual. Returns bool mask of selected primal edges.
    # fixed_edges/fixed_state: primal edges that must (1) or must not (0) be selected. A selected one
    # takes its two corners out of the matching problem, a not selected one is just left out.
    n_edges = len(corners)
    free = np.ones(n_edges, dtype=bool)
    removed = np.zeros(n_exp, dtype=bool)
    mask = np.zeros(n_edges, dtype=bool)
    if fixed_edges is not None:
        free[fixed_edges] = False
        mask[fixed_edges] = fixed_state
        removed[corners[fixed_edges][np.asarray(fixed_state, dtype=bool)].ravel()] = True
    # all perfect matchings have the same size, so the shift does not change the optimal one
    shift = max(0, -int(weights.min()))
    GE = nx.Graph()
    GE.add_nodes_from(np.flatnonzero(~removed).tolist())
    GE.add_weighted_edges_from((a, b, shift) for a, b in zero_edges if not (removed[a] or removed[b]))
    for ecnt in np.flatnonzero(free):
        GE.add_edge(int(corners[ecnt, 0]), int(corners[ecnt, 1]), weight=int(weights[ecnt]) + shift, primal=int(ecnt))
    matching = nx.algorithms.max_weight_matching(GE, maxcardinality=True, weight='weight')
    assert 2 * len(matching) == GE.number_of_nodes(), 'no perfect matching in the expanded dual'
    for a, b in matching:
        ecnt = GE[a][b].get('primal')
        if ecnt is not None:
            mask[ecnt] = True
    return mask


def paint_nodes_bfs(adj, mask, root=0):
    # Node groups from a set of cut edges. Also returns False if the edges are not a cut.
    N = len(adj)
    grps = -np.ones(N, dtype=int)
    is_cut = True
    for v_start in [root] + list(range(N)):
        if grps[v_start] >= 0:
            continue
        grps[v_start] = 0
        queue = deque([v_start])
        while queue:
            v0 = queue.popleft()
            for v1, ecnt in adj[v0]:
                g = grps[v0] ^ int(mask[ecnt])
                if grps[v1] < 0:
                    grps[v1] = g
                    queue.append(v1)
                elif grps[v1] != g:
                    is_cut = False
    return grps, is_cut


def find_generator_cycles(adj, edges_array, edge_faces, n_faces, root):
    # Edges of the shortest non-contractible cycles through root, one per direction (2 for a toroid, none
    # for a planar graph): BFS tree of the graph + spanning tree of the dual over the remaining edges,
    # taking the edges with the longest cycles first; every edge that is in none of the trees closes a cycle.
    N = len(adj)
    dist = -np.ones(N, dtype=int)
    parent_edge = -np.ones(N, dtype=int)
    parent = -np.ones(N, dtype=int)
    dist[root] = 0
    queue = deque([root])
    while queue:
        v0 = queue.popleft()
        for v1, ecnt in adj[v0]:
            if dist[v1] < 0:
                dist[v1] = dist[v0] + 1
                parent[v1] = v0
                parent_edge[v1] = ecnt
                queue.append(v1)
    in_tree = np.zeros(len(edges_array), dtype=bool)
    in_tree[parent_edge[parent_edge >= 0]] = True
    cycle_len = dist[edges_array[:, 0]] + dist[edges_array[:, 1]] + 1
    uf = list(range(n_faces))

    def find(f):
        while uf[f] != f:
            uf[f] = uf[uf[f]]
            f = uf[f]
        return f

    cycles_edges = set()
    for ecnt in np.argsort(-cycle_len, kind='stable'):
        if in_tree[ecnt]:
            continue
        f0, f1 = find(edge_faces[ecnt, 0]), find(edge_faces[ecnt, 1])
        if f0 != f1:
            uf[f0] = f1
        else:
            cycle = {int(ecnt)}
            for v in edges_array[ecnt]:
                v = int(v)
                while v != root:
                    cycle ^= {int(parent_edge[v])}   # the common part of the two tree paths drops out
                    v = parent[v]
            cycles_edges |= cycle
    return np.array(sorted(cycles_edges), dtype=int)


def get_toroid_max_cut_lp(H, edges_array, faces_array, max_iters=50, max_stall=10, seed=None):
    # Returns: cut value, upper bound, node groups, cut edges mask, info
    # The cut is optimal if cut value == upper bound.
    t0 = time.time()
    rng = np.random.default_rng(seed)
    N = len(H)
    edges_array = np.asarray(edges_array).astype(int)
    weights = H[edges_array[:, 0], edges_array[:, 1]]
    assert np.all(weights == np.round(weights)), 'integer weights are expected'
    weights = weights.astype(int)
    adj = [[] for _ in range(N)]
    for ecnt, (v0, v1) in enumerate(edges_array):
        adj[v0].append((int(v1), ecnt))
        adj[v1].append((int(v0), ecnt))

    n_exp, zero_edges, corners, edge_faces = build_expanded_dual(edges_array, faces_array)
    mask = match_expanded_dual(n_exp, zero_edges, corners, weights)
    mc_bound = int(weights[mask].sum())
    grps, is_cut = paint_nodes_bfs(adj, mask)
    print(f' LP: MCBound={mc_bound}, matched edges are {"" if is_cut else "not "}a cut, time={time.time() - t0:.1f} sec.')

    def cut_of(g):
        m = g[edges_array[:, 0]] != g[edges_array[:, 1]]
        return m, int(weights[m].sum())

    iters = 0
    mask, mc_val = cut_of(grps)
    if not is_cut:
        stall = 0
        root = 0     # same BFS root as the painting above: the first fix changes the matching only on the 2 closing edges
        while iters < max_iters and stall < max_stall and mc_val < mc_bound:
            iters += 1
            fixed_edges = find_generator_cycles(adj, edges_array, edge_faces, len(faces_array), root)
            new_mask = match_expanded_dual(n_exp, zero_edges, corners, weights, fixed_edges, mask[fixed_edges])
            new_grps, new_is_cut = paint_nodes_bfs(adj, new_mask)
            assert new_is_cut
            new_val = int(weights[new_mask].sum())
            assert new_val >= mc_val
            stall = 0 if new_val > mc_val else stall + 1
            print(f' LP: iter {iters}, cycles at node {root}, #fixed edges={len(fixed_edges)}, cut={new_val}, time={time.time() - t0:.1f} sec.')
            # an equally good cut is also taken, to move over plateaus
            mask, grps, mc_val = new_mask, new_grps, new_val
            root = int(rng.integers(N))
    info = {'iters': iters, 'optimal': mc_val == mc_bound, 'time': time.time() - t0}
    return mc_val, mc_bound, grps, mask, info


def get_toroid_max_cut_bound_lp(H, edges_array, faces_array):
    # Upper bound of Max-Cut only (the first matching of get_toroid_max_cut_lp), does not depend on MinPathThreshhold.
    # Also returns True if the matched edges are a cut, i.e. the bound is the Max-Cut value.
    edges_array = np.asarray(edges_array).astype(int)
    weights = H[edges_array[:, 0], edges_array[:, 1]]
    assert np.all(weights == np.round(weights)), 'integer weights are expected'
    weights = weights.astype(int)
    adj = [[] for _ in range(len(H))]
    for ecnt, (v0, v1) in enumerate(edges_array):
        adj[v0].append((int(v1), ecnt))
        adj[v1].append((int(v0), ecnt))
    n_exp, zero_edges, corners, edge_faces = build_expanded_dual(edges_array, faces_array)
    mask = match_expanded_dual(n_exp, zero_edges, corners, weights)
    grps, is_cut = paint_nodes_bfs(adj, mask)
    return int(weights[mask].sum()), is_cut
