import numpy as np
import networkx as nx


def find_inner_cycle_path_of_toroid_simple(H, cycle_start=0):
    ii_connected = np.argwhere(H[cycle_start, :] != 0)
    ii_connected = ii_connected[ii_connected>cycle_start]
    if len(ii_connected) >1:
        cycle_end = int(ii_connected[1])
        G = nx.from_numpy_array(H, create_using=nx.DiGraph)
        in_cycle_path = []
        for vi in range(cycle_start,cycle_end):
            sh_paths = nx.shortest_path(G, source=vi, target=vi+1 , method='dijkstra')
            for v in sh_paths[:-1]:
                in_cycle_path.append(v)
        in_cycle_path.append(cycle_end)
    else:
        in_cycle_path = []
    return in_cycle_path


def find_inner_cycle_of_toroid(H, cycle_start = 0):
    ii_connected = np.argwhere(H[cycle_start, :] != 0)
    ii_connected = ii_connected[ii_connected>cycle_start]
    if len(ii_connected) >1:
        cycle_end = int(ii_connected[1])

        G = nx.from_numpy_array(H, create_using=nx.DiGraph)
        in_cycle = [[cycle_start, cycle_end]]
        for vi in range(cycle_start,cycle_end):
            sh_paths = nx.shortest_path(G, source=vi, target=vi+1 , method='dijkstra')
            for ecnt in range(len(sh_paths)-1):
                current_edge = sh_paths[ecnt:ecnt+2]
                in_cycle.append(current_edge)
        # print(in_cycle)
    else:
        in_cycle = [[-1 ,-1]]
    return in_cycle


def cut_toroid(H, cycle_start, D):
    inner_cycle_path_1 = find_inner_cycle_path_of_toroid_simple(H, cycle_start)
    inner_cycle_path_2 = []
    for v0 in inner_cycle_path_1:
        for v1 in range(v0 + 1, v0 + 2 * D):
            if not v1 in inner_cycle_path_1:
                if H[v0, v1] != 0:
                    H[v0, v1] = 0
                    H[v1, v0] = 0
                    inner_cycle_path_2.append(v1)
    inner_cycle_path_1 = np.sort(inner_cycle_path_1)
    inner_cycle_path_2 = np.sort(inner_cycle_path_2)
    return inner_cycle_path_1, inner_cycle_path_2


def find_inner_cycle_path_of_toroid(H, edges, neighbors, cycle_start = 0):
    ii_connected = np.argwhere(H[cycle_start, :] != 0)
    ii_connected = ii_connected[ii_connected>cycle_start]
    if len(ii_connected) >1:
        cycle_end = int(ii_connected[1])
        G = nx.from_numpy_array(H, create_using=nx.DiGraph)
        in_cycle_path = []
        for vi in range(cycle_start,cycle_end):
            sh_paths = nx.shortest_path(G, source=vi, target=vi+1 , method='dijkstra')
            for ecnt in range(len(sh_paths)-1):
                current_v = sh_paths[ecnt]
                in_cycle_path.append(current_v)
        # print(in_cycle)
        in_cycle_path.append(cycle_end)

        cycle_edges = []
        cycle_dual_edges = []
        for cvcnt in range(len(in_cycle_path)):
            cv1, cv2 = in_cycle_path[cvcnt], in_cycle_path[(cvcnt + 1) % len(in_cycle_path)]
            if cv1 > cv2:
                cv1, cv2 = cv2, cv1
            ecnt = -1
            for edge in edges:
                ecnt += 1
                assert edge[1]>edge[0]
                if edge[0] == cv1 and edge[1] == cv2:
                    prime_edge = [int(x) for x in edge]
                    dual_edge = [int(x) for x in neighbors[ecnt]]
                    cycle_edges.append(prime_edge)
                    cycle_dual_edges.append(dual_edge)
                    break
    else:
        in_cycle_path = []
        cycle_dual_edges = []
        cycle_edges = []
    return in_cycle_path, cycle_dual_edges, cycle_edges



