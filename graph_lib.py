# General graph functions library
import numpy as np
import networkx as nx
import matplotlib.pyplot as plt

def plot_graph(G, pos, title="Graph", non_planar_edges=None):
    """Plot a graph with optional highlighting of non-planar edges."""
    plt.figure(figsize=(8, 8))
    nx.draw(G, pos, with_labels=True, node_color='lightblue', edge_color='gray',
            node_size=500, font_size=10, font_color='black')
    if non_planar_edges:
        # Highlight non-planar edges in red if any
        nx.draw_networkx_edges(G, pos, edgelist=non_planar_edges, edge_color='red', style='dashed')
    plt.title(title)
    plt.show()


def draw_planar_representation(G, planar_embedding):
    """Draw the graph in a planar layout."""
    if planar_embedding is None:
        pos = nx.planar_layout(G)
    else:
        pos = nx.planar_layout(G, planar_embedding)
    plot_graph(G, pos, title="Planar Representation")


def separate_planar_non_planar(G, quiet_en=0):
    """
    Separate a graph into planar and non-planar parts.
    Returns the planar embedding, planar subgraph, and non-planar edges.
    """
    is_planar, planar_embedding = nx.check_planarity(G)

    if is_planar:
        print("The graph is planar.")
        return True, planar_embedding, G, []

    if not quiet_en:
        print("The graph is not planar. Separating planar and non-planar parts...")

    planar_subgraph = nx.DiGraph()

    planar_subgraph.add_nodes_from(G.nodes)

    # Build a maximal planar subgraph
    for edge in G.edges:
        planar_subgraph.add_edge(*edge, weight=G.edges[edge]['weight'])
        if not nx.check_planarity(planar_subgraph)[0]:
            planar_subgraph.remove_edge(*edge)

    # Identify non-planar edges
    non_planar_edges = [edge for edge in G.edges if not planar_subgraph.has_edge(*edge)]

    return False, planar_embedding, planar_subgraph, non_planar_edges


def graph_to_adjacency_matrix(G):
    """Convert a graph to its adjacency matrix."""
    return nx.to_numpy_array(G, weight='weight')


def create_non_planar_matrix(G, non_planar_edges):
    """Create an adjacency matrix for the non-planar edges."""
    non_planar_matrix = np.zeros((len(G.nodes), len(G.nodes)))
    for edge in non_planar_edges:
        i, j = edge
        non_planar_matrix[i, j] = G.edges[edge]['weight']
    return non_planar_matrix


def plot_graph_from_adjacency_matrix(matrix, plot_en=0):
    """Create and plot graphs from a given weighted adjacency matrix."""
    # Ensure the matrix is square
    if matrix.shape[0] != matrix.shape[1]:
        raise ValueError("The adjacency matrix must be square.")

    # Create graph from adjacency matrix
    G = nx.from_numpy_array(matrix, create_using=nx.DiGraph)

    # Check and separate planar and non-planar parts
    is_planar, planar_embedding, planar_subgraph, non_planar_edges = separate_planar_non_planar(G)
    print(is_planar, planar_embedding, planar_subgraph, non_planar_edges)
    if is_planar:
        if plot_en:
            draw_planar_representation(G, None)
    else:
        # Plot the planar subgraph
        pos_planar = nx.spring_layout(planar_subgraph)
        if plot_en:
            plot_graph(planar_subgraph, pos_planar, title="Planar Subgraph")

        # Highlight non-planar edges
        pos_full = nx.spring_layout(G)
        if plot_en:
            plot_graph(G, pos_full, title="Non-Planar Edges Highlighted in Red", non_planar_edges=non_planar_edges)

    # Convert planar subgraph to a matrix
    planar_matrix = graph_to_adjacency_matrix(planar_subgraph)
    print("Planar Subgraph Matrix:")
    print(planar_matrix)

    # Convert non-planar edges to a matrix
    non_planar_matrix = create_non_planar_matrix(G, non_planar_edges)
    print("\nNon-Planar Edges Matrix:")
    print(non_planar_matrix)
    return planar_matrix, non_planar_matrix


def get_cut(edges, b, M=np.atleast_2d([-1])):
    cut_val = 0
    for edge in edges:
        if b[edge[0]] != b[edge[1]]:
            if len(M)==1:
                cut_val += 1
            else:
                cut_val += M[edge[0],edge[1]]
    return int(cut_val)


def get_cut_simple(grp, wei):
    cut_val = 0
    for cnt0 in range(len(wei)):
        for cnt1 in range(cnt0+1,len(wei)):
            if grp[cnt0] != grp[cnt1]:
                cut_val += wei[cnt0,cnt1]
    return int(cut_val)



def get_cut_adjmat(b, adj_mat):
    cut_val = 0
    for ii in range(len(adj_mat)):
        for jj in range(ii+1,len(adj_mat)):
            if b[ii] != b[jj]:
                cut_val += adj_mat[ii,jj]
    return cut_val


def max_cut_exhaustive_search(N, edges, adj_mat=np.atleast_2d([-1])):
    max_cut_optimal_val = 0
    max_cut_optimal_grp = []
    for cnt in range(2**N):
        n = cnt
        b = np.zeros(N, dtype=int)
        for k in range(N):
            b[k] = n % 2
            # print(N, cnt, n, k, b)
            n = n >> 1
        cut_val = get_cut(edges, b, adj_mat)
        # print(cnt, b, cut_val)
        if cut_val > max_cut_optimal_val:
            max_cut_optimal_val = cut_val
            max_cut_optimal_grp = b
    return max_cut_optimal_grp, max_cut_optimal_val

def get_clock_wise_dir_edges(edges_in):
    nodes_num = edges_in.max() + 1
    edges_dir = edges_in.copy()
    for ecnt in range(len(edges_dir)):
        n0, n1 = edges_dir[ecnt, 0].astype(np.int16), edges_dir[ecnt, 1].astype(np.int16)
        d0, d1 = (n1 - n0) % nodes_num, (n0 - n1) % nodes_num
        if d0 > d1:
            edges_dir[ecnt, :] = edges_dir[ecnt, ::-1]
    return edges_dir

