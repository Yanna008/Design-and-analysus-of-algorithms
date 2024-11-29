import networkx as nx
from mtsp_dp import mtsp_dp
# from student_utils import *

def path_floyd(G):
    nNodes = G.number_of_nodes()
    
    # floyd 
    f = [[float('inf')] * nNodes for _ in range(nNodes)]
    spath_via = [[[] for j in range(nNodes)] for i in range(nNodes)]
    for i in range(nNodes):
        for j in range(nNodes):
            if not G.has_edge(i, j):
                continue
            f[i][j] = G[i][j]['weight']
    
    for k in range(nNodes):
        for i in range(nNodes):
            for j in range(nNodes):
                t = f[i][k] + f[j][k]
                if t >= f[i][j]:
                    continue
                f[i][j] = t
                spath_via[i][j] = spath_via[i][k] + [k] + spath_via[k][j]
    # print("finish path_floyd")
    return f, spath_via

def path_dij(G):
    nNodes = G.number_of_nodes()
    all_pairs_paths = dict(nx.all_pairs_dijkstra_path(G, weight='weight'))
    all_pairs_lengths = dict(nx.all_pairs_dijkstra_path_length(G, weight='weight'))
    f = [[all_pairs_lengths[i][j] for j in range(nNodes)] for i in range(nNodes)]
    spath_via = [[all_pairs_paths[i][j][1:-1] for j in range(nNodes)] for i in range(nNodes)]
    return f, spath_via


def pthp_solver_from_tsp(G, H):
    """
    PTHP sovler via reduction to Euclidean TSP.
    Input:
        G: a NetworkX graph representing the city.\
        This directed graph is equivalent to an undirected one by construction.
        H: a list of home nodes that you must vist.
    Output:
        tour: a list of nodes traversed by your car.

    All nodes are reprented as integers.

    You must solve the question by first transforming a PTHP\
    problem to a TSP problem. After solving TSP via the dynammic\
    programming algorithm introduced in lectures, construct a solution\
    for the original PTHP problem.
    
    The tour must begin and end at node 0.
    It can only go through edges that exist in the graph..
    It must visit every node in H.
    """
    f, spath_via = path_floyd(G)
    
    H_prime = [0] + H
    reduced_graph = nx.Graph()
    reduced_graph.add_nodes_from(list(range(len(H_prime))))
    for i1, h1 in enumerate(H_prime):
        for i2, h2 in enumerate(H_prime):
            if h1 == h2:
                continue
            reduced_graph.add_edge(i1, i2, weight=f[h1][h2])

    # reduction

    tsp_tour = mtsp_dp(reduced_graph)
    # print('?')
    # reduction
    tour = [H_prime[tsp_tour[0]]]
    for u, v in zip(tsp_tour[:-1], tsp_tour[1:]):
        tour += spath_via[H_prime[u]][H_prime[v]]
        tour.append(H_prime[v])

    return tour

if __name__ == "__main__":
    # G, H, alpha = input_file_to_instance('inputs/1.in')
    with open('inputs/1.in') as fp:
        lines = [line.strip().split() for line in fp.readlines()]
    alpha = float(lines[0][0])
    nNodes = int(lines[1][0])
    H = list(map(int, lines[2]))
    G = nx.Graph()
    G.add_nodes_from(list(range(nNodes)))
    nowid = 3
    for i in range(nNodes):
        idx, d = int(lines[nowid][0]), int(lines[nowid][1])
        for did in range(d):
            t, w = int(lines[nowid+did+1][0]), int(lines[nowid+did+1][1])
            G.add_edge(idx, t, weight=w)
        nowid += d + 1
    
    tour = pthp_solver_from_tsp(G, H)
    print(tour)
