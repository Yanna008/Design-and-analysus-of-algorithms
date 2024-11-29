import networkx as nx
from student_utils import *

import pulp
pulp.LpSolverDefault.msg = False
import itertools

def ptp_solver(G:nx.DiGraph, H:list, alpha:float):
    """
    PTP sovler.
    Input:
        G: a NetworkX graph representing the city.\
        This directed graph is equivalent to an undirected one by construction.
        H: a list of home nodes that you must vist.
        alpha: the coefficient of calculating cost.
    Output:
        tour: a list of nodes traversed by your car.
        pick_up_locs_dict: a dictionary of (pick-up-locations, friends-picked-up) pairs\
        where friends-picked-up is a list/tuple containing friends who get picked up at\
        that sepcific pick-up location. Friends are represented by their home nodes.

    All nodes are reprented as integers.
    
    The tour must begin and end at node 0.
    It can only go through edges that exist in the graph..
    Pick-up locations must be in the tour.
    Everyone should get picked up exactly once
    """

    dis = dict(nx.all_pairs_dijkstra_path_length(G, weight='weight'))

    # decision variables:
    X = pulp.LpVariable.dicts('X', G.edges, cat='Binary')
    Y = pulp.LpVariable.dicts('Y', list(itertools.product(H, G.nodes)), cat='Binary')
    R = pulp.LpVariable.dicts('Residual', G.edges, 0, cat='Integer')

    problem = pulp.LpProblem('PartyTogetherProblem', pulp.LpMinimize)

    # minimize:
    problem += pulp.lpSum([X[(u, v)] * G[u][v]['weight'] for u, v in G.edges]) * alpha + pulp.lpSum(Y[(i, j)] * dis[i][j] for i in H for j in G.nodes())

    # subject to:

    for i in G.nodes:
        problem += pulp.lpSum(X[(j, i)]for j in G.predecessors(i)) == pulp.lpSum(X[(i, j)] for j in G.successors(i))
    
    for i in H:
        problem += pulp.lpSum(Y[(i, j)] for j in G.nodes) == 1

    for i in G.nodes:
        if i == 0:
            continue
        problem += pulp.lpSum(R[(i, u)] for u in G.successors(i)) - pulp.lpSum(R[(u, i)] for u in G.predecessors(i)) == pulp.lpSum(Y[(j, i)] for j in H)
    
    # Eliminate subtours
    lenH = len(H)
    problem += pulp.lpSum(R[(0, u)] for u in G.successors(0)) == pulp.lpSum(Y[(i, 0)] for i in H)
    problem += pulp.lpSum(R[(u, 0)] for u in G.predecessors(0)) == lenH

    for u, v in G.edges():
        problem += X[(u, v)] * lenH >= R[(u, v)]
    
    # solve
    problem.solve(pulp.PULP_CBC_CMD(options=['sec=1400'], msg=False))
    # result = pulp.value(problem.objective)
    # print(result)

    # extract tour & locs
    tour = []
    pick_up_locs_dict = {}

    tour = [0]

    nNodes = G.number_of_nodes()
    Xij = [[0] * nNodes for _ in range(nNodes)]
    for u, v in G.edges:
        Xij[u][v] = int(X[(u, v)].varValue)
    
    while True:
        now = tour[-1]
        new = None

        for v in G.successors(now):
            if Xij[now][v] >= 1:
                Xij[now][v] -= 1
                new = v
                break

        if new is None:
            break
        tour.append(new)
    
    for i in H:
        pi = None
        for j in G.nodes:
            if Y[(i, j)].varValue > 0.99:
                pi = j
                if j not in tour:
                    pi = min(((p, dis[i][p]) for p in tour), key=lambda x: x[1])[0]
                    # print("new pi:", pi)
                break

        if pi in pick_up_locs_dict:   
            pick_up_locs_dict[pi].append(i)  
        else:   
            pick_up_locs_dict[pi] = [i]  

    return tour, pick_up_locs_dict


if __name__ == "__main__":
    import sys
    G, H, alpha = input_file_to_instance(sys.argv[1])
    tour, pick_up_locs_dict = ptp_solver(G, H, alpha)
    print(' '.join(map(str, tour)))
    print(len(pick_up_locs_dict))
    for k, v in pick_up_locs_dict.items():
        print(f"{k} {' '.join((map(str, v)))}")
