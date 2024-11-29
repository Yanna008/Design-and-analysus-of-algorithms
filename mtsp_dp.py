import networkx as nx
import itertools

def mtsp_dp(G):
    """
    TSP solver using dynamic programming.
    Input:
        G: a NetworkX graph representing the city.\
        This directed graph is equivalent to an undirected one by construction.
    Output:
        tour: a list of nodes traversed by your car.

    All nodes are reprented as integers.

    You must solve the problem using dynamic programming.
    
    The tour must begin and end at node 0.
    It can only go through edges that exist in the graph..
    It must visit every node in G exactly once.
    """
    nNodes = G

    nNodes = G.number_of_nodes()
    nSet = 2 ** (nNodes - 1)

    # dp[i][S] 0-> S -> i
    dp = [[float('inf')] * nSet for _ in range(nNodes)]
    for i in range(1, nNodes):
        dp[i][0] = G[0][i]['weight']
    
    for node_set in range(1, nSet):
        for i in range(nNodes):
            
            if i != 0 and ((1 << (i - 1)) & node_set) != 0:
                # S contains i, error
                continue

            for j in range(1, nNodes):

                set_j = 1 << (j - 1)

                if (set_j & node_set) == 0:
                    # S does not contains j, error
                    continue

                # S contains j, S doesn't contain i
                # dp[i][S] = min_j { dp[j][S\j] + w[j][i] }
                dp[i][node_set] = min(dp[i][node_set], dp[j][node_set ^ set_j] + G[j][i]['weight'])

    ans_value = dp[0][nSet - 1]

    rev_path = [0, ]
    rem_state = nSet - 1
    while rem_state != 0:
        now = rev_path[-1]
        for prev in range(1, nNodes):
            set_prev = 1 << (prev - 1)
            if (rem_state & set_prev) != 0 \
                and dp[now][rem_state] == dp[prev][rem_state ^ set_prev] + G[prev][now]['weight']:

                rev_path.append(prev)
                rem_state ^= set_prev
                break
    rev_path.append(0)
    rev_path.reverse()

    return rev_path

# mtsp_dp(G)