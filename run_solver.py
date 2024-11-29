import os
import networkx as nx
import json
from ptp_solver import ptp_solver

INPUT_FILE_DIRECTORY = 'inputs'

def list_all_files(directory, extension):
    """
    List all files under path
    """
    files = get_files_with_extension(directory, extension)
    for file in files:
        print(file)

def get_files_with_extension(directory, extension):
    """
    Get all files end with specified extension under directory
    """
    files = []
    for file in os.listdir(directory):
        if file.endswith(extension):
            files.append(file)
    return sorted(files)

def input_file_names_to_file_path(user_in_files_str, in_files_all):
    user_in_files = user_in_files_str.split()
    file_paths = []
    message = ''
    for file in user_in_files:
        if file in in_files_all:
            file_paths.append(os.path.join(os.path.join(os.getcwd(), INPUT_FILE_DIRECTORY), file))
        else:
            file_con = file + '.in'
            if file_con in in_files_all:
                file_paths.append(os.path.join(os.path.join(os.getcwd(), INPUT_FILE_DIRECTORY), file_con))
            else:
                message += f'{file} '
    if message:
        message = 'Input ' + message + 'Not Exist'
    return file_paths, message

def read_file(file):
    """
    Read all lines in file 
    Store the data in a list of lines
    where each line is splited to list of words   
    """
    with open(file, 'r') as f:
        lines = f.readlines()
    return [line.strip().split() for line in lines]

def write_to_file(file, data, mode='w'):
    """
    Write data into file
    Default mode: 'w'
    """
    with open(file, mode) as f:
        f.write(data)

def data_parser(input_data):
    """
    Parsing input data
    """
    alpha = float(input_data[0][0])
    number_of_nodes = int(input_data[1][0])
    number_of_homes = int(input_data[1][1])
    home_node_list = [int(node) for node in input_data[2]]
    multiline_adjlist = input_data[3:]
    edge_list = []
    i = 0
    while i < len(multiline_adjlist):
        line = multiline_adjlist[i]
        u, d = int(line[0]), int(line[1])
        i += 1
        for _ in range(d):
            line = multiline_adjlist[i]
            v, w = int(line[0]), float(line[1])
            i += 1
            edge_list.append((u, v, w))
    return alpha, number_of_nodes, number_of_homes, home_node_list, edge_list

def weighted_edge_list_to_graph(edge_list):
    """
    Create a graph from a weighted edge list
    For technical reasons, the graph is made directed here.
    But it is equivalent to an undirected one by construction.
    """
    G = nx.DiGraph()
    G.add_weighted_edges_from(edge_list)
    return G

def input_file_to_instance(file):
    """
    Create an instance of PTP problem from a specific file
    Input: file: path of the file
    Output: 
        G: the graph, 
        H: a list of home nodes
        alpha: the coeeficient alpha
        
    I also store H and graph in the graph G.
    You can get access to them via G.graph['alpha'], G.graph['H']
    """
    input_data = read_file(file)
    alpha, _, _, H, edge_list = data_parser(input_data)
    G = weighted_edge_list_to_graph(edge_list)
    G.graph['H'] = H
    G.graph['alpha'] = alpha
    return G, H, alpha

def is_metric(G):
    """
    Check whether a given graph G is metric or not,
    i.e., whether triangle inequality holds.
    """
    d = nx.floyd_warshall(G)
    for u, v, data in G.edges(data=True):
        if abs(d[u][v] - data['weight']) >= 0.0000001:
            return False
    return True

def is_connected(G):
    """
    Check whether a graph G is connected or not
    """
    return nx.is_connected(nx.to_undirected(G))

def is_valid_input(file):
    """
    Check if an input is valid or not
    Input: file: path to the input file
    Output:
        is_valid: bool, input valid or not
        message: log message
    """
    is_valid = True
    message = ''

    input_data = read_file(file)
    file_name = os.path.basename(file)
    file_name = file_name.split(".")[0]
    file_name = file_name.split("_")
    max_number_of_nodes = int(file_name[0])
    max_number_of_friends = int(max_number_of_nodes/2)
    expected_alpha = float(file_name[1])/10
    try:
        alpha, number_of_nodes, number_of_friends, H, edge_list = data_parser(input_data)
    except:
        return False, "Cannot parse data"

    if abs(alpha - expected_alpha) >= 0.001:
        is_valid = False
        message += 'alpha is not as required\n'

    if number_of_nodes > max_number_of_nodes:
        is_valid = False
        message += 'maximum number of nodes exceeded\n'

    if number_of_friends > max_number_of_friends:
        is_valid = False
        message += 'maximum number of friends exceeded\n'

    if len(H) != number_of_friends:
        is_valid = False
        message += 'number of friends not correct\n'

    if len(set(H)) != len(H):
        is_valid = False
        message += 'friend homes not distinct\n'

    try:
        G = weighted_edge_list_to_graph(edge_list)
    except:
        message += 'Cannot create graph'
        return False, message

    if G.number_of_nodes() != number_of_nodes:
        is_valid = False
        message += 'number of nodes not correcr\n'

    if not is_connected(G):
        is_valid = False
        message += 'graph is not connected\n'

    if not is_metric(G):
        is_valid = False
        message += 'graph does not have triangle inequality\n'

    nodes = sorted(G.nodes())

    if nodes != list(range(G.number_of_nodes())):
        is_valid = False
        message += 'nodes not indexed from 0 to n-1\n'

    if any([h not in nodes for h in H]):
        is_valid = False
        message += 'some home nodes not in the graph\n'

    for u, v in G.edges():
        if u == v:
            is_valid = False
            message += 'exist edge conneting a node with itself\n'
            break
        if G[u][v]['weight'] != G[v][u]['weight']:
            is_valid = False
            message += 'edge weights not symmetric\n'
            break
    
    return is_valid, message


def analyze_solution(G, H, alpha, tour, pick_up_locs_dict):
    """
    Analyze solutuon for a given instance of the problem
    Input:
        G: the graph G
        H: a list of home nodes
        alpha: the cost coefficient
        tour: the tour of the car
        pick_up_locs_dict: a dictionary of (pick-up location, friends picked up) pair
                        for pthp, this would be empty({})
    Output:
        is_legitimate: bool, whether the solution is legitimate or not
        driving_cost: folat, total cost(unhappniess) of the driving car
        walking_cost: float, total cost of walking friends

    A solution is legitimate iff following coditions hold.
        The indices must be in the graph, i.e., integers from 0 to |V| - 1. 
        The tour must begin and end at node 0. It can only go through edges that exist in the graph. 
        The pick-up locations must be in the tour . Everyone should get picked up.
    Total cost would be the sum of driving_cost and walking_cost
    An infeasible solution would have positive infinity cost

    Example:
        For the example in the description file, a possible solution would be
            tour = [0, 1, 0]
            pick_up_locs_dict = {1: (1, 2 3)}
        The output would thus be, (assume alpha = 2/3)
            True, float(10/3)
    """
    driving_cost = 0
    walking_cost = 0
    # the tour must start and end at node 0
    if not (tour[0] == 0 and tour[-1] == 0):
        print("not cycle")
        return False, float('infinity'), float('infinity') 
    # every road in the tour must exist in the graph
    for i in range(1, len(tour)):
        if not G.has_edge(tour[i-1], tour[i]):
            print(f"edge{tour[i-1], tour[i]} not exist")
            return False, float('infinity'), float('infinity')
        driving_cost += float(alpha * G.get_edge_data(tour[i-1], tour[i])['weight'])
    # every should get picked up exactly once    
    if pick_up_locs_dict:
        all_shortest_path_lengths = nx.floyd_warshall(G)
        friends_get_picked_up = []
        for pick_up_loc in pick_up_locs_dict:
            friends = pick_up_locs_dict[pick_up_loc]
            # pick up locations must be in the tour
            if not pick_up_loc in tour:
                print(f"Pick up location {pick_up_loc} not in the tour")
                return False, float('infinity'), float('infinity')
            for friend in friends:
                walking_cost += all_shortest_path_lengths[pick_up_loc][friend]
                friends_get_picked_up.append(friend)
        friends_get_picked_up = sorted(friends_get_picked_up)
        if len(H) == len(friends_get_picked_up) \
                and all([home == friends_get_picked_up[i]] for i, home in enumerate(sorted(H))):
            return True, driving_cost, walking_cost
        else:
            print("Not all friends picked up")
            return False, float('infinity'), float('infinity')
        # return np.all(np.sort(H) == np.sort(friends_get_picked_up))
    # PTHP solution, no pick-up locaitons. Every node in H should be in tour
    if all([i in tour for i in H]):
        return True, driving_cost, walking_cost
    else:
        print("Not visit every node in H")
        return False, float('infinity'), float('infinity')
    
def run_solver(in_files, ptp_solver):
    count = 0
    costs = {}
    for file in in_files:
        print(f"\nReading file {os.path.basename(file)}...")
        G, H, alpha = input_file_to_instance(file)
        print(f"n = {G.number_of_nodes()}, |H| = {len(H)}, alpha = {alpha}")
        # print('Graph constructed...')
        try:
            ptp_tour, ptp_pick_up_locs_dict = ptp_solver(G, H, alpha)
        # print('Tour generated by PTP solver...')
        # print('Analyzing the solution...')
            ptp_is_legitimate, ptp_driving_cost, ptp_walking_cost = analyze_solution(G, H, alpha, ptp_tour, ptp_pick_up_locs_dict)
            ptp_cost = ptp_driving_cost + ptp_walking_cost
            print(f"Your PTP solution is{' NOT' if not ptp_is_legitimate else ''} legimitate.")
            print(f"Total cost of your PTP solution: {ptp_cost:.6f}")
        # print(f"Total driving cost of your PTP solution: {ptp_driving_cost:.6f}")
        # print(f"Total walking cost of your PTP solution: {ptp_walking_cost:.6f}")
        except:
            print("Error Executing Solver")
            ptp_cost = float("infinity")
        costs[os.path.basename(file)] = ptp_cost
        count += 1
    print(f"\nSucessfully test {count} files")
    return costs


in_file_dir = os.path.join(os.getcwd(), INPUT_FILE_DIRECTORY)
in_files = get_files_with_extension(in_file_dir, extension=".in")
in_files = [os.path.join(in_file_dir, file) for file in in_files]
costs = run_solver(in_files, ptp_solver)
out_file = os.path.join(os.getcwd(), "costs.txt")
with open(out_file, 'w') as file:
    json.dump(costs, file, indent=4)
print(f"Costs written to {out_file}")