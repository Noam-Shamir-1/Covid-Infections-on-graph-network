from typing import List, Set, Dict
import matplotlib.pyplot as plt
import pandas as pd
import numpy as np
import networkx
import copy


def LTM(graph: networkx.Graph, patients_0: List, iterations: int) -> Set:
    total_infected = set(patients_0)
    #initialize
    for i in graph.nodes:
        graph.nodes[i]["concerned"] = 0

    #iteration
    for i in range(iterations):
        temp_infected=set()
        for t in graph.nodes:
            if t in total_infected:
                continue
            sum_weight = 0
            for k in graph.neighbors(t):
                    if k in total_infected:
                        sum_weight += graph[t][k]["weight"]
            if CONTAGION * sum_weight >= 1 + graph.nodes[t]["concerned"]:
                temp_infected.add(t)

        for j in graph.nodes:
            counter = 0
            for k in graph.neighbors(j):
                if k in total_infected:
                    counter = counter + 1
            if len(graph[j].keys()) > 0:
                graph.nodes[j]["concerned"] = counter/len(graph[j].keys())

        for infect in temp_infected:
            total_infected.add(infect)

    return total_infected

def ICM(graph: networkx.Graph, patients_0: List, iterations: int) -> [Set, Set]:
    total_deceaded = set()
    total_infected = set()
    s = set(graph.nodes) - set(patients_0)
    NI = set()

    #initialize
    for i in graph.nodes:
        graph.nodes[i]["concern"] = 0

    #first - kill iteration 0
    for p0 in patients_0:
        rand = np.random.rand()
        if rand <= LETHALITY:
            total_deceaded.add(p0)
        else:
            NI.add(p0)

    for i in range(iterations):
        temp_deceaded = set()
        temp_infected = set()
        for t in s:
            for k in set.intersection(set(graph.neighbors(t)), NI):
                rand1 = np.random.rand()
                if rand1 <= min(1., CONTAGION * graph[k][t]["weight"] * (1 - graph.nodes[t]["concern"])):# infected
                    rand2 = np.random.rand()
                    if rand2 <= LETHALITY:  # removed
                        temp_deceaded.add(t)
                    else:
                        temp_infected.add(t)
                    break

            # compute concerned for all healthy nodes
            infected_neighbor = 0
            death_neighbor = 0
            count_neighbor = 0
            for k in graph.neighbors(t):
                count_neighbor += 1
                if k in total_infected or k in NI:
                    infected_neighbor += 1
                elif k in total_deceaded:
                    death_neighbor += 1
            long_expression = 1.0
            if count_neighbor > 0:
                long_expression = (infected_neighbor + (3 * death_neighbor)) / count_neighbor
            graph.nodes[t]["concern"] = min(1.0, long_expression)

        #update: total_infected, new_infected_ total_deceaded
        total_infected = set.union(total_infected, NI)
        total_deceaded = set.union(total_deceaded,temp_deceaded)
        NI = temp_infected
        s = s - temp_infected - temp_deceaded

    return total_infected.union(NI), total_deceaded


def plot_degree_histogram(histogram: Dict):
    plt.bar(list(histogram.keys()), histogram.values(), color='g')
    plt.show()


def calc_degree_histogram(graph: networkx.Graph) -> Dict:
    """
    Example:
    if histogram[1] = 10 -> 10 nodes have only 1 friend
    """
    G = graph
    histogram = {}
    for i in networkx.nodes(G):
        if len(list((networkx.neighbors(G,i)))) in histogram:
            histogram[len(list((networkx.neighbors(G, i))))] +=1
        else:
            histogram[len(list((networkx.neighbors(G, i))))] = 1
    histogram = dict(sorted(histogram.items()))
    return histogram


def build_graph(filename: str) -> networkx.Graph:
    G = networkx.Graph()
    # test_list = []
    if filename!='PartB-C.csv':
        df = pd.read_csv(filename, usecols=['from', 'to'])
        for i in df.values:
            G.add_edge(i[0], i[1])
    else:
        df = pd.read_csv(filename, usecols=['from', 'to','w'])
        for i in df.values:
            G.add_weighted_edges_from([(i[0], i[1], i[2])])
            # test_list.append(i[2])
    # print("mean edge weight is: ", np.mean(test_list), " and median edge weight is: ", np.median(test_list))
    return G


def clustering_coefficient(graph: networkx.Graph) -> float:
    G=graph
    triangles_dic = networkx.triangles(G)
    count_triangles = 0
    for t in triangles_dic.values():
        count_triangles += t

    connected_triplets = 0
    for node in G.nodes():
        degree_node = len(list((networkx.neighbors(G, node))))
        if (degree_node>1):
            connected_triplets += np.math.factorial(degree_node) / (2 * np.math.factorial(degree_node - 2))
    cc = count_triangles/connected_triplets
    return cc


def compute_lethality_effect(graph: networkx.Graph, t: int) -> [Dict, Dict]:
    global LETHALITY
    mean_deaths = {}
    mean_infected = {}

    for l in (.05, .15, .3, .5, .7):
        LETHALITY = l
        infected_list = []
        deceaded_list = []
        for iteration in range(30):
            G = copy.deepcopy(graph)
            patients_0 = np.random.choice(list(G.nodes), size=50, replace=False, p=None)
            infected, dead = ICM(G, patients_0, t)
            infected_list.append(len(infected))
            deceaded_list.append(len(dead))
        mean_deaths[l] = np.mean(infected_list)
        mean_infected[l] = np.mean(deceaded_list)

    return mean_deaths, mean_infected


def plot_lethality_effect(mean_deaths: Dict, mean_infected: Dict):
    plt.plot(mean_deaths.keys(), mean_deaths.values(), label="death", color="red")
    plt.plot(mean_infected.keys(), mean_infected.values(), label="infected", color="blue")
    plt.xlabel("lethality value")
    plt.ylabel("mean number of peoples")
    plt.legend()
    plt.show()


def choose_who_to_vaccinate(graph: networkx.Graph) -> List:
#test 1
    node2degree = dict(graph.degree(weight="weight"))
    sorted_nodes = sorted(node2degree.items(), key=lambda item: item[1], reverse=True)[:50]
    people_to_vaccinate1 = [node[0] for node in sorted_nodes]

#test 2
    node2degree = dict(graph.degree)
    sorted_nodes = sorted(node2degree.items(), key=lambda item: item[1], reverse=True)[:50]
    people_to_vaccinate2 = [node[0] for node in sorted_nodes]


    patients0 = list(pd.read_csv("patients0.csv", header=None)[0])

    g1 = copy.deepcopy(graph)
    g2 = copy.deepcopy(graph)
    for i in people_to_vaccinate1:
        g1.remove_node(i)
    len_g1 = len(ICM(g1, patients0, 6)[0])

    for i in people_to_vaccinate2:
        g2.remove_node(i)
    len_g2 = len(ICM(g2, patients0, 6)[0])
    # print("len1 is ", len_g1, " len 2 is",len_g2)
    if len_g1 < len_g2:
        return people_to_vaccinate1
    else:
        return people_to_vaccinate2

def choose_who_to_vaccinate_example(graph: networkx.Graph) -> List:
    """
    The following heuristic for Part C is simply taking the top 50 friendly people;
     that is, it returns the top 50 nodes in the graph with the highest degree.
    """
    node2degree = dict(graph.degree)
    sorted_nodes = sorted(node2degree.items(), key=lambda item: item[1], reverse=True)[:50]
    people_to_vaccinate = [node[0] for node in sorted_nodes]
    return people_to_vaccinate


"Global Hyper-parameters"
CONTAGION = 0.8
LETHALITY = .2

if __name__ == "__main__":
    filename1 = "PartA1.csv"
    G1 = build_graph(filename1)
    dict_histogram1 = calc_degree_histogram(G1)
    plot_degree_histogram(dict_histogram1)

    filename2 = "PartA2.csv"
    G2 = build_graph(filename2)
    dict_histogram2 = calc_degree_histogram(G2)
    plot_degree_histogram(dict_histogram2)


    filename3 = "PartB-C.csv"
    G3 = build_graph(filename3)
    dict_histogram3 = calc_degree_histogram(G3)
    plot_degree_histogram(dict_histogram3)

    print("clustering coefficient of graph A1 is:", clustering_coefficient(G1))
    print("clustering coefficient of graph A2 is:", clustering_coefficient(G2))


    # dict_infected, dict_deceaded = compute_lethality_effect(G3, 6)
    # plot_lethality_effect(dict_deceaded, dict_infected)


    # patients0 = list(pd.read_csv("patients0.csv", header=None)[0])
    # print(len(LTM(G3, patients0[0:20], 6)))
    # print("normal graph")
    # print(len(ICM(G3, patients0, 6)[0]))

    print(choose_who_to_vaccinate(G3))
