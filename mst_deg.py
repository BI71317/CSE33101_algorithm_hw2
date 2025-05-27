import random
import networkx as nx
from collections import Counter

def generate_random_complete_graph(n, weight_range=(1, 100)):
    G = nx.complete_graph(n)
    for u, v in G.edges():
        G[u][v]['weight'] = random.randint(*weight_range)
    return G

def analyze_mst_degree_parity(G):
    mst = nx.minimum_spanning_tree(G, algorithm='kruskal')
    degrees = dict(mst.degree())
    parity_count = Counter(['odd' if deg % 2 == 1 else 'even' for deg in degrees.values()])
    return parity_count['odd'], parity_count['even']

def run_experiments(n=100, trials=100):
    total_odd = 0
    total_even = 0

    for _ in range(trials):
        G = generate_random_complete_graph(n)
        odd, even = analyze_mst_degree_parity(G)
        total_odd += odd
        total_even += even

    total = total_odd + total_even
    odd_ratio = total_odd / total
    even_ratio = total_even / total

    print(f"After {trials} trials with {n}-node MSTs:")
    print(f"Avg Odd-degree nodes: {total_odd / trials:.2f}")
    print(f"Avg Even-degree nodes: {total_even / trials:.2f}")
    print(f"Odd-degree ratio: {odd_ratio:.2%}")
    print(f"Even-degree ratio: {even_ratio:.2%}")

run_experiments(n=100, trials=100)