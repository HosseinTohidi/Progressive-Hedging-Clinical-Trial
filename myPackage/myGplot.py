# -*- coding: utf-8 -*-
"""
Created on Tue Dec  3 18:16:48 2019

@author: atohidi
"""

import networkx as nx
import matplotlib.pyplot as plt
def myGplot(G, node_label = False, edge_label = False):
    pos = nx.get_node_attributes(G, 'pos')
    node_colors = [G.nodes[u]['nColor'] for u in G.nodes]    
    edge_colors = [G[u][v]['edge_color'] for u,v in G.edges]
    nx.draw(G, pos, edge_color = edge_colors, node_color = node_colors)
    if edge_label:
        weight = nx.get_edge_attributes(G, 'weight')
        nx.draw_networkx_edge_labels(G, pos,edge_labels = weight, label_pos = 0.7)
    if node_label:
        nodes = nx.get_node_attributes(G, 'name')
        nx.draw_networkx_labels(G, pos, nodes)
    plt.show()