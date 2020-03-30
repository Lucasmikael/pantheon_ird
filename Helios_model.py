from matplotlib.ticker import NullFormatter  # useful for `logit` scale
import matplotlib.pyplot as plt
import numpy as np
from matplotlib.backends.backend_tkagg import FigureCanvasTkAgg
import PySimpleGUI as sg
import matplotlib
matplotlib.use('TkAgg')
from pythonis_model import *
import networkx as nx



def drawGraph(genes_names_list,network_as_list, flow) :
    G = nx.MultiDiGraph()


    # G.add_nodes_from(genes_names_list)
    # val_map = {'ARR23': "v"}
    #
    # values = [val_map.get(node, "o") for node in G.nodes()]




    positive_gene, negative_gene,normal_gene = getAutoregulate(network_as_list,genes_names_list)

    G = addNodes(positive_gene,negative_gene,normal_gene, G)
    G = addEdges(network_as_list, G)

    # val_map = {'ARR23': "v"}
    #
    # values = [val_map.get(node, "o") for node in G.nodes()]


    source, target, value_source = getFlow1 (flow)



    activation_table = activationFlowTable (value_source, genes_names_list)
    print(activation_table)

    #ajout des aretes et des noeuds

    for i in range (len(network_as_list)):
        duet = "0"
        gene_source = network_as_list[i][0]
        gene_target = network_as_list[i][2]
        interaction = network_as_list[i][1]
        for k in range (len(network_as_list)) :
             if network_as_list[k][0] == gene_target and network_as_list[k][2] == gene_source :
                duet = "1"
        G.add_edge(gene_source, gene_target, arrowstyle = interaction, duet = duet)


    pos = nx.circular_layout(G)

    for node in G.nodes(data=True):
        print(node[0], ": ",node[1]['forme'])
        if node[1]['forme'] == "1":
            node_shape = "^"
        if node[1]['forme'] == "-1":
            node_shape = "v"
        else :
            node_shape = "o"
        nx.draw_networkx_nodes(G, pos, nodelist= [node[0]], node_size=1500, node_shape = node_shape)



    nx.draw_networkx_labels(G, pos)
    for edge in G.edges(data=True):
        if edge[2]['arrowstyle'] == "-1":
            arrowstyle = "-["
        if edge[2]['arrowstyle'] == "1":
            arrowstyle = "-|>"
        if edge[2]['duet'] == "1":
            form_arrow = 'arc3, rad = 0.2'
        if edge[2]['duet'] == "0":
            form_arrow = 'arc3, rad = 0.0'
        for i in range (len(activation_table)):
            if activation_table[i][0] == edge[0]:
                if activation_table[i][1] == 1:
                    activateedge = "r"
                    print(activateedge)
                else:
                    activateedge = "b"
                    print(activateedge)
        nx.draw_networkx_edges(G, pos, edgelist=[(edge[0],edge[1])], arrowstyle = arrowstyle, node_size=1900, connectionstyle=form_arrow, edge_color= activateedge)
    # plt.show()

    plt.gca().yaxis.set_minor_formatter(NullFormatter())
    plt.subplots_adjust(top=1, bottom=0, left=0, right=1, hspace=0.25,
                        wspace=0.35)
    fig = plt.gcf()


    return fig


# def stateActiveInteraction (flow) :


def activationFlowTable (value_source, genes_name_list):
    all_values = []
    for i in range(len(genes_name_list)):
        pair_values = []
        pair_values.append(genes_name_list[i])
        pair_values.append(value_source[i])
        all_values.append(pair_values)

    print(all_values)
    return all_values





def getFlow1 (flow):
    for key, value in flow.items():
        value = 0
        value_source = []
        all_values = []
        source = key
        target = value
        print ("cle : ",source)
        print("valeur : ", target)
        value += 1
        for number in source :
            value_source.append(number)
        if value == 1 :
            return source, target, value_source




def getAutoregulate(network_as_list, genes_names_list):
    positive_gene = []
    negative_gene = []
    normal_gene = []

    for i in range (len(network_as_list)):
        gene_source = network_as_list[i][0]
        gene_target = network_as_list[i][2]
        interaction = network_as_list[i][1]
        if gene_target == gene_source:
            if interaction == "1":
                positive_gene.append(gene_source)

            else :
                negative_gene.append(gene_source)

    for gene in genes_names_list:
        if gene not in positive_gene and gene not in negative_gene:
            normal_gene.append(gene)






    return positive_gene, negative_gene,normal_gene

def addEdges(network_as_list, G):
    for i in range (len(network_as_list)):
        duet = "0"
    gene_source = network_as_list[i][0]
    gene_target = network_as_list[i][2]
    interaction = network_as_list[i][1]
    for k in range (len(network_as_list)) :
        if network_as_list[k][0] == gene_target and network_as_list[k][2] == gene_source :
            duet = "1"
    G.add_edge(gene_source, gene_target, arrowstyle = interaction, duet = duet)

    return G


def addNodes(positive_gene,negative_gene,normal_gene, G):
    if len(positive_gene) != 0:
        interaction = "1"
        for i in range (len(positive_gene)):
            G.add_node(positive_gene[i], forme = interaction)


    if len(normal_gene) != 0:
        interaction = "0"
        for i in range (len(normal_gene)):
            G.add_node(normal_gene[i], forme = interaction)

    if len(negative_gene) != 0:
        interaction = "-1"
        for i in range (len(negative_gene)):
            G.add_node(negative_gene[i], forme = interaction)

    return G
# ------------------------------- END OF YOUR MATPLOTLIB CODE -------------------------------

# ------------------------------- Beginning of Matplotlib helper code -----------------------
######Ca c'est bon#########

def draw_figure(canvas, figure, loc=(0, 0)):
    figure_canvas_agg = FigureCanvasTkAgg(figure, canvas)
    figure_canvas_agg.draw()
    figure_canvas_agg.get_tk_widget().pack(side='top', fill='both', expand=1)
    return figure_canvas_agg
# ------------------------------- Beginning of GUI CODE -------------------------------
def on_pick(event):
    artist = event.artist
    xmouse, ymouse = event.mouseevent.xdata, event.mouseevent.ydata
    x, y = artist.get_xdata(), artist.get_ydata()
    ind = event.ind
    print ('Artist picked:', event.artist)

