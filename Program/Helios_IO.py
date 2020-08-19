import csv
from Helios_model import *
from pythonis_model import *
import networkx as nx


def saveData(network_as_list, flow, genes_names_list,name_saved_file):
    file = open(name_saved_file+".csv", 'w')
    interaction = 0
    print(genes_names_list)
    with file:
        writer = csv.writer(file)
        for i in range(len(flow)):
            G = nx.MultiDiGraph()
            value_source = getFlow(flow, i)
            global_gene_state = getRegulationActivation(network_as_list, genes_names_list, value_source)
            G = addNodes(global_gene_state, G)
            G = addEdges(network_as_list, G)
            if interaction == 0:
                writer.writerow(["Interaction"])
                for edge in G.edges(data=True):
                    source = edge[0]
                    target = edge[1]
                    arrow = edge[2]['arrowstyle']
                    duet = edge[2]['duet']
                    writer.writerow([source, target, arrow, duet])
            interaction += 1
            writer.writerow(["Etat", i])
            for node in G.nodes(data=True):
                name = node[0]
                form = node[1]['forme']
                state = node[1]['active_state']
                writer.writerow([name, form, state])
            writer.writerow(["Fin"])


def openData(load_saved_file):
    file = open(load_saved_file, "r")
    reader = csv.reader(file, delimiter=",")
    interaction = False
    state = False
    number_state = -1
    interaction_table = []
    node_table = []
    list_panel_loading = []
    start_loading_list_gene = 0


    for row in reader:
        if row[0] == "Fin":
            node_table.append(node_state_table)
            state = False
            # print(node_table)
        if row[0] == "Etat" and row[1] != number_state:
            state = False
            node_state_table = []
            number_state = row[1]
        # if row[0] == "Fin":
        #     node_table.append(node_state_table)
        if state:
            row_table = []
            row_table.append(row[0])
            row_table.append(row[1])
            row_table.append(row[2])
            node_state_table.append(row_table)
        if row[0] == "Etat":
            interaction = False
            state = True
        if interaction:
            row_table = []
            row_table.append(row[0])
            row_table.append(row[1])
            row_table.append(row[2])
            row_table.append(row[3])
            interaction_table.append(row_table)
        if row[0] == "Interaction":
            interaction = True
    print(node_table[0][0][0])

    for i in range (len(node_table)) :
        list_panel_loading.append(i)

    print(list_panel_loading)
    file.close()
    return interaction_table, node_table, list_panel_loading


def drawLoadGraph(layout_selected, interaction_table, node_table, graph_selected, color_activate_node,
                  color_inactivate_node, color_active_edge, color_inactivate_edge, activate_widthedge, genes_selected_visu, state):
    G = nx.MultiDiGraph()
    G, global_gene_state = addNodesImport(node_table, graph_selected, G)
    G = addEdgesImport(interaction_table, G)
    fig, G = drawFig(G, global_gene_state, layout_selected, color_activate_node, color_inactivate_node,
                     color_active_edge, color_inactivate_edge,activate_widthedge, genes_selected_visu, state)

    return fig, G


def addEdgesImport(interaction_table, G):
    for edge in interaction_table:
        G.add_edge(edge[0], edge[1], arrowstyle=edge[2], duet=edge[3])
        # print(edge[0],edge[1],edge[2],edge[3])
    return G


def addNodesImport(node_table, graph_selected, G):
    global_gene_state = []
    for graph in range(len(node_table)):
        if graph == graph_selected:
            for node in range(len(node_table[graph])):
                new_node = []
                # print(node_table[graph][node][0],node_table[graph][node][1],node_table[graph][node][2])
                G.add_node(node_table[graph][node][0], forme=node_table[graph][node][1],
                           active_state=int(node_table[graph][node][2]))
                new_node.append(node_table[graph][node][0])
                new_node.append(int(node_table[graph][node][2]))
                global_gene_state.append(new_node)
    return G, global_gene_state
