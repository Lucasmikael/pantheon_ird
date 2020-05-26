import csv
from Helios_model import *
from pythonis_model import *
import networkx as nx

def saveData(network_as_list,flow, genes_names_list):
    file = open("nodegraph.csv", 'w')
    interaction = 0
    with file :
        writer = csv.writer(file)
        for i in range (len(flow)):
            G = nx.MultiDiGraph()
            value_source = getFlow(flow, i)
            global_gene_state = getRegulationActivation(network_as_list, genes_names_list, value_source)
            G = addNodes(global_gene_state, G)
            G = addEdges(network_as_list, G)
            if interaction == 0 :
                writer.writerow(["Interaction"])
                for edge in G.edges(data=True):
                    source = edge[0]
                    target = edge[1]
                    arrow = edge[2]['arrowstyle']
                    duet = edge[2]['duet']
                    writer.writerow([source,target,arrow,duet])
            interaction += 1
            writer.writerow(["Etat",i])
            for node in G.nodes(data=True):
                name= node[0]
                form = node[1]['forme']
                state = node[1]['active_state']
                writer.writerow([name,form,state])
            writer.writerow(["Fin"])




def openData():
    file = open("nodegraph.csv","r")
    reader = csv.reader(file, delimiter = ",")
    interaction = False
    state = False
    number_state = -1
    interaction_table = []
    node_table = []

    for row in reader :
        if row[0] == "Fin":
            node_table.append(node_state_table)
            state = False
        if row[0] == "Etat" and row[1] != number_state:
            state = False
            node_state_table = []
            number_state = row[1]
        if row[0] == "Fin":
            node_table.append(node_state_table)
        if state :
            row_table = []
            row_table.append(row[0])
            row_table.append(row[1])
            row_table.append(row[2])
            node_state_table.append(row_table)
        if row[0] == "Etat":
            interaction = False
            state = True
        if interaction :
            row_table = []
            row_table.append(row[0])
            row_table.append(row[1])
            row_table.append(row[2])
            row_table.append(row[3])
            interaction_table.append(row_table)
        if row[0] == "Interaction":
            interaction = True
    file.close()
    return interaction_table, node_table
