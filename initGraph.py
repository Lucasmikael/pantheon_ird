import networkx as nx
from pythonis_model import *
import matplotlib.pyplot as plt

def drawGraph(genes_names_list, network_as_list):

    G = nx.MultiDiGraph()

    G.add_nodes_from(genes_names_list)
    val_map = {'ARR23': 1.0,
               'CRF2': 0.5714285714285714}

    values = [val_map.get(node, 0.85) for node in G.nodes()]

    #ajout des aretes et des noeuds
    for i in range (len(network_as_list)):
        gene_source = network_as_list[i][0]
        gene_target = network_as_list[i][2]
        interaction = network_as_list[i][1]
        G.add_edge(gene_source, gene_target, arrowstyle = interaction)

        # if gene_source not in G.nodes :
        #     if gene_target == gene_source :
        #         G.add_node(gene_source, regulation = "auto")
        #     else:
        #         G.add_node(gene_source, regulation = "none")


    pos = nx.circular_layout(G)
    nx.draw_networkx_nodes(G, pos, node_size=1500, node_color= values)
    nx.draw_networkx_labels(G, pos)
    for edge in G.edges(data=True):
        if edge[2]['arrowstyle'] == "-1":
            arrowstyle = "-["
        else :
            arrowstyle = "-|>"
        nx.draw_networkx_edges(G, pos, edgelist=[(edge[0],edge[1])], arrowstyle = arrowstyle, node_size=1900)

    plt.show()


# def getReseau(network_dictionary):
#
if __name__ == "__main__":

    # Path des fichiers en Import
    genes_list_filename = "l_gnp-190515.txt"
    network_filename = "jeu_essai_helios.txt"

    # Values en dur pour le run du modèle booléen
    boolean_model_selected = "logical"
    boundary_model_selected = "transient"
    initial_state_selected = "random"
    number_initial_states = "1"
    KO_genes_selected = ['foo']
    OA_genes_selected = ['foo']
    genes_selected = ['foo']

    aretes = []


    genes_names_list, network_dictionary, unused_genes, network_as_list = ImportBooleanModel(genes_list_filename,
                                                                                             network_filename)


    drawGraph(genes_names_list,network_as_list)





    # #Identifier les duets du fichier
    # for i in range (len(network_as_list)) :
    #         gene_source = network_as_list[i][0]
    #         gene_target = network_as_list[i][2]
    #         for k in range (len(network_as_list)) :
    #              if network_as_list[k][0] == gene_target and network_as_list[k][2] == gene_source :
    #                 print (gene_source, "and ", gene_target)

    # for i in range(len(network_as_list)):
    #     # source_target = (network_as_list[i][0], network_as_list[i][2], network_as_list[i][1])
    #     # aretes.append(source_target)
    #     print("source : ", network_as_list[i][0])
    #     print("target : ", network_as_list[i][2])
    #     print("interaction : ", network_as_list[i][1])
    #     g.add_edge(network_as_list[0], network_as_list[2], arrowstyle= network_as_list[1])
    #
    # g.add_nodes_from(genes_names_list)
    # # g.add_edges_from(aretes)
    # print(g.edges)
    #
    #
    # # Defini les couleurs de noeuds
    # val_map = {'ARR23': 1.0,
    #            'CRF2': 0.5714285714285714}
    #
    # values = [val_map.get(node, 0.85) for node in g.nodes()]
    #
    # for i in range (len(g.edges)) :
    #     if aretes[i][2] == '-1':
    #         arrowstyle.append( "-[")
    #
    #     else:
    #         arrowstyle.append("-|>")
    #
    #
    # print(g.edges)
    # pos = nx.circular_layout(g)
    # for node in g.nodes(data = True):
    #     nx.draw_networkx_nodes(g, pos, node_shape='s', node_size=2500, node_color=values)
    # nx.draw_networkx_labels(g, pos)
    # # nx.draw_networkx_edges(g, pos, arrows=True,arrowsize= 20, arrowstyle='-|>', connectionstyle='arc3, rad = 0.1',node_size= 2500)
    # for edge in g.edges(data = True):
    #     if edge[2] == '-1' :
    #         arrowstyle = "-["
    #     else:
    #         arrowstyle = "-|>"
    #     nx.draw_networkx_edges(g, pos, edgelist=[(edge[0], edge[1])], arrows=True,arrowsize= 20, arrowstyle=arrowstyle, connectionstyle='arc3, rad = 0.1', node_size = 2600)
    #
    #
    #
    #
    #
    #
    #
    # # test separation du Draw
    # # for edge in g.edges(data=True):
    # #     # w = edge[2]['w']
    # #     nx.draw_networkx_edges(g, pos, edgelist=[(edge[0], edge[1])], arrowsize=800, node_size=300)
    # #
    # #
    # # nx.draw(g, pos=nx.circular_layout(g), node_shape='s', node_size=2000, node_color=values, with_labels=True,
    # #         arrows=True, arrowstyle='-|>', connectionstyle='arc3, rad = 0.1')
    #
    #
    # # Pour la visualistaion
    # # nx.draw(g, pos=nx.circular_layout(g),with_labels=True,arrows = True, connectionstyle='arc3)
    # plt.show()

    # flow, stable, start, stable_states, initial_states, genes_names, data, network, time = RunBooleanModel(
    #     genes_names=genes_names_list, genes_network=network_dictionary,
    #     initial_state_number=number_initial_states, initial_state_choice=initial_state_selected,
    #     model=boolean_model_selected, stimulus=boundary_model_selected, initial_state_genes=genes_selected,
    #     KO_genes=KO_genes_selected, OA_genes=OA_genes_selected)
    #
    # g = nx.Graph()
    #
    # g.add_nodes_from(genes_names_list)

