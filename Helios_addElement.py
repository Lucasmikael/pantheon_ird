import re
from Helios_model import *
from pythonis_model import *


def addElement(genes_network, genes_names_list, network_as_list):
    regex = re.compile('[^a-zA-Z0-9 ]')
    new_element = ["ARR23", "1", "NEWELEMENT"]

    source_gene, interaction, target_gene = new_element[0], new_element[1], new_element[2]
    source_gene, target_gene = regex.sub('.', source_gene), regex.sub('.', target_gene)

    print("source :", source_gene, " interaction", interaction, " Target:", target_gene)
    network_as_list.append([source_gene, interaction, target_gene])
    try:
        genes_network[source_gene].append([interaction, target_gene])
    except KeyError:
        genes_network[source_gene] = [[interaction, target_gene]]
    if (source_gene) not in set(genes_names_list):
        genes_names_list.append(source_gene)
    if (target_gene) not in set(genes_names_list):
        genes_names_list.append(target_gene)
    print(genes_names_list)
    print(genes_network)
    return genes_names_list, genes_network, network_as_list
