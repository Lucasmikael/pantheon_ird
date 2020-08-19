import re
from Helios_model import *
from pythonis_model import *


def addElementfromCSV(name_file_addelement):
    file = open(name_file_addelement, "r")
    reader = csv.reader(file, delimiter=",")
    # interaction = False
    # state = False
    # number_state = -1
    next(reader)
    new_element = []
    source_gene_list = []
    interaction_list= []
    target_gene_list = []

    for row in reader:
        row_table = []
        source_gene_list.append(row[0])
        interaction_list.append(row[1])
        target_gene_list.append(row[2])

    file.close()
    return source_gene_list,interaction_list, target_gene_list


def catchGUIElement(source_gene_add,interaction_add, target_gene_add):
    source_gene_list = []
    interaction_list = []
    target_gene_list = []
    source_gene_list.append(source_gene_add)
    interaction_list.append(interaction_add)
    target_gene_list.append(target_gene_add)

    return source_gene_list, interaction_list, target_gene_list



def addElement(genes_network, genes_names_list, network_as_list, source_gene_list, interaction_list,target_gene_list):
    regex = re.compile('[^a-zA-Z0-9 ]')
    for i in range (len(source_gene_list)):
        source_gene, target_gene = regex.sub('.', source_gene_list[i]), regex.sub('.', target_gene_list[i])
        interaction_gene = interaction_list[i]
        network_as_list.append([source_gene, interaction_gene, target_gene])
        try:
            genes_network[source_gene].append([interaction_gene, target_gene])
        except KeyError:
            genes_network[source_gene] = [[interaction_gene, target_gene]]
        if (source_gene) not in set(genes_names_list):
            genes_names_list.append(source_gene_list[i])
        if (target_gene) not in set(genes_names_list):
            genes_names_list.append(target_gene_list[i])
    print(genes_names_list)
    print(genes_network)
    return genes_names_list, genes_network, network_as_list
