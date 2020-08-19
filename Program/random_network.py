# -*- coding: utf-8 -*-

import os
import networkx as nx
import random
import matplotlib.pyplot as plt
import itertools

## Creation d'interactions aléatoires entre 1 et -1
## J'ai un une boucle avec une condition pour eviter d'avoir le chiffre 0

def get_random_interactions():

   k=0
   while k >= 0:
    random_number = random.randint(-1, 1)
    if random_number == 0:
       k = k+1
    else:
        return random_number

## Création d'associations aléatoires entre les genes
## Par défaut, le nombre de paires (d'interactions) souhaité est 'all'

def get_random_pairs(genes, number_of_pairs = 'all'):
    # Genere toutes les associations possibles
    a = list(itertools.permutations(genes, 2))    #itertools.permutations implique les associations du type (6,5) et (5,6)
    b = list(itertools.combinations_with_replacement(genes, 2))       #.combinations implique les associations du type (6,6) , (5,5) etc
    c = a + b
    c = list(set(c))      ##Enlever les doublons
    # Mélange alééatoirement les associations
    random.shuffle(c)
    if number_of_pairs == 'all':
         pairs = c
    else:
         pairs = random.sample(c, number_of_pairs)      # Prendre un echantillon des associations quand c'est spécifié

    return pairs


def RandomGRN1(working_directory, genes_list_file, network_structure_file, number_of_nodes, number_of_edges):
    os.chdir(working_directory)



    ## Creation d'une liste de genes en fonction du nombre de noeuds donne

    j = 1
    list_genes = []
    string = 'G'
    while j <= number_of_nodes:
        list_genes.append(string + str(j))  # test[j] n'existe pas en python. On utilise append pour ajouter en fin de liste.
        j = j + 1

    ## Creation du fichier contenant la liste de genes qui est variable

    grn = open(genes_list_file, 'w')
    for i in list_genes:
        grn.write(i)
        grn.write('\n')

    grn.close()

    ## Creation du fichier contenant les interactions de maniere aleatoire

        ## Creation de la liste d'associations entre les gènes en fonction du nombre d'aretes souhaité

    paires = get_random_pairs(list_genes, number_of_pairs=number_of_edges)
    print(str(paires))
        ## Creation du fichier avec les interactions

    ntwk = open(network_structure_file, 'w')
    ntwk.write("source\tinter\ttarget\n")
    for j in paires:
           str_pairs = str(j)     ## les listes ne peuvent pas être splitées. On convertit en str
           a = str_pairs.split("'")      ## On sépare en fonction des ' ex= ('G1','G2') splité donne ( G1 , G2 et )
           ntwk.write((a[1]))                                      ## a[1] correspond au 2 élément splité, ici c'est G1
           ntwk.write("\t")
           ntwk.write(str(get_random_interactions()))
           ntwk.write("\t")
           ntwk.write((a[3]))
           ntwk.write('\n')

    ntwk.close()

        ## Creation et visualisation  du graphe dirigé avec les associations de gènes

    g = nx.MultiDiGraph()
    g.add_nodes_from(list_genes)
    g.add_edges_from(paires)

    inter = get_random_interactions()
    print(inter)
    print("Actual number of edges: " + str(g.number_of_edges()))

    # Pour la visualistaion
    nx.draw(g, pos=nx.circular_layout(g))
    plt.show()

    return


if __name__ == "__main__":

    RandomGRN1("C:\Pythonis\Files_to_run".decode("utf-8"), genes_list_file='TEST.txt',
               network_structure_file='TEST1.txt', number_of_nodes=50, number_of_edges=1000)
