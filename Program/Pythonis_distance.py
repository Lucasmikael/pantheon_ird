import os
import csv
import re
import itertools
import time
import datetime
import numpy as np
import scipy as sp
import random as rd
import matplotlib.pyplot as plt
from numpy import array
from pythonis_tools import *
from pythonis_filesIO import *
from pythonis_model import *
from scipy.spatial import distance


def fonctionMain(working_directory, genes_list_file, network_structure_file, nb_columns_genes=2, name_index=0,
                 nb_columns_network=3, network_headers=1, starting_sample_size=100):
    N = starting_sample_size

    os.chdir(working_directory)
    a, b, c = ImportBooleanModel(working_directory, genes_list_file, network_structure_file, nb_columns_genes,
                                 name_index,
                                 nb_columns_network, network_headers)

    gene_list = a

    parameter_list = [('logical', 'transient'), ('logical', 'constant'), ('algebraic', 'transient'),
                      ('algebraic', 'constant')]

    starting_states = []
    reference_stable_state = []
    altered_stable_states = []

    for i in range(N):
        state = InitializeState(gene_list)  # by default run as random state choice
        N.append(state)

    for (param1, param2) in parameter_list:
        return

    flow, stable, start, stable_states, initial_states, genes_names, data, network, time = RunBooleanModel(a, b, 't',
                                                                                                           't',
                                                                                                           initial_state_number=100,
                                                                                                           initial_state_choice='random',
                                                                                                           stimulus='transient',
                                                                                                           initial_state_genes=[
                                                                                                               'foo'],
                                                                                                           model='algebraic')

    resMod(a, b, start, flow, stable, genes_list_file, network_structure_file, genes_names, network, initial_states,
           stable_states, data)
    # listing(stable_states, data, stable)
    # performance(time, initial_state_number=50000)

    return a, b, flow, stable, start


#######
#
# Main program
#
#######

if __name__ == "__main__":
    a, b, flow, stable, start = fonctionMain("C:\Pythonis\Files_to_run".decode("utf-8"), "genes_list_os.txt",
                                             "ReseauConsolideLitterature.sif", nb_columns_genes=1, name_index=1,
                                             nb_columns_network=3)

"""

Script to find most influential genes or intteractions in a network for a collection of models 

    generate a collection of N starting state (as a list) for the given network
    
    for each set of parameters in 'list of parameters', 
    
        for each starting state S in the list, 
        
            run the simulation from state S and record the reference final state R (single or loop) in a table
            
            for each gene in the network :
                
                set the gene to 0 permanently, then run the simulation from state S again, find the final state F, 
                compute the distance between F and R and record this distance as well as the identity of the mutation as an additional entry in a list in a dictionary (key being name of gene and type of mutation)
                
                set the gene to 1 permanently and do the same thing as previously
                
            for each interaction in the network :
                
                remove the interaction, then run the simulation from state S again, find the final state F, 
                compute the distance between F and R and record this distance as well as the identity of the mutation as an additional entry in a list in a dictionary (key being name of gene and type of mutation)       
 
        
        for gene mutation and for edge mutation 
            compute the mean distance for each mutation tested (mean of the list of distance found from N different starting states)
            sort the list from max to min distance
            
        Record the sorted list of distances and mutation as well as the set of parameters used for the computation  

"""
