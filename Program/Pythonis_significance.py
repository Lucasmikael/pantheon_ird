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


def fullKOOAStudy(working_directory, genes_list_file, network_structure_file, nb_columns_genes=2, name_index=0,
                  nb_columns_network=3, network_headers=0, starting_sample_size=100):
    start_time = time.time()

    N = starting_sample_size

    verbose_output = False
    os.chdir(working_directory)
    gene_list, network, unused_genes = ImportBooleanModel(working_directory, genes_list_file, network_structure_file,
                                                          nb_columns_genes, name_index, nb_columns_network,
                                                          network_headers)

    parameter_list = [('logical', 'transient'), ('logical', 'constant'), ('algebraic', 'transient'),
                      ('algebraic', 'constant')]

    starting_states_set = []
    reference_stable_state = []
    altered_stable_states = {}
    output = {}
    meansKO = {}
    meansOA = {}
    for gene in gene_list:
        output[gene] = []
        meansKO[gene] = 0
        meansOA[gene] = 0

    for i in range(N):
        state = InitializeState(gene_list)  # by default run as random state choice
        starting_states_set.append(state)

    print
    '########################\n#'
    for (param1, param2) in parameter_list:

        run_number = 0
        runcount = 0
        tenpercentrun = len(starting_states_set) / 10

        print
        '# Running model set : ', param1, param2
        print
        '#'

        for starting_point in starting_states_set:

            run_number += 1
            runcount += 1
            if (len(starting_states_set) > 10):
                print
                "# Running - %.0f%% of all runs done" % (100 * float(runcount) / float(len(starting_states_set)))
            else:
                print
                "# Running - run", runcount, "out of", len(starting_states_set)

            flow, stable, start, stable_states, initial_states, genes_names, data, networkoutput, runtime = RunBooleanModel(
                gene_list, network, 't', 't',
                initial_state_number='specified', initial_state_choice='random',
                stimulus=param2, initial_state_genes=['foo'], model=param1, specified_starting_states=[starting_point],
                verbose=verbose_output)

            for k, v in stable.items():
                reference_stable_state = k  # only a single pair of key,value in the stable states dictionary when the model is run from a single starting state

            genecount = 0
            tenpercentgene = len(gene_list) / 10
            for gene in gene_list:
                genecount += 1
                if (len(gene_list) > 10):
                    if not (genecount % tenpercentgene):
                        print
                        "# Mutating - %.0f%% of all genes done" % (100 * float(genecount) / float(len(gene_list)))
                else:
                    print
                    "# Mutating - gene", genecount, "out of", len(gene_list)

                # run all possible mutated network for gene KO or overexpressor
                flow, stable_KO, start, stable_states, initial_states, genes_names, data, networkoutput, runtime = RunBooleanModel(
                    gene_list, network, 't', 't',
                    initial_state_number='specified', initial_state_choice='random',
                    stimulus=param2, initial_state_genes=['foo'], model=param1,
                    specified_starting_states=[starting_point], KO_genes=gene, verbose=verbose_output)

                flow, stable_OA, start, stable_states, initial_states, genes_names, data, networkoutput, runtime = RunBooleanModel(
                    gene_list, network, 't', 't',
                    initial_state_number='specified', initial_state_choice='random',
                    stimulus=param2, initial_state_genes=['foo'], model=param1,
                    specified_starting_states=[starting_point], OA_genes=gene, verbose=verbose_output)

                stable_KO_state = stable_KO.keys()  # only a single key from starting from a single state
                stable_OA_state = stable_OA.keys()  # only a single key from starting from a single state

                altered_stable_states[gene] = ('KO', stable_KO_state[0]), ('OA', stable_OA_state[0])

                # compute the distance between the reference stable state and the ones reached after mutations of the network
                # first check if either stable states are cycles or not

                stable_ref_is_cyclic = any(isinstance(i, tuple) for i in reference_stable_state)
                stable_KO_is_cyclic = any(isinstance(i, tuple) for i in stable_KO_state[0])
                stable_OA_is_cyclic = any(isinstance(i, tuple) for i in stable_OA_state[0])

                if stable_ref_is_cyclic:

                    if stable_KO_is_cyclic:  # two cycles, we compute the distance between each state pair and record the biggest distance
                        dist_list = []
                        for state in reference_stable_state:
                            for m_state in stable_KO_state[0]:
                                dist_list.append(distance.hamming(state, m_state))
                        dist_KO = max(dist_list)

                    else:  # single stable mutant state against cyclic refence stable state - we compare the mutant to each state of the cycle and record to biggest distance
                        dist_list = []
                        for state in reference_stable_state:
                            dist_list.append(distance.hamming(state, stable_KO_state[0]))
                        dist_KO = max(dist_list)

                    if stable_OA_is_cyclic:
                        dist_list = []
                        for state in reference_stable_state:
                            for m_state in stable_OA_state[0]:
                                dist_list.append(distance.hamming(state, m_state))
                        dist_OA = max(dist_list)

                    else:  # single stable mutant state against cyclic refence stable state - we compare the mutant to each state of the cycle and record to biggest distance
                        dist_list = []
                        for state in reference_stable_state:
                            dist_list.append(distance.hamming(state, stable_OA_state[0]))
                        dist_OA = max(dist_list)

                else:  # single reference stable state

                    if stable_KO_is_cyclic:  # cyclic mutant stable state against single refence stable state - we compare the reference to each state of the mutant cycle and record to biggest distance
                        dist_list = []
                        for state in stable_KO_state[0]:
                            dist_list.append(distance.hamming(state, reference_stable_state))
                        dist_KO = max(dist_list)

                    else:  # both ref and altered states are single
                        dist_KO = distance.hamming(reference_stable_state, stable_KO_state[0])

                    if stable_OA_is_cyclic:  # cyclic mutant stable state against single refence stable state - we compare the reference to each state of the mutant cycle and record to biggest distance
                        dist_list = []
                        for state in stable_OA_state[0]:
                            dist_list.append(distance.hamming(state, reference_stable_state))
                        dist_OA = max(dist_list)

                    else:  # both ref and altered states are single
                        dist_OA = distance.hamming(reference_stable_state, stable_OA_state[0])

                # record the run result in the output dictionaries
                meansKO[gene] = meansKO[gene] + dist_KO
                meansOA[gene] = meansOA[gene] + dist_OA
                output[gene] = output[gene] + [(param1, param2, run_number, dist_KO, dist_OA)]

            # resMod(a, b, start, flow, stable, genes_list_file, network_structure_file, genes_names, network, initial_states, stable_states, data)

    print
    '############# Results ###################\n#'
    for k, v in meansKO.items():
        total_nb_runs = len(parameter_list) * N
        meansKO[k] = meansKO[k] / total_nb_runs
        meansOA[k] = meansOA[k] / total_nb_runs
        print
        '# gene :', k, '- mean KO distance :', meansKO[k], '- mean OA distance :', meansOA[
            k], ' - over', total_nb_runs, 'run combinations.'

    end_time = time.time() - start_time
    print
    '# Running time (seconds) %s', end_time
    print
    '#\n#########################################'

    return gene_list, output, meansKO, meansOA


def targetedKOStudy(working_directory, genes_list_file, network_structure_file, nb_columns_genes=2, name_index=0,
                    nb_columns_network=3, network_headers=0, starting_sample_size=100, KO_genes_list=['foo']):
    start_time = time.time()

    N = starting_sample_size

    verbose_output = False
    os.chdir(working_directory)
    gene_list, network, unused_genes = ImportBooleanModel(working_directory, genes_list_file, network_structure_file,
                                                          nb_columns_genes, name_index, nb_columns_network,
                                                          network_headers)

    parameter_list = [('logical', 'transient'), ('logical', 'constant'), ('algebraic', 'transient'),
                      ('algebraic', 'constant')]

    starting_states_set = []
    reference_stable_state = []
    altered_stable_states = {}
    output = {}
    meansKO = {}

    for KO_set in KO_genes_list:
        output[str(KO_set)] = []
        meansKO[str(KO_set)] = 0

    for i in range(N):
        state = InitializeState(gene_list)  # by default run as random state choice
        starting_states_set.append(state)

    print
    '########################\n#'
    for (param1, param2) in parameter_list:

        run_number = 0
        runcount = 0
        tenpercentrun = len(starting_states_set) / 10

        print
        '# Running model set : ', param1, param2
        print
        '#'

        for starting_point in starting_states_set:

            run_number += 1
            runcount += 1
            if (len(starting_states_set) > 10):
                print
                "# Running - %.0f%% of all runs done" % (100 * float(runcount) / float(len(starting_states_set)))
            else:
                print
                "# Running - run", runcount, "out of", len(starting_states_set)

            flow, stable, start, stable_states, initial_states, genes_names, data, networkoutput, runtime = RunBooleanModel(
                gene_list, network, 't', 't',
                initial_state_number='specified', initial_state_choice='random',
                stimulus=param2, initial_state_genes=['foo'], model=param1, specified_starting_states=[starting_point],
                verbose=verbose_output)

            for k, v in stable.items():
                reference_stable_state = k  # only a single pair of key,value in the stable states dictionary when the model is run from a single starting state

            genecount = 0
            tenpercentmutants = len(KO_genes_list) / 10

            for KO_set in KO_genes_list:
                genecount += 1
                if (len(KO_genes_list) > 10):
                    if not (genecount % tenpercentmutants):
                        print
                        "# Mutating - %.0f%% of all KO set done" % (100 * float(genecount) / float(len(KO_genes_list)))
                else:
                    print
                    "# Mutating - KO set", genecount, "out of", len(KO_genes_list)

                # run all gene set KO
                flow, stable_KO, start, stable_states, initial_states, genes_names, data, networkoutput, runtime = RunBooleanModel(
                    gene_list, network, 't', 't',
                    initial_state_number='specified', initial_state_choice='random',
                    stimulus=param2, initial_state_genes=['foo'], model=param1,
                    specified_starting_states=[starting_point], KO_genes=KO_set, verbose=verbose_output)

                stable_KO_state = stable_KO.keys()  # only a single key from starting from a single state

                altered_stable_states[str(KO_set)] = ('KO', stable_KO_state[0])

                # compute the distance between the reference stable state and the ones reached after mutations of the network
                # first check if either stable states are cycles or not

                stable_ref_is_cyclic = any(isinstance(i, tuple) for i in reference_stable_state)
                stable_KO_is_cyclic = any(isinstance(i, tuple) for i in stable_KO_state[0])

                if stable_ref_is_cyclic:

                    if stable_KO_is_cyclic:  # two cycles, we compute the distance between each state pair and record the biggest distance
                        dist_list = []
                        for state in reference_stable_state:
                            for m_state in stable_KO_state[0]:
                                dist_list.append(distance.hamming(state, m_state))
                        dist_KO = max(dist_list)

                    else:  # single stable mutant state against cyclic refence stable state - we compare the mutant to each state of the cycle and record to biggest distance
                        dist_list = []
                        for state in reference_stable_state:
                            dist_list.append(distance.hamming(state, stable_KO_state[0]))
                        dist_KO = max(dist_list)


                else:  # single reference stable state

                    if stable_KO_is_cyclic:  # cyclic mutant stable state against single refence stable state - we compare the reference to each state of the mutant cycle and record to biggest distance
                        dist_list = []
                        for state in stable_KO_state[0]:
                            dist_list.append(distance.hamming(state, reference_stable_state))
                        dist_KO = max(dist_list)

                    else:  # both ref and altered states are single
                        dist_KO = distance.hamming(reference_stable_state, stable_KO_state[0])

                # record the run result in the output dictionaries
                meansKO[str(KO_set)] = meansKO[str(KO_set)] + dist_KO

                output[str(KO_set)] = output[str(KO_set)] + [(param1, param2, run_number, dist_KO)]

            # resMod(a, b, start, flow, stable, genes_list_file, network_structure_file, genes_names, network, initial_states, stable_states, data)

    print
    '############# Results ###################\n#'
    for k, v in meansKO.items():
        total_nb_runs = len(parameter_list) * N
        meansKO[k] = meansKO[k] / total_nb_runs
        print
        '# gene :', k, '- mean KO distance :', meansKO[k], ' - over', total_nb_runs, 'run combinations.'

    end_time = time.time() - start_time
    print
    '# Running time (seconds) %s', end_time
    print
    '#\n#########################################'

    return gene_list, output, meansKO


def generate_pairs(source):
    pairs = []
    for p1 in range(len(source)):
        for p2 in range(p1 + 1, len(source)):
            pairs.append([source[p1], source[p2]])
    return pairs


# check itertools to generate combinations list
# >>> import itertools
# >>> bills = [20, 20, 20, 10, 10, 10, 10, 10, 5, 5, 1, 1, 1, 1, 1]
# >>> list(itertools.combinations(bills, 3))
# [(20, 20, 20), (20, 20, 10), (20, 20, 10), ... ]


#######
#
# Main program
#
#######

if __name__ == "__main__":

    rice_filenames = ['ReseauConsolideLitterature.sif', 'unionnetwork.sif',
                      'liste215_sansprior_moy3rep.sif']  # rice run - output in C:\Pythonis\Files_to_run\significance and gene list genes_list_os.txt
    arabido_filenames = ['TDCor6.32_output_221018_parallel.txt']  # arabido run
    PLTARFKO = [['PLT1'], ['PLT2'], ['PLT3'], ['PLT5'], ['PLT7'], ['ARF6'], ['ARF8'], ['PLT1', 'PLT2'],
                ['PLT5', 'PLT7'], ['ARF6', 'ARF8'], ['PLT1', 'PLT2', 'PLT3'], ['PLT3', 'PLT5', 'PLT7']]
    KO_genes_input = ['PLT1', 'ARF6', 'LRP1', 'PHB', 'TMO5', 'SHR', 'SCR', 'SHP1', 'ATML1', 'PID2']
    OA_genes_input = ['PLT7', 'PUCHI', 'CRF1', 'ARF2', 'PLT5', 'ARF9', 'ARF17', 'U.box', 'ARF19', 'WRKY43']
    sample_size = 100

    for filename in arabido_filenames:
        # gene_list, output, meansKO, meansOA = fullKOOAStudy("C:\Work\Papers\PYTHONIS\LRPnetwork-201903".decode("utf-8"), "l_gnp-221018-161genes-Ju+PI-prior.txt", filename, nb_columns_genes=1, name_index=1, nb_columns_network=3, starting_sample_size=sample_size)
        # WriteDistanceToCsv('GS_' + filename, gene_list, meansKO, meansOA, sample_size)

        gene_list, output, meansKO = targetedKOStudy("C:\Work\Papers\PYTHONIS\LRPnetwork-201903".decode("utf-8"),
                                                     "l_gnp-221018-161genes-Ju+PI-prior.txt", filename,
                                                     nb_columns_genes=1, name_index=1, nb_columns_network=3,
                                                     starting_sample_size=sample_size, KO_genes_list=PLTARFKO)

        # gene_list, output, meansKO, meansOA = fonctionMain("C:\Pythonis\Files_to_run\significance".decode("utf-8"), "test_gene_list.txt", "test_gene_network.txt", nb_columns_genes=1, name_index=1, nb_columns_network=3, starting_sample_size=sample_size)

        WriteTargetedKODistanceToCsv('GS_' + filename, PLTARFKO, meansKO, sample_size)

# if __name__ == "__main__":
#   output = fonctionMain("C:\Pythonis\Files_to_run\significance".decode("utf-8"), "test_gene_list.txt",
#                                     "test_gene_network.txt", nb_columns_genes=1, name_index=1, nb_columns_network=3, starting_sample_size=50)


"""

Script to find most influential genes or interactions in a network for a collection of models 

    generate a collection of N starting state (as a list) for the given network
    
    for each set of parameters in 'list of parameters', 
    
        for each starting state S in the list, 
        
            run the simulation from state S and record the reference final state R (single or loop) in a table
            
            for each gene in the network :
                
                set the gene to 0 permanently, then run the simulation from state S again, find the final state F, 
                compute the distance between F and R and record this distance as well as the identity of the mutation as an additional entry in a list in a dictionary (key being name of gene and type of mutation)
                
                set the gene to 1 permanently and do the same thing as previously
                
            for each interaction in the network :  # TO BE DONE
                
                remove the interaction, then run the simulation from state S again, find the final state F, 
                compute the distance between F and R and record this distance as well as the identity of the mutation as an additional entry in a list in a dictionary (key being name of gene and type of mutation)       
 
        
        for gene mutation (and for edge mutation - TO BE DONE)
            compute the mean distance for each mutation tested (mean of the list of distance found from N different starting states)
            sort the list from max to min distance
            
        Record the sorted list of distances and mutation as well as the set of parameters used for the computation  

"""
