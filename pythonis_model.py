import os
import csv
import re
import itertools
import time
import datetime
import numpy as np
import scipy as sp
import random as rd
from numpy import array
from pythonis_tools import Flatten


def ImportBooleanModel(genes_list_file, network_structure_file, nb_columns_genes=3, name_index=1,
                       nb_columns_network=6, network_headers=1):
    """
    Create the basic datastructure for the boolean model. Designed to work from TDCor output files as base.

    working_directory : the directory where the text files are located
    genes_list_file : multi-columns text file including a gene list. Here used for a text file with c1 : At ID code and c2 : common gene name
    nb_columns_genes : nb of columns in genes_list_file
    name_index : id of the column to use as gene name source from gene_list file - start at 0 for first column
    network_structure_file : tab separated multi-column text file with at least 3 columns which are c1 : source gene name, c2 : interaction weight and nature (number), c3 : target gene name.
    nb_columns_network : nb of columns in genes_network_structure file
    network_headers : presence or absence of a header line in the network file
    The gene names used in the network structure must match those imported from gene_list.

    Return the list of genes names used, the network dictionary {source gene:[interaction, target gene]} and
    the list of genes from the initial list which are not in the network
    """

    #####
    # TODO : assess whether to use numpy array as an alternative for the data loading process
    # w/ genfromtxt(usecols=(x,y,...), names=True or skip_header=(a,b,...)) arguments
    #####

    # get in directory where text files are located - deprecated with gui as we get the absolute file location directly
    # os.chdir(working_directory)

    # open genes list file, read its content and close it to clear memory
    # create a list from the file content and keep only the (name_index)th elements starting at the (name_index)ed one
    # BEWARE - once a file has been read, its pointer sit on empty (end of file)

    f_genes = open(genes_list_file, "r")
    content_genes = f_genes.read()
    f_genes.close()
    content_genes = content_genes.split()  # create a list from each word in the file
    genes_names = []
    genes_names = content_genes[
                  name_index::nb_columns_genes]  # from name_indexed element, take every nb_columned elements
    # use regular expression filtering to replace non-alphanumeric characters in genes names by '.' to avoid name-formatting issue
    regex = re.compile('[^a-zA-Z0-9 ]')
    for gene in genes_names:
        genes_names[genes_names.index(gene)] = regex.sub('.', gene)

    # open network structure file, read its content as a csv table and close it to clear memory
    # BEWARE - the rUniversal option on open call is to take care of windows/mac line termination convention difference
    # that may arise depending on where the txt file was created and may cause issues with csv reader

    f_network = open(network_structure_file, 'rU')
    network_table = list(csv.reader(f_network, delimiter='\t'))  # copy file content as a list
    f_network.close()

    # create a dictionary that will stock the network structure for each of the genes
    # additionnaly manage a list of initial genes and cross-out all genes present in the network
    # BEWARE - to avoid border effects on the initial genes_names list, a copy is created throught the list() call

    genes_network = {}
    genes_network_for_gui = []
    used_genes = []

    if (network_headers):
        del network_table[0]  # remove the header line if there is one

    for line in network_table:
        source_gene, interaction, target_gene = line[0], line[1], line[2]

        # use regular expression filtering to replace non-alphanumeric characters in genes names by '.' to avoid name-formatting issue
        source_gene, target_gene = regex.sub('.', source_gene), regex.sub('.', target_gene)

        genes_network_for_gui.append([source_gene, interaction, target_gene])

        # fill the network dictionnary
        # use a try/except statement to avoid testing the existence of each key and speed up the building process
        # equivalent code in comments afterward for reference
        try:
            genes_network[source_gene].append([interaction, target_gene])
        except KeyError:
            genes_network[source_gene] = [[interaction, target_gene]]
        # if source_gene in genes_network:
        #	genes_network[source_gene].append([interaction, target_gene])
        # else :
        #	genes_network[source_gene] = [[interaction, target_gene]]

        # keep track of used genes
        used_genes.append(source_gene)
        used_genes.append(target_gene)
        # add to the list of genes any that might be present in the network file but not in the initial gene list
        if (source_gene) not in set(genes_names):
            genes_names.append(source_gene)
        if (target_gene) not in set(genes_names):
            genes_names.append(target_gene)

    # substract the set of unique used genes to the the set of unique initial gene names and return the list of unused genes
    absent_genes = list(set(genes_names) - set(used_genes))

    present_genes = list(set(used_genes))
    present_genes.sort()
    genes_network_for_gui.sort()

    return present_genes, genes_network, absent_genes, genes_network_for_gui


def InitializeState(genes_names, initial_state_choice='random', initial_state_genes=['foo']):
    """
    Initialize the gene list state according to the chosen option.

    initial_state_choice can be 'random' (random fill of initial state with 0 and 1), 'all_zeros' (fill initial state with 0), 'all_ones' (fill initial state with 1), 'specified' (set the given genes to 1 and the rest to 0, 'random-specified' (set the given genes to 1 and the rest to random 0 or 1

    Return a 1-row numpy array.

    MODIFICATIONS FROM JEREMY (21/04/2017)
    for initial_state_choice == 'specified'
    retrieves the initial_state_genes variable, which contains the IDs of the genes that we want set to 1
    returns a state variable in accordance with the considered gene

    """

    if initial_state_choice == 'random':
        state = np.random.random_integers(0, 1, len(genes_names))

    elif initial_state_choice == 'all_zeros':
        state = np.zeros(len(genes_names), dtype=np.int)

    elif initial_state_choice == 'all_ones':
        state = np.ones(len(genes_names), dtype=np.int)

    elif initial_state_choice == 'specified':
        if initial_state_genes == ['foo']:
            print("\nList of active initial genes is empty - ending computation")
            return False

        else:
            count = 0
            res = [0] * len(genes_names)
            for x in genes_names[:]:
                if x in initial_state_genes:
                    res[count] = 1
                    count += 1
                else:
                    count += 1
            state = array(res)

    elif initial_state_choice == 'random-specified':
        a = np.random.random_integers(0, 1, len(genes_names))
        state = np.array(a)
        all_ones = np.ones(len(genes_names), dtype=np.int)
        if initial_state_genes == ['foo']:
            print("\nGene non specified - ending computation")
            return False
        else:
            count = 0
            res = [0] * len(genes_names)
            for x in genes_names[:]:
                if x in initial_state_genes:
                    res[count] = 1
                    count += 1
                else:
                    res[count] = rd.randint(0, 1)
                    count += 1
            state = array(res)

    else:
        print("\nInvalid initialization choice - ending computation")
        return False

    return state


def ComputeNextState(genes_names, genes_network, state_flow, model='logical', stimulus='transient', KO_genes=['foo'],
                     OA_genes=['foo']):
    """
    Compute the next state and append it to the state flow according to the chosen model and given gene network.

    genes_network : genes network in the format { source gene name : [interaction, target gene name] },
    state_flow : numpy array with one column per gene in the gene network - lines are the successive states, last state is a copy of next-to-last state that will be updated
    model : choice of model to run for the computation of the next state, default is 'logical' - pure logical without memory, other options are 'algebraic'
    stimulus : choice of model for the behavior of network boundary root nodes, default is 'transient' - switch to off after being active, other option is 'constant' - will stay active indefinitely
    KO_genes : list of genes to be inactivated for the computation (always at 0)
    OA_genes : list of genes to be overactivated for the computation (always at 1)


    Return the updated numpy array state_flow.
    """

    # note : keep track of gene name and respective state using genes_names.index('name_of_gene_of_interest') call
    # note : check array dimensions with array_name.ndim
    # note : use map(function,sequence) to apply a function to all elements of the given sequence

    # logical memoryless model -
    # activator present --> activation
    # activator absent --> do not activate (can be overriden by the presence of other activators)
    # repressor present --> repression regardless of number of activators presents
    # repressor absent : either there are other activators, then they have control of target gene activity,
    #					 OR, there is no other positive source, then it is assumed that the target gene is auto-activating
    # root genes : if stimulus is set to transient, root genes will switch off after the fisrt time step, if stimulus is set to constant, root genes state will not change during the computation

    # logical algebraic model -
    # make the sum of activator and repressor signal, if the end result is positive, the regulated gene is switched on, otherwise it is switched off
    # root genes : if stimulus is set to transient, root genes will switch off after the fisrt time step, if stimulus is set to constant, root genes state will not change during the computation

    if model == 'logical':
        # start by blanking the copied last state - no memory in the pure logical model
        state_flow[-1] = np.zeros(len(state_flow[-1]))

        # list all genes that receive input from the network to use as a reference to identify 'root' genes (i.e. genes with no input)
        all_targets = genes_network.values()
        regulated_genes = []
        for target_list in all_targets:
            for target in target_list:
                regulated_genes.append(target[1])
        regulated_genes = list(set(regulated_genes))

        for source_gene in genes_network:

            if source_gene not in regulated_genes:
                # root genes depend on input external to the network, we use the specified model to determine their behavior when we encounter them since they are not otherwise regulated
                # transient - gene is only stimulated once
                # constant - any stimulation is constant
                if stimulus == 'transient':
                    state_flow[-1][genes_names.index(source_gene)] = 0  # stimulus only applied once if applied
                elif stimulus == 'constant':
                    state_flow[-1][genes_names.index(source_gene)] = state_flow[-2][genes_names.index(
                        source_gene)]  # could have left the line blank as no other input will change the root genes
                else:
                    print(
                        "\nWarning - wrong model name given for the behavior of network root genes - using default model (transient) instead")
                    state_flow[-1][genes_names.index(source_gene)] = 0  # stimulus only applied once if applied

            # and then we apply the proper regulations to their targets

            for [interaction, target_gene] in genes_network[source_gene]:
                # in pure logic, a single False value trumps all True values, no use to continue computing if gene has been inactivated
                if state_flow[-1][genes_names.index(
                        target_gene)] >= 0:  # state >= 0 means that no inhibition flag has been raised yet for this gene
                    if float(interaction) < 0:  # source gene is an inhibitor
                        if state_flow[-2][genes_names.index(source_gene)]:  # inhibitor is present
                            state_flow[-1][genes_names.index(
                                target_gene)] = -1  # flag target with inescapable inactivation for the rest of this network state update
                        else:  # inhibitor is absent
                            if state_flow[-1][genes_names.index(
                                    target_gene)] == 0:  # if state is 0, no positive regulator (active or inactive) has been encountered yet in the network exploration
                                state_flow[-1][genes_names.index(
                                    target_gene)] = 2  # flag for autoactivation should no further positive regulator been found during the network exploration
                    elif float(interaction) > 0:  # source gene is an activator
                        if state_flow[-2][genes_names.index(source_gene)]:  # activator is present
                            state_flow[-1][genes_names.index(
                                target_gene)] = 1  # activate the target gene, override existing flags
                        else:  # activator absent - raise non-autoactivation flag if gene is not already activated because an activator source exist for target gene, dropping any existing auto-activation flag
                            if (state_flow[-1][genes_names.index(target_gene)] == 0) or (
                                    state_flow[-1][genes_names.index(target_gene)] == 2):
                                state_flow[-1][genes_names.index(target_gene)] = 3
                    else:  # null interaction - should not happen
                        print("\nWarning - null interaction detected between "), source_gene, target_gene

        # set KO or overactivated genes to their constant values - do it in a single pass rather than checking during each previous computation
        for gene in genes_names:
            if gene in KO_genes:
                state_flow[-1][genes_names.index(gene)] = 0  # always off
            elif gene in OA_genes:
                state_flow[-1][genes_names.index(gene)] = 1  # always on

        # binarize gene activity to remove repression, autoactivation and non-autoactivation flags after all operations are done
        # binarize = lambda x: (0 if (x < 0) or (x == 3) else 1 if x > 1 else x)
        # state_flow[-1] = map(binarize, state_flow[-1])
        state_flow[-1] = [(0 if (x < 0) or (x == 3) else 1 if x > 1 else x) for x in state_flow[-1]]



    # logical algebraic model - gene activity is dependant on the sum of the activity of its regulators
    elif model == 'algebraic':

        # start by blanking the copied last state - no memory in the pure algebraic model
        state_flow[-1] = np.zeros(len(state_flow[-1]))

        # list all genes that receive input from the network to use as a reference to identify 'root' genes (i.e. genes with no input)
        all_targets = genes_network.values()
        regulated_genes = []
        for target_list in all_targets:
            for target in target_list:
                regulated_genes.append(target[1])
        regulated_genes = list(set(regulated_genes))

        for source_gene in genes_network:

            if source_gene not in regulated_genes:
                # root genes depend on input external to the network, we use the specified model to determine their behavior once we encounter them since they are not otherwise regulated
                # transient - gene is only stimulated once
                # constant - any stimulation is constant
                if stimulus == 'transient':
                    state_flow[-1][genes_names.index(source_gene)] = 0  # stimulus only applied once if applied
                elif stimulus == 'constant':
                    state_flow[-1][genes_names.index(source_gene)] = state_flow[-2][genes_names.index(
                        source_gene)]  # could have left the line blank as no other input will change the root genes
                else:
                    print(
                        "\nWarning - wrong model name given for the behavior of network root genes - using default model (transient) instead")
                    state_flow[-1][genes_names.index(source_gene)] = 0  # stimulus only applied once if applied

            # then we apply the proper regulation to their targets

            for [interaction, target_gene] in genes_network[source_gene]:
                state_flow[-1][genes_names.index(target_gene)] += (
                            float(interaction) * float(state_flow[-2][genes_names.index(source_gene)]))

        # set KO or overactivated genes to their constant values - do it in a single pass rather than checking during each previous computation
        for gene in genes_names:
            if gene in KO_genes:
                state_flow[-1][genes_names.index(gene)] = 0  # always off
            elif gene in OA_genes:
                state_flow[-1][genes_names.index(gene)] = 1  # always on

        # binarize gene activity after all operations are done
        # binarize = lambda x: (0 if x < 0 else 1 if x > 1 else x)
        # state_flow[-1] = map(binarize, state_flow[-1])

        state_flow[-1] = [(0 if x < 0 else 1 if x > 1 else x) for x in state_flow[-1]]


    # additive inertial model - target gene activity at time (t-1) in the copied state is added to (interaction*source gene(t-1) in original state)
    elif model == 'additive_inertial':
        for source_gene in genes_network:
            for [interaction, target_gene] in genes_network[source_gene]:
                state_flow[-1][genes_names.index(target_gene)] += (
                        float(interaction) * float(state_flow[-2][genes_names.index(source_gene)]))
        # binarize gene activity after all operations are done
        binarize = lambda x: (0 if x < 0 else 1 if x > 1 else x)
        state_flow[-1] = map(binarize, state_flow[-1])

    else:
        print("\nInvalid model name - ending computation")
        return False

    return state_flow


def HarvestStableStates(genes_names, genes_network, state_flow_network, verbose=True):
    """
    Take a gene list, a gene network and its corresponding state flow network, extract the stable states and analyze them.

    genes_names: list of names for the genes in the network
    genes_network: genes network in the format { source gene name : [interaction, target gene name] },
    state_flow_network: state flow in the format { initial state : next state } where a state is a list of 0 and 1 (state of each gene in the network)

    Return a dictionnary of {stable state : attraction basin size}.
    """

    # create the dictionnary of stable state
    stable_states = {}

    # create a list of all possible starting states to explore
    # rationalize the harvest process by only starting from the origin states in the flow, i.e. the states that are present in the keys but not in the values

    if len(state_flow_network) == 0:
        print('\n invalid state flow network (empty) - ending computation')
        return false

    elif len(state_flow_network) == 1:  # only a single stable state in the flow
        for k in state_flow_network.keys():
            single = k
        stable_states[tuple(single)] = 1

    else:  # more than 1 state in the flow
        starts = set(state_flow_network.keys())
        ends = set(state_flow_network.values())
        starting_points = list(starts - ends)

        if starting_points == []:  # the state flow network was only composed of the states forming a cyclic dynamic stable state
            state = rd.choice(state_flow_network.keys())
            explored_states = []
            path = []
            count = 0
            total_state_nb = len(state_flow_network)

            while (state != state_flow_network[state]) and (state not in explored_states):
                path.append(state)
                explored_states.append(state)
                count += 1
                if verbose:
                    if (total_state_nb > 10):
                        if not (count % tenpercent):
                            print("Harvesting stable states - %.0f%% of all given states explored") % (
                                    100 * float(count) / float(total_state_nb))
                    else:
                        print("Harvesting stable state from state", count, "out of", total_state_nb)

                state = state_flow_network[state]

                if (state in path):  # we have come back to the start of the loop
                    start = path.index(
                        state)  # find starting point of this new loop - index returns the position of the first occurence of 'state' in the path
                    loop = path[start:]  # slice the part of the path that correspond to the new loop
                    stable_states[tuple(loop)] = len(
                        path)  # create an entry for the loop in the stable states record with initial basin size = len path (states leading to the loop + size of the loop)


        else:  # flow start outside of any loop
            explored_states = []

            #####
            # BEWARE : need to create a copy of the starting points list with slicing [:] or list(starting_points) because removing elements from this list while iterating on it
            # generate strange behaviour where not all elements are dealt with
            # NB : we do not remove element for the list anyway ; the slicing for keeps all element in mmemory from the first slicing and do not update its list... it is
            # too dangerous to use it like this --> we just check whether an element has already been encountered
            #####

            count = 0
            total_state_nb = len(state_flow_network)
            tenpercent = total_state_nb / 10

            for state in starting_points[:]:
                # for each possible starting state left, initialize a new path of attraction to record potential loops and basin size count
                path = []
                # if the current state is not stable and has not been previously encountered
                # then increase the path size, record the passage through this state and progress to the next state in the state flow
                while (state != state_flow_network[state]) and (state not in explored_states):
                    path.append(state)
                    explored_states.append(state)
                    count += 1
                    if verbose:
                        if (total_state_nb > 10):
                            if not (count % tenpercent):
                                print("Harvesting stable states -", 100 * count / total_state_nb,
                                      "of all given states explored")
                        else:
                            print("Harvesting stable state from state", count, "out of", total_state_nb)
                    # remove this new state from the possible starting points if it is there to avoid charting the same path multiple times -
                    # if state in starting_points:
                    #	starting_points.remove(state)
                    state = state_flow_network[state]

                # either we started from an already accounted for state, or the while loop has ended after at least one iteration
                # if the while loop iterated, either the last encountered state is part of a dynamic stable state loop, or a stable state,
                # or just a junction point to an already explored path (n.b. the already explored path can be a single point - i.e. a stable state)

                if (state in path):  # the path contain a dynamical stable state loop that was not already accounted for
                    start = path.index(
                        state)  # find starting point of this new loop - index returns the position of the first occurence of 'state' in the path
                    loop = path[start:]  # slice the part of the path that correspond to the new loop
                    stable_states[tuple(loop)] = len(
                        path)  # create an entry for the loop in the stable states record with initial basin size = len path (states leading to the loop + size of the loop)

                elif (
                        state in explored_states):  # either started from an already treated state, or found a junction point - means that what is forward of there has already been dealt with - find the downward stable point and add to its basin size the len of newly uncovered path
                    if path != []:  # we did at least one iteration of the while loop
                        terminal_flag = False
                        while not terminal_flag:  # while we have not reached one of the stable states
                            for final_state_list in stable_states.keys():  # check all stable states one by one
                                if not any(isinstance(sub_state, tuple) for sub_state in
                                           final_state_list):  # there is only a unique stable state, not a loop with sub-states
                                    if (state == final_state_list):
                                        terminal_flag = True
                                        stable_states[final_state_list] += len(
                                            path)  # add the len of path to the basin of attraction size for this stable state
                                        break  # stop searching through the stable states
                                elif state in final_state_list:  # if the current state is part of a stable loop
                                    terminal_flag = True
                                    stable_states[final_state_list] += len(
                                        path)  # add the len of path to the basin of attraction size for this stable state
                                    break  # stop searching through the stable states
                            if not terminal_flag:  # if the current state was not part of a stable state, progress one state further
                                state = state_flow_network[state]
                    else:  # there was no iteration of the while loop, we started right from an already explored state - do nothing
                        pass

                elif (state == state_flow_network[state]):  # unique stable state not already encountered
                    explored_states.append(state)  # check it as explored now
                    # if state in starting_points:     		# remove it from possible starting points
                    #	starting_points.remove(state)
                    stable_states[tuple(state)] = len(
                        path) + 1  # create an entry for the stable state as a tuple with 1 element and basin size = len path leading there + 1 for the stable state
                    # we found a new attractor which jumped us out of the while loop without counting it - so count it for the progress tracker
                    count += 1
                    if verbose:
                        if (total_state_nb > 10):
                            if not (count % tenpercent):
                                print("Harvesting stable states - %.0f%% of all given states explored" % (
                                        100 * float(count) / float(total_state_nb)))
                        else:
                            print("Harvesting stable state from state", count, "out of", total_state_nb)
                else:
                    print(
                        'Something strange is happening that should not happen during harvesting. Ending computation.')
                    return False

    return stable_states


def RunBooleanModel(genes_names, genes_network, initial_state_number='all', initial_state_choice='random',
                    initial_state_genes=['foo'],
                    model='logical', stimulus='transient', verbose=True, KO_genes=['foo'], OA_genes=['foo']):
    """
    Run the chosen boolean model over the given network.

    genes_names: list of names for the genes in the network
    genes_network: genes network in the format { source gene name : [interaction, target gene name] },
    output_stateflow_filename: name of the file to write the state flow output in
    output_stablestates_filename : name of the file to write the stable states output in

    initial_state_number: 'all' (default) / 'xxx' (int) - either run the model from all possible states or a given xxx number of states
    initial_state_choice: 'random' (default) / 'all_zeros' / 'all_ones' / 'specific' / 'random-specific' - specific set a given set of genes to 1 and the rest to 0, random-specific set a given set of genes to 1 and the rest randomly to 0 or 1
    model : 'logical' (default) / 'algebraic' / ... - choice of model to run for the computation of the next state, pure logical without memory by default
    stimulus : choice of model for the behavior of network root nodes, default is 'transient' - switch to off after being active once, other option is 'constant' - will stay active indefinitely

    verbose : boolean - determine if information regarding computation progress are displayed or not

    specified_starting_states : set of pre-defined network states to be used as starting points for the computation
    KO_genes : list of genes to be inactivated for the computation (always at 0)
    OA_genes : list of genes to be overactivated for the computation (always at 1)

    depending on the model choice, interaction is either :
    - a number quantifying the degree of positive or negative regulation of source over target
    - ...

    Return a dictionary of all the state flow for the given initial state, a list of stable states sequences + size of attraction basin and either the starting state (run from a single state) or 'all_state_run' if running from all possible states
    """

    # note : keep track of gene name and respective state using genes_names.index('name_of_gene_of_interest') call
    # note : use np.vstack((array1,array2,...)) to add new state in the state flow
    # note : can use zero_row = np.zeros(len(genes_names)) to create a line of 0

    # Create the collection of starting states according to the user choice
    if verbose:
        print("\nInitializing network starting state collection")
    start_time = time.time()

    initial_states = []

    if initial_state_number == 'all':
        print(
            "\nWarning - if the number of genes is large, generating all possible network states can be an extensive process")
        # create all possible 01 configurations for the given genes names space
        all_states_sequences = map(''.join, itertools.product('01', repeat=len(genes_names)))
        # convert those string sequence to numpy arrays with all gene states individualized
        for seq in all_states_sequences:
            initial_states.append(np.array(map(int, seq)))
        print("\nWarning - the total number of possible initial states is", len(initial_states),
              tenpercent=len(initial_states) / 10)
    else:  # generate a given number of initial states
        try:
            number = int(initial_state_number)
        except:
            print("\nInvalid number of initial states - please enter either 'single', 'all' or a valid integer")
        else:
            if number <= 0:
                print("\nNull or negative number of starting states - adjusting computation to a single starting state")
                number = 1
            if number == 1:
                starting_state = InitializeState(genes_names, initial_state_choice, initial_state_genes)
                initial_states.append(starting_state)
            else:
                for i in range(number):  # using xrange instead of range to speed up allocating process
                    initial_states.append(InitializeState(genes_names, initial_state_choice, initial_state_genes))
                tenpercent = len(initial_states) / 10

    # create a flag for stable state reaching (or already encountered state reaching), global records of the state flow network and stable states
    stop_flag = False
    state_flow_network = {}
    stable_states = {}

    # start the loop for flow progression from each non already checked starting state
    if verbose:
        print("\nComputing state progression from each starting state")

    count = 0

    # BEWARE : as it is dangerous to edit a list you are looping over, create a slicing copy for the loop to iterate over and edit the original
    for start in initial_states[:]:

        # keep user informed of the computation progress
        count += 1
        if verbose:
            if (initial_state_number != 'single') and (len(initial_states) > 10):
                if not (count % tenpercent):
                    print("Exploring - %.0f%% of specified initial states done" % (
                            100 * float(count) / float(len(initial_states))))
            else:
                print("Exploring from state", count, "out of", len(initial_states))

        # if this state has not already been encountered previously in the flow, follow this path for a first iteration
        # otherwise skip to next possible initial state
        if not (tuple(start) in state_flow_network):

            # reference the vstack method once outside all loops to speed up computation
            # testing show a 30% gain on computation time for large number of initial states that way
            stack = np.vstack

            # Initialize the new flow
            state_flow = start
            # add a copy of the initial state to the local state flow for first computation of the next state
            state_flow = stack((state_flow, state_flow))

            # make the first state computation outside of the while loop to avoid having to loop over
            # an additional 'if' choice for each of the first state copy instruction (can't use vstack w/state_flow[-1] if there is only one row)

            # compute the new state
            state_flow = ComputeNextState(genes_names, genes_network, state_flow, model, stimulus, KO_genes, OA_genes)

            # add the successive states in the global state flow network dictionary as {state(t-1) -> state(t)} entry
            state_flow_network[tuple(state_flow[-2])] = tuple(state_flow[-1])

            # check if the last generated state has not already been encountered (i.e. is present in state_flow_network keys)
            stop_flag = tuple(state_flow[-1]) in state_flow_network.keys()

            # generate new states until an already encountered state is found (i.e. a stable state, a new loop or an already encountered branch)

            while not stop_flag:
                # add a copy of the last state to the local state flow as a basis for the computation of the next state
                # if model is 'logical', this added state will be wiped (no memory)
                # in other model cases, it may be used as the memory of network state for computation purpose
                state_flow = stack((state_flow, state_flow[-1]))
                state_flow = ComputeNextState(genes_names, genes_network, state_flow, model, stimulus, KO_genes,
                                              OA_genes)
                state_flow_network[tuple(state_flow[-2])] = tuple(state_flow[-1])
                stop_flag = tuple(state_flow[-1]) in state_flow_network.keys()

    if verbose:
        print("\nHarvesting stable states")

    stable_states = HarvestStableStates(genes_names, genes_network, state_flow_network, verbose)

    # create the (stable state size, basin size) list - have to use isinstance check to avoid unique stable state size being estimated as the number of genes in them
    data = []
    for k in stable_states.keys():
        if not any(isinstance(sub_state, tuple) for sub_state in k):
            data.append((1, stable_states[k]))
        else:
            data.append((len(k), stable_states[k]))

    time_passed = (time.time() - start_time)
    if verbose:
        print("\n####### Results summary #######\n#\n# Run parameters :", len(genes_names), "genes in the network ;",
              len(
                  Flatten(list(genes_network.values()))) / 2, "interactions in the network ;\n# ", model, "model ;",
              stimulus, "boundary condition ; ", "KO genes : ", KO_genes, " ; OA genes : ", OA_genes)
        print("# Number of initial states :", len(initial_states))
        print("# Number of stable states :", len(stable_states))
        print("# Size of stable states and respective attraction basins :")
        print("#", data)
        print("# Running time (seconds) %s", time_passed)
        print("#\n########################")
        print("\n")
        print("# Full run result saved to csv files.")

    if (initial_state_number == 'single') or (int(initial_state_number) <= 1):
        return state_flow_network, stable_states, starting_state, len(stable_states), len(initial_states), \
               len(genes_names), data, len(Flatten(list(genes_network.values()))) / 2, time_passed
    else:
        return state_flow_network, stable_states, 'several_state_run', len(stable_states), len(initial_states), \
               len(genes_names), data, len(Flatten(list(genes_network.values()))) / 2, time_passed
