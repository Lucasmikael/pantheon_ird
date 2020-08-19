import matplotlib.pyplot as plt


#############################################################
# 		network manipulation and analyses tools 		    #
#############################################################

def ExtractCoreNetwork(genes_network):
    """
    Takes a genes network and recursively prune it from the leaf up and root down until only the core network remain (i.e. cluster or set of clusters of genes).
    Returns the core network.
    """

    # to pop a key from a dictionary and return its value or 'default_value' if it is not in the dict : name_of__dict.pop(name_of_key, default_value)
    # to check if a key is in a dictionary and return its value of 'default_value' if it is not in the dict : name_of_dict.get(name_of_key, default_value)
    # BEWARE : you cannot change the size of a dictionary while iterating over it -> first create a list of keys to be removed, then iterate over this list to delete these entries from the dictionary
    # OR create a copy of the dict to iterate over while you delete the entries from the original dictionary name_of_dict2 = dict(name_of_dict)
    # to iterate efficiently over a dict key/value pairs, use .items() : for key, values in name_of_dict.items(): whatever...

    genes_network_copy = dict(genes_network)
    source_genes = genes_network.keys()
    target_genes = []
    for inter_target_list in genes_network.values():
        for inter_target in inter_target_list:
            target_genes.append(inter_target[1])
    target_genes = set(target_genes)

    network_downcopy = dict(genes_network_copy)
    down_pruning_not_done = 1

    # pruning the network down the roots first
    print ('Extracting core network - root pruning in progress...')
    while down_pruning_not_done :
        # set the flag to stop if no genes are removed from the network during the iteration
        down_pruning_not_done = 0
        for source, targets in network_downcopy.items():
            if source not in target_genes:
                # remove all genes with no regulators from the network
                genes_network_copy.pop(source,None)
                # update the list of genes being targeted
                target_genes = []
                for inter_target_list in genes_network_copy.values():
                    for inter_target in inter_target_list:
                        target_genes.append(inter_target[1])
                target_genes = set(target_genes)
                # raise the flag for an additional pass each time we find genes to remove
                down_pruning_not_done = 1
        # update the network copy for the next pass if there is one
        if down_pruning_not_done:
            network_downcopy = dict(genes_network_copy)


    network_upcopy = dict(genes_network_copy)
    up_pruning_not_done = 1

    #pruning the network up the leaves next
    print ('Extracting core network - leaf pruning in progress...')
    while up_pruning_not_done :
        # set the flag to stop if no genes are removed from the network during the iteration
        up_pruning_not_done = 0
        for source, targets in network_upcopy.items():
            target_list = list(targets)
            for current_target in target_list:
                # current_target is in the format [interaction, target_gene_name]
                if current_target[1] not in network_upcopy.keys():
                    # remove target genes which have no target themselves
                    target_list.remove(current_target)
                    # update the list of targets from the network and if it was the last target of the source gene, remove its entry from the network
                    if target_list != []:
                        genes_network_copy[source]=target_list
                    else:
                        genes_network_copy.pop(source,None)
                    # raise the flag for an additional pass each time we find genes to remove
                    up_pruning_not_done = 1
        # update the network copy for the next pass if there is one
        if up_pruning_not_done:
            network_upcopy = dict(genes_network_copy)

    return genes_network_copy


def Flatten(l, ltypes=(list, tuple)):
    """
    Utility function to flatten a nested list.
    Take a list and return a list.
    """
    ltype = type(l)
    l = list(l)
    i = 0
    while i < len(l):
        while isinstance(l[i], ltypes):
            if not l[i]:
                l.pop(i)
                i -= 1
                break
            else:
                l[i:i + 1] = l[i]
        i += 1
    return ltype(l)

### check this to replace Flatten
#>>> a = [[1, 2], [3, 4], [5, 6]]
#>>> list(itertools.chain.from_iterable(a))
#[1, 2, 3, 4, 5, 6]

def IterPlot(genes_names, genes_network, nb_initial_states, repeat, stimulus='transient'):
    """
    Iterate the Boolean model nb_initial_states times and plot the graph 'basin size'=f('stable state size')

    genes_names: list of names for the genes in the network
    genes_network: genes network in the format { source gene name : [interaction, target gene name] },
    nb_initial_states: number of random initial states to be tested at each iteration
    repeat: number of iteration
    stimulus: model for behavior of entry nodes to the network - either 'constant' or 'transient'. Default is 'transient'.
    """


    data_to_plot = []
    max_x = 0
    max_y = 0
    plt.figure(1)
    plt.subplot(111)  # Note : subplot(xyz) == subplot(x,y,z) if x*y <10 - x:nb row, y:nb col, z: fig id from 1 to x*y
    for i in xrange(repeat):
        print ("\nIteration"), i + 1
        flow, stable = RunBooleanModel(genes_names, genes_network, 't', 't', nb_initial_states, stimulus)
        x = [len(k) for k in stable.keys()]
        y = [stable[k] for k in stable.keys()]
        max_x = max(max_x, x)
        max_y = max(max_y, y)
        plt.scatter(x, y, c=np.random.rand(nb_initial_states), alpha=0.5)
    plt.ylabel('Size of attraction basin')
    plt.xlabel('Size of stable state')
    plt.axis([0, max(x) + 5, 0, max(y) + 5])
    plt.show()


def AccessGene(gene, genes_names, network_state, edit=False, value=0):
    """
    Access the state of given gene in a given global network states or an ensemble of state (e.g. for cyclical stable states).
    If a single network state is given and edit is True, will edit the gene state in network_state to 'value'

    gene: the name of the gene of interest as it would appear in the genes_names listing
    genes_names: genes names listing generated by the ImportBooleanModel function and used for RunBooleanModel call
    network_state: a single gene network state

    Return either the single gene state, a list of the gene state or the edited global network state
    """

    # reading gene state in the given set of global network states
    if not edit:
        if not any(isinstance(sub_state, tuple) for sub_state in network_state):  # single network state
            try:
                state = network_state[genes_names.index(gene)]
                return state
            except:
                print ("\nInvalid gene name")
        else:  # multiple network state
            state = []
            try:
                for sub_state in network_state:
                    state.append(sub_state[genes_names.index(gene)])
                return state
            except:
                print ("\nInvalid gene name")
    # editing gene state in the given global unique network state
    else:
        try:
            if not any(isinstance(sub_state, tuple) for sub_state in network_state):  # single network state
                net = list(network_state)  # dict key must be tuple - have to make it a list for edition
                net[(genes_names.index(gene))] = value
                return tuple(net)
            else:
                print ("Cannot edit multiple network state simultaneously, please enter a single network state for edition.")
        except:
            print ("\nInvalid gene name")


def CheckGeneStableStates(gene, genes_names, stable_states):
    """
    Check the state of given gene in a dictionnary of stable network states {state:basin size}.

    gene: the name of the gene of interest as it would appear in the genes_names listing
    genes_names: genes names listing generated by the ImportBooleanModel function and used for RunBooleanModel call
    stable_states: set of gene network global states

    Return a list of states for given gene.
    """

    result = []
    for states in stable_states.keys():
        result.append(AccessGene(gene, genes_names, states))
    return result


def FilterGenes(genes_names, states_set, criterion='changing'):
    """
    Search the gene list and set of states to find all genes that fit the filter.

    genes_names: genes names listing generated by the ImportBooleanModel function and used for RunBooleanModel call
    states_set: set of gene network global states
    criterion: 'changing', 0 or 1 - select genes that change states, are always off (0) or always on (1)

    Return a set of lists of gene names (1 list per list of states analyzed).
    """

    filtered_genes = []

    if criterion == 'changing':
        # if there is only one state in the given set, cannot search for changing gene activity
        if not any(isinstance(sub_state, tuple) for sub_state in states_set):  # single network state
            print ("\nSingle network state - cannot filter for gene activity change")

        else:  # multiple network state - can search for changes in gene activity
            init_state = states_set[0]
            if len(init_state) != len(genes_names):
                print ("\nWarning - gene name list length is not equal to the number of gene states - returning empty list as N/A result")
                return filtered_genes
            for gene in genes_names:  # for each gene
                init_gene = init_state[genes_names.index(gene)]  # get the first state
                for sub_state in states_set:  # go through the states while the gene state does not change
                    if sub_state[genes_names.index(
                            gene)] != init_gene:  # when gene state changes, record gene name and go to next gene
                        filtered_genes.append(gene)
                        break

    elif criterion == 0:
        if not any(isinstance(sub_state, tuple) for sub_state in states_set):  # single network state
            if len(states_set) != len(genes_names):
                print ("\nWarning - gene name list length is not equal to the number of gene states - returning empty list as N/A result")
                return filtered_genes
            for gene in genes_names:  # for each gene
                init_gene = states_set[genes_names.index(gene)]  # get the gene state
                if (init_gene == 0):  # if it is 0
                    filtered_genes.append(gene)  # record the gene name

        else:  # multiple network state - can search for changes in gene activity
            init_state = states_set[0]
            if len(init_state) != len(genes_names):
                print ("\nWarning - gene name list length is not equal to the number of gene states - returning empty list as N/A result")
                return filtered_genes
            for gene in genes_names:  # for each gene
                init_gene = init_state[genes_names.index(gene)]  # get the first state
                if (init_gene == 0):  # if its 0
                    for sub_state in states_set:  # go through the states while the gene state does not change
                        if sub_state[
                            genes_names.index(gene)] != init_gene:  # as soon as it changes, stop and go to next gene
                            break

    elif criterion == 1:
        if not any(isinstance(sub_state, tuple) for sub_state in states_set):  # single network state
            if len(states_set) != len(genes_names):
                print ("\nWarning - gene name list length is not equal to the number of gene states - returning empty list as N/A result")
                return filtered_genes
            for gene in genes_names:  # for each gene
                init_gene = states_set[genes_names.index(gene)]  # get the gene state
                if (init_gene == 1):  # if it is 1
                    filtered_genes.append(gene)  # record the gene name

        else:  # multiple network state - can search for changes in gene activity
            init_state = states_set[0]
            if len(init_state) != len(genes_names):
                print ("\nWarning - gene name list length is not equal to the number of gene states - returning empty list as N/A result")
                return filtered_genes
            for gene in genes_names:  # for each gene
                init_gene = init_state[genes_names.index(gene)]  # get the first state
                if (init_gene == 1):  # if its 1
                    for sub_state in states_set:  # go through the states while the gene state does not change
                        if sub_state[
                            genes_names.index(gene)] != init_gene:  # as soon as it changes, stop and go to next gene
                            break
    else:
        print ("\nInvalid filter criterion")

    return filtered_genes


def FilterStableStates(genes_names, stable_states, criterion='changing'):
    """
    Apply the selected filter to the dictionnary of stable state to find all genes that fit.

    genes_names: genes names listing generated by the ImportBooleanModel function and used for RunBooleanModel call
    stable_states: set of gene network global states
    criterion: 'changing', 0 or 1 - select genes that change states, are always off (0) or always on (1)

    Return a set of list of gene names (1 list per stable state analyzed).
    """

    result = []
    for states in stable_states.keys():
        result.append(FilterGenes(genes_names, states, criterion))
    return result

def generate_pairs(source):
    pairs = []
    for p1 in range(len(source)):
        for p2 in range(p1+1,len(source)):
            pairs.append([source[p1],source[p2]])
    return pairs
