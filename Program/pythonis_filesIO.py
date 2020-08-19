import csv
import datetime


def WriteOutputToCsv(filename, genes_list, stable_states_list):
    """
    Write the result of a run in csv file format. First line is the genes names list, lines afterward are the sequence of stable states. The line immediately after a stable state flow hold the size of its attraction basin.
    :param filename:
    :param gene_list:
    :param stable_states:
    :return:
    """
    run_result_file = open(filename, 'w')

    # quick an dirty data drop in a CSV compliant format to a text file
    for gene in genes_list[:-1]:  # write all names followed by ',' except the last one
        run_result_file.write(gene + ',')
    run_result_file.write(genes_list[-1] + "\n")  # write the last one and go to next line

    for stable_state, basin_size in stable_states_list.items():
        # if stable state is not a nested list then it is a single state and will be written on a single line
        if not any(isinstance(i, tuple) for i in stable_state):
            strip_state = str(stable_state).replace('(', '')
            strip_state = strip_state.replace(')', '')
            run_result_file.write(strip_state + '\n')
            run_result_file.write(str(basin_size) + '\n')
        else:
            for state in stable_state:
                strip_state = str(state).replace('(', '')
                strip_state = strip_state.replace(')', '')
                run_result_file.write(strip_state + '\n')
            run_result_file.write(str(basin_size) + '\n')

    run_result_file.close()


def WriteFlowToCsv(filename, genes_list, state_flow, stable_state):
    """
    Write the state flow from a single initial state run in csv file format. First line is the genes names list, lines afterward are the sequence of state until the stable state. Last line hold the size of the stable state.
    """
    run_result_file = open(filename, 'w')

    # quick an dirty data drop in a CSV compliant format to a text file
    for gene in genes_list[:-1]:  # write all names followed by ',' expect the last one
        run_result_file.write(gene + ',')
    run_result_file.write(genes_list[-1] + "\n")  # write the last one and go to next line

    # find and write the starting state
    if len(state_flow) == 1:  # flow with a single starting and ending state
        for k in state_flow.keys():
            lone_state = k
        strip_state = str(lone_state).replace('(', '')
        strip_state = strip_state.replace(')', '')
        strip_state = strip_state.replace('set', '')
        strip_state = strip_state.replace('[', '')
        strip_state = strip_state.replace(']', '')
        run_result_file.write(strip_state + '\n' + '1')
        run_result_file.close()


    else:  # more than 1 state in the flow
        starts = set(state_flow.keys())
        ends = set(state_flow.values())
        lone_state = starts - ends
        strip_state = str(lone_state).replace('(', '')
        strip_state = strip_state.replace(')', '')
        strip_state = strip_state.replace('set', '')
        strip_state = strip_state.replace('[', '')
        strip_state = strip_state.replace(']', '')
        run_result_file.write(strip_state + '\n')

        # write the flow
        start = lone_state.pop()
        current_state = start
        encountered_states = []
        encountered_states.append(current_state)

        while current_state in starts and state_flow[current_state] not in encountered_states:
            next_state = state_flow[current_state]
            strip_state = str(next_state).replace('(', '')
            strip_state = strip_state.replace(')', '')
            run_result_file.write(strip_state + '\n')
            current_state = next_state
            encountered_states.append(current_state)

        # add the size of the stable state for reference
        stable = stable_state.keys()
        for i in stable:
            run_result_file.write(str(len(i)) + '\n')

        run_result_file.close()


def WriteDistanceToCsv(filename, genes_list, meansKO, meanOA, sample_size):
    """
    Write the result of a run in csv file format. First line is the genes names list, lines afterward are the sequence of stable states. The line immediately after a stable state flow hold the size of its attraction basin.
    :param filename:
    :param gene_list:
    :param stable_states:
    :return:
    """

    # time stamp for the current run
    run_ref = str(datetime.datetime.now())
    ref_for_filename = run_ref.replace(':', '.')
    ref_for_filename = ref_for_filename.replace(' ', '.')
    csv_file = filename + '_' + ref_for_filename + '_' + str(sample_size) + 'IS.csv'

    run_result_file = open(csv_file, 'w')

    # quick data drop in a CSV compliant format to a text file
    run_result_file.write('Genes,mean distance for KO,mean distance for OA\n')
    for gene in genes_list:  # write all names followed by ',' and the values of their mean distances for KO and OA mutations
        run_result_file.write(gene + ',' + str(meansKO[gene]) + ',' + str(meanOA[gene]) + '\n')

    run_result_file.close()


def WriteTargetedKODistanceToCsv(filename, KO_genes_list, meansKO, sample_size):
    """
    Write the result of a run in csv file format. First line is the genes names list, lines afterward are the sequence of stable states. The line immediately after a stable state flow hold the size of its attraction basin.
    :param filename:
    :param gene_list:
    :param stable_states:
    :return:
    """

    # time stamp for the current run
    run_ref = str(datetime.datetime.now())
    ref_for_filename = run_ref.replace(':', '.')
    ref_for_filename = ref_for_filename.replace(' ', '.')
    csv_file = filename + '_' + ref_for_filename + '_' + str(sample_size) + 'IS.csv'

    run_result_file = open(csv_file, 'w')

    # quick data drop in a CSV compliant format to a text file
    run_result_file.write('Genes,mean distance for KO\n')
    for KO_set in KO_genes_list:  # write all KO set names followed by ',' and the values of their mean distances for the mutations
        run_result_file.write(str(KO_set) + ',' + str(meansKO[str(KO_set)]) + '\n')

    run_result_file.close()


def resMod(r_a, r_b, r_start, r_flow, r_stable, r_gene_list, r_network_structure, nb_genes, KO_genes, OA_genes,
           model_type, boundary_model, network_values, nb_initial_states,
           nb_stable_states, data, output_file='output/Run_output'):
    """
    This function record the parameters and results from the boolean simulation into a text file.
    :param r_a:
    :param r_b:
    :param r_flow:
    :param r_stable:
    :param r_gene_list:
    :param r_network_structure:
    :param nom_genes:
    :param network_values:
    :param initial_states:
    :param stable_states:
    :param data:
    :return:
    """

    sep = "\n\n"

    # time stamp for the current run
    run_ref = str(datetime.datetime.now())
    ref_for_filename = run_ref.replace(':', '.')
    ref_for_filename = ref_for_filename.replace(' ', '.')

    # opening a csv file and filling it with the run results in easy to read format
    csv_file = output_file + '_' + ref_for_filename + '.csv'
    WriteOutputToCsv(csv_file, r_a, r_stable)

    # writing the flow as a csv file if starting from a single state
    if nb_initial_states == 1:
        csv_file_flow = output_file + '_' + ref_for_filename + '_stateflow_' + '.csv'
        WriteFlowToCsv(csv_file_flow, r_a, r_flow, r_stable)

    # opening text file for recording parameters and results summary
    filename = output_file + '_' + ref_for_filename + '.txt'
    f = open(filename, 'a')

    # filling content
    f.write('Date and time of run : ')
    f.write(run_ref)
    f.write(sep)
    f.write('\nStarting parameters :\n')
    f.write(r_gene_list)
    f.write(', \n')
    f.write('\nModel type : ')
    f.write(model_type)
    f.write(', \n')
    f.write('\nBoundary model : ')
    f.write(boundary_model)
    f.write(', \n')
    f.write('KO: ')
    f.write(str(KO_genes))
    f.write(sep)
    f.write('OA: ')
    f.write(str(OA_genes))
    f.write(sep)
    f.write(r_network_structure)
    f.write(', nb_columns_genes=1, name_index=0, nb_columns_network=3, network_headers=1')
    f.write(sep)
    f.write('List of genes in the network : ')
    f.write(str(r_a))
    f.write(sep)
    f.write('Interactions: ')
    f.write(str(r_b))
    f.write(sep)

    f.write("\n####### Results  #######\n#\n# Run parameters : ")
    f.write(str(nb_genes))
    f.write(" genes in the network ; ")
    f.write(str(network_values))
    f.write(" interactions in the network")
    f.write("\n# Number of initial states : ")
    f.write(str(nb_initial_states))
    f.write("\n# Number of stable states :")
    f.write(str(nb_stable_states))
    f.write("\n# Size of stable states and respective attraction basins :")
    f.write("\n# ")
    f.write(str(data))
    f.write("\n#\n########################")
    f.write(sep)
    f.write('Flow: ')
    f.write(str(r_flow))
    f.write(sep)
    f.write('Stable states: ')
    f.write(str(r_stable))
    f.write(sep)
    f.write('_________________________________________________________________________________________________________'
            '_____________________________________________________________')
    f.write(sep)
    f.close()
    return


def listing(stable_states, data, stable):
    g = open('Listing_new.txt', 'a')
    g.write('\n')
    g.write(str(stable))
    g.write('\n')
    return


def performance(time, initial_state_number):
    p = open('graphe.txt', 'a')
    p.write('\n')
    p.write(str(initial_state_number))
    p.write('\t')
    p.write(str(time))
    return
