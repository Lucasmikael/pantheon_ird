# -*- coding: utf-8 -*-

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
import PySimpleGUI as sg
from numpy import array
from pythonis_tools import *
from pythonis_filesIO import *
from pythonis_model import *
from random_network import *
from Helios_model import *
from Helios_Refactor_Pythonis import *
from Helios_IO import *


# def fonctionMain(working_directory, genes_list_file, network_structure_file, nb_columns_genes=1, name_index=1,
#                  nb_columns_network=3, network_headers=0, nb_starting_states='10', KO_genes_param=['foo'],
#                  OA_genes_param=['foo'], model_type='logical', boundary_model='transient'):
#     os.chdir(working_directory)
#     a, b, c = ImportBooleanModel(working_directory, genes_list_file, network_structure_file, nb_columns_genes,
#                                  name_index,
#                                  nb_columns_network, network_headers)
#
#     flow, stable, start, stable_states, initial_states, genes_names, data, network, time = RunBooleanModel(a, b, 't',
#                                                                                                            't',
#                                                                                                            initial_state_number=nb_starting_states,
#                                                                                                            initial_state_choice='random',
#                                                                                                            stimulus=boundary_model,
#                                                                                                            initial_state_genes=[
#                                                                                                                'foo'],
#                                                                                                            model=model_type,
#                                                                                                            KO_genes=KO_genes_param,
#                                                                                                            OA_genes=OA_genes_param)
#
#     resMod(a, b, start, flow, stable, genes_list_file, network_structure_file, genes_names, KO_genes_param,
#            OA_genes_param, model_type, boundary_model, network, initial_states,
#            stable_states, data)
#     # listing(stable_states, data, stable)
#     # performance(time, initial_state_number=50000)
#
#     return a, b, flow, stable, start
#
#
# def generate_pairs(source):
#     pairs = []
#     for p1 in range(len(source)):
#         for p2 in range(p1 + 1, len(source)):
#             pairs.append([source[p1], source[p2]])
#     return pairs
#

#######
#
# Main program
#
#######


folder = "F:\paperPYTHONIS\LRPnetwork-201903"
gene_list_filename = "l_gnp-221018-161genes-Ju+PI-prior.txt"
network_filename = "TDCor6.32_output_221018_parallel.txt"
KO_genes_input = ['PLT1', 'ARF6', 'LRP1', 'PHB', 'TMO5', 'SHR', 'SCR', 'SHP1', 'ATML1', 'PID2']
OA_genes_input = ['PLT7', 'PUCHI', 'CRF1', 'ARF2', 'PLT5', 'ARF9', 'ARF17', 'U.box', 'ARF19', 'WRKY43']
more_KO_genes_input = [['PLT1'], ['PLT2'], ['PLT3'], ['PLT5'], ['PLT7'], ['ARF6'], ['ARF8'], ['PLT1', 'PLT2'],
                       ['PLT5', 'PLT7'], ['ARF6', 'ARF8'], ['PLT1', 'PLT2', 'PLT3'], ['PLT3', 'PLT5', 'PLT7']]
KO_mutation_type_list = ['none', 'single KO', 'multiple KO', 'combination KO']
OA_mutation_type_list = ['none', 'single OA', 'multiple OA', 'combination OA']
model_type_list = ['logical', 'algebraic']
boundary_model_list = ['transient', 'constant']
initial_states_choice_list = ['random', 'all_zeros', 'all_ones', 'specified', 'random-specified']

if __name__ == "__main__":

    #############
    # GUI setup
    #############

    # set background globally to dark grey and text to light grey
    sg.SetOptions(text_color="#e4e4e4", font='trajanpro 11', background_color="#343434")

    # layout for the title window - as it does nothing by itself, it has to be set to autoclose after 4sec when invoked
    layout_title = [[sg.Image(r'pantheon_front_50percent.gif')], ]
    window_title = sg.Window("title", layout=layout_title, auto_close=True, auto_close_duration=4, no_titlebar=True)
    event_title, values_title = window_title.Read()
    window_title.Close()

    # layout for the first window - import of gene list and gene network files
    layout_start = [[sg.Frame(layout=[[sg.Text('Gene List File', background_color='#343434', text_color='#e4e4e4')],
                                      [sg.Input(), sg.FileBrowse()],
                                      [sg.Text('Network Structure File', background_color='#343434',
                                               text_color='#e4e4e4')],
                                      [sg.Input(), sg.FileBrowse()]],
                              title='Import network', title_color='#e4e4e4', relief=sg.RELIEF_SUNKEN,
                              background_color='#343434',
                              tooltip='Choose gene list file and network structure file to import a predefined network model')],
                    [sg.Frame(layout=[[sg.Text('Number of nodes', background_color='#343434', text_color='#e4e4e4')],
                                      [sg.Input()],
                                      [sg.Text('Number of interactions', background_color='#343434',
                                               text_color='#e4e4e4')],
                                      [sg.Input()]],
                              title='Gene network Automated Initiation Algorithm (GAIA)', title_color='#e4e4e4',
                              relief=sg.RELIEF_SUNKEN,
                              background_color='#343434',
                              tooltip='Create a random network to work with')],
                    [sg.Frame(layout=[[sg.Radio('Batch Computation Mode', 'RADIO1', default=True, size=(18, 1),
                                                background_color='#343434', text_color='#e4e4e4'),
                                       sg.Radio('Visual Network Mode', 'RADIO1', background_color='#343434',
                                                text_color='#e4e4e4')]],
                              title='Run mode', title_color='#e4e4e4', relief=sg.RELIEF_SUNKEN,
                              background_color='#343434',
                              tooltip='Batch computation for whole network simulation and genes significance study / Interactive network to quick test hypothesis manually and visually')],
                    [sg.CloseButton('Import'), sg.CloseButton('Cancel')]]
    window_start = sg.Window("Initialize Network", alpha_channel=0.95, layout=layout_start)
    event_start, values_start = window_start.Read()
    window_start.Close()

    genes_list_filename = values_start[0]
    network_filename = values_start[1]
    GAIA_genes = values_start[2]
    GAIA_interactions = values_start[3]
    launchBatch = values_start[4]
    launchVisu = values_start[5]



    #########mettre condition If pour l'ouverture du batch mode ou de la visualisation
    if launchBatch :
        if (GAIA_genes != '') and (GAIA_interactions != ''):
            RandomGRN1("F:/paperPYTHONIS/LRPNetwork-201903/files_to_run".decode("utf-8"), genes_list_file='GAIA_genes.txt',
                       network_structure_file='GAIA_network.txt', number_of_nodes=GAIA_genes,
                       number_of_edges=GAIA_interactions)
            genes_names_list, network_dictionary, unused_genes, network_as_list, genes_non_sort = ImportBooleanModel(
                "F:/paperPYTHONIS/LRPNetwork-201903/files_to_run/GAIA_genes.txt",
                "F:/paperPYTHONIS/LRPNetwork-201903/files_to_run/GAIA_network.txt")
        else:
            genes_names_list, network_dictionary, unused_genes, network_as_list, genes_non_sort = ImportBooleanModel(genes_list_filename,
                                                                                                     network_filename)

        # Columns layout for the batch mode window
        col_batch1 = [[sg.Text('Input : List of genes', text_color='#e4e4e4', background_color='#343434')],
                      [sg.Listbox(values=genes_names_list,
                                  select_mode=sg.LISTBOX_SELECT_MODE_MULTIPLE, size=(20, 35),
                                  tooltip='Highlighted genes will be set to 1 at the start of each simulation run if the choice of genes initial states parameter is << specified >>')],
                      ]

        col_batch2 = [[sg.Text('Input : Interactions in the network', text_color='#e4e4e4', background_color='#343434')],
                      [sg.Listbox(values=network_as_list,
                                  select_mode=sg.LISTBOX_SELECT_MODE_MULTIPLE, size=(30, 35),
                                  tooltip='Format is [Source gene, interaction (-1: repression, 1 : activation), target gene]')],
                      ]

        col_subbatch3 = [[sg.Listbox(values=model_type_list, select_mode=sg.LISTBOX_SELECT_MODE_MULTIPLE, size=(15, 4),
                                     tooltip='Highlight the models to run through ARGOS'),
                          sg.Listbox(values=boundary_model_list, select_mode=sg.LISTBOX_SELECT_MODE_MULTIPLE, size=(15, 4),
                                     tooltip='Highlight the models to run through ARGOS')]]

        col_batch3 = [[sg.Frame(layout=[[sg.Text('Model type', text_color='#e4e4e4', background_color='#343434',
                                                 tooltip='Choice of rules to apply to the network structure'),
                                         sg.Combo(model_type_list)],
                                        [sg.Text('Boundary conditions', text_color='#e4e4e4', background_color='#343434',
                                                 tooltip='Choice of behavior for nodes without regulators'),
                                         sg.Combo(boundary_model_list)],
                                        [sg.Text('Genes initial states', text_color='#e4e4e4', background_color='#343434',
                                                 tooltip='In <<specified>> mode, use the gene list panel to highlight the genes to be set to 1 at start,\n all other genes will be set to 0 (0 or 1 in <<random-specified>> mode)'),
                                         sg.Combo(initial_states_choice_list)],
                                        [sg.Text('Number of initial states :', text_color='#e4e4e4',
                                                 background_color='#343434',
                                                 tooltip='Choice of how many initial states the model will run from'),
                                         sg.Radio('All possible states', 'RADIO2', background_color='#343434',
                                                  text_color='#e4e4e4'),
                                         sg.Radio('Given number of states', 'RADIO2', background_color='#343434',
                                                  text_color='#e4e4e4', default=True),
                                         sg.InputText(size=(10, 1), default_text='1')],
                                        [sg.Text('KO mutation type', text_color='#e4e4e4', background_color='#343434',
                                                 tooltip='Choice of Knocked Out mutation to apply to the network \n- single KO will run the model mutating each chosen genes one at a time'
                                                         '\n- multiple KO will run the model once with all chosen genes as KO'
                                                         '\n- combination KO will go through all combinations of double KO for the chosen genes and run the model for each'),
                                         sg.Combo(KO_mutation_type_list)],
                                        [sg.Text('List of KO genes', text_color='#e4e4e4', background_color='#343434',
                                                 tooltip='Enter gene names separated by a colon e.g. ABC, LMN, XYZ'),
                                         sg.InputText()],
                                        [sg.Text('OA mutation type', text_color='#e4e4e4', background_color='#343434',
                                                 tooltip='Choice of Over Activated mutation to apply to the network \n- single OA will run the model mutating each chosen genes one at a time'
                                                         '\n- multiple OA will run the model once with all chosen genes as KO'
                                                         '\n- combination OA will go through all combinations of double KO for the chosen genes and run the model for each'),
                                         sg.Combo(OA_mutation_type_list)],
                                        [sg.Text('List of OA genes', text_color='#e4e4e4', background_color='#343434',
                                                 tooltip='Enter gene names separated by a colon e.g. ABC, LMN, XYZ'),
                                         sg.InputText()],
                                        [sg.Button('Run PYTHONIS')]
                                        ],
                                title='Predict System Fates', title_color='#e4e4e4', background_color='#343434',
                                relief=sg.RELIEF_GROOVE)],
                      [sg.Text('', background_color='#343434')],
                      [sg.Frame(layout=[[sg.Text('Number of iterations', text_color='#e4e4e4', background_color='#343434',
                                                 tooltip='Size of the sample of random states that will be used as starting points for ARGOS'),
                                         sg.InputText(size=(8, 1), default_text='1')],
                                        [sg.Text('Choose the models to run through ARGOS :', text_color='#e4e4e4',
                                                 background_color='#343434')],
                                        [sg.Column(col_subbatch3, background_color='#343434')],
                                        [sg.Button('Run ARGOS')]],
                                title='Automatically Research Genes Of Significance', title_color='#e4e4e4',
                                background_color='#343434', relief=sg.RELIEF_GROOVE)],
                      [sg.Text('', background_color='#343434')],
                      [sg.Frame(layout=[[sg.Text('Placeholder - WIP algo', text_color='#e4e4e4', background_color='#343434',
                                                 tooltip='This will house the APOLLO algo call')],
                                        [sg.Text('Current Fate :', text_color='#e4e4e4', background_color='#343434'),
                                         sg.InputText()],
                                        [sg.Text('Targeted Fate :', text_color='#e4e4e4', background_color='#343434'),
                                         sg.InputText()],
                                        [sg.Button('Run APOLLO')]],
                                title='Challenge Fates', title_color='#e4e4e4', background_color='#343434',
                                relief=sg.RELIEF_GROOVE)]
                      ]

        layout_batchmode = [
            [sg.Column(col_batch1, background_color='#343434'), sg.Column(col_batch2, background_color='#343434'),
             sg.Column(col_batch3, background_color='#343434')],
            [sg.Text('Extract core network :', text_color='#e4e4e4', background_color='#343434'),
             sg.Button('Run HEPHAISTOS')],
            [sg.Button('Back to files import'), sg.Exit()]]

        window_batch = sg.Window("Batch Mode", alpha_channel=0.95, layout=layout_batchmode)

        # output window layout
        layout_pythonis_output = [[sg.Output(size=(80, 20))],
                                  [sg.Button('Update prediction', button_color=(sg.YELLOWS[0], sg.BLUES[0]))]
                                  ]

        window_pythonis_output = sg.Window('Pythonis Prediction', alpha_channel=0.95, layout=layout_pythonis_output,
                                           default_element_size=(30, 2))

        # TODO update this call to read the working directory from GUI file import
        # os.chdir("LRPNetwork-201903/files_to_run/output")

        while True:
            event_batch, values_batch = window_batch.Read()

            # catch all parameters from the GUI
            genes_selected = values_batch[0]
            interactions_selected = values_batch[1]
            boolean_model_selected = values_batch[2]
            boundary_model_selected = values_batch[3]
            initial_state_selected = values_batch[4]
            all_initial_states_bool = values_batch[5]
            subset_initial_states_bool = values_batch[6]
            number_initial_states = values_batch[7]
            KO_type_selected = values_batch[8]
            if KO_type_selected != 'none':
                KO_genes_selected = values_batch[9]
                KO_genes_selected = KO_genes_selected.rsplit(',')
            else:
                KO_genes_selected = ['foo']
            OA_type_selected = values_batch[10]
            if OA_type_selected != 'none':
                OA_genes_selected = values_batch[11]
                OA_genes_selected = OA_genes_selected.rsplit(',')
            else:
                OA_genes_selected = ['foo']
            number_iteration_ARGOS = values_batch[12]
            boolean_model_ARGOS = values_batch[13]
            boundary_model_ARGOS = values_batch[14]
            current_fate = values_batch[15]
            target_fate = values_batch[16]

            # launch modules or exit depending on button pushed
            if event_batch is None or event_batch == 'Exit':
                break
            elif event_batch == 'Run PYTHONIS':
                if all_initial_states_bool:
                    flow, stable, start, stable_states, initial_states, genes_names, data, network, time = RunBooleanModel(
                        genes_names=genes_names_list, genes_network=network_dictionary, initial_state_number='all',
                        initial_state_choice=initial_state_selected,
                        model=boolean_model_selected, stimulus=boundary_model_selected, initial_state_genes=genes_selected,
                        KO_genes=KO_genes_selected, OA_genes=OA_genes_selected)

                    resMod(genes_names_list, network_dictionary, start, flow, stable, genes_list_filename, network_filename,
                           genes_names, KO_genes_selected, OA_genes_selected, boolean_model_selected,
                           boundary_model_selected, network, initial_states, stable_states, data)

                elif subset_initial_states_bool:
                    event, value = window_pythonis_output.Read()
                    flow, stable, start, stable_states, initial_states, genes_names, data, network, time = RunBooleanModel(
                        genes_names=genes_names_list, genes_network=network_dictionary,
                        initial_state_number=number_initial_states, initial_state_choice=initial_state_selected,
                        model=boolean_model_selected, stimulus=boundary_model_selected, initial_state_genes=genes_selected,
                        KO_genes=KO_genes_selected, OA_genes=OA_genes_selected)

                    resMod(genes_names_list, network_dictionary, start, flow, stable, genes_list_filename, network_filename,
                           genes_names, KO_genes_selected, OA_genes_selected, boolean_model_selected,
                           boundary_model_selected, network, initial_states, stable_states, data)

                else:
                    sg.Popup(
                        'Please select value for number of initial states ( <<all>> or enter numerical value ) and run PYTHONIS again',
                        no_titlebar=True, background_color='#343434')

        window_batch.Close()





    ###################Visualization#######################

    if launchVisu :
        if (GAIA_genes != '') and (GAIA_interactions != ''):
            RandomGRN1("F:/paperPYTHONIS/LRPNetwork-201903/files_to_run".decode("utf-8"), genes_list_file='GAIA_genes.txt',
                       network_structure_file='GAIA_network.txt', number_of_nodes=GAIA_genes,
                       number_of_edges=GAIA_interactions)
            genes_names_list, network_dictionary, unused_genes, network_as_list, genes_non_sort = ImportBooleanModel(
                "F:/paperPYTHONIS/LRPNetwork-201903/files_to_run/GAIA_genes.txt",
                "F:/paperPYTHONIS/LRPNetwork-201903/files_to_run/GAIA_network.txt")
        else:
            genes_names_list, network_dictionary, unused_genes, network_as_list, genes_non_sort = ImportBooleanModel(genes_list_filename,
                                                                                                     network_filename)

        genes_names_list, network_dictionary, network_as_list = addElement(network_dictionary, genes_names_list,network_as_list)
        col_visu1 = [[sg.Text('Input : List of genes', text_color='#e4e4e4', background_color='#343434')],
                      [sg.Listbox(values=genes_names_list,
                                  select_mode=sg.LISTBOX_SELECT_MODE_MULTIPLE, size=(20, 35),
                                  tooltip='Highlighted genes will be set to 1 at the start of each simulation run if the choice of genes initial states parameter is << specified >>')],
                      ]

        col_visu2 = [[sg.Text('Input : Interactions in the network', text_color='#e4e4e4', background_color='#343434')],
                      [sg.Listbox(values=network_as_list,
                                  select_mode=sg.LISTBOX_SELECT_MODE_MULTIPLE, size=(30, 35),
                                  tooltip='Format is [Source gene, interaction (-1: repression, 1 : activation), target gene]')],
                      ]

        # col_subbatch3 = [[sg.Listbox(values=model_type_list, select_mode=sg.LISTBOX_SELECT_MODE_MULTIPLE, size=(15, 4),
        #                              tooltip='Highlight the models to run through ARGOS'),
        #                   sg.Listbox(values=boundary_model_list, select_mode=sg.LISTBOX_SELECT_MODE_MULTIPLE, size=(15, 4),
        #                              tooltip='Highlight the models to run through ARGOS')]]

        col_visu3 = [[sg.Frame(layout=[[sg.Text('Model type', text_color='#e4e4e4', background_color='#343434',
                                                 tooltip='Choice of rules to apply to the network structure'),
                                         sg.Combo(model_type_list)],
                                        [sg.Text('Boundary conditions', text_color='#e4e4e4', background_color='#343434',
                                                 tooltip='Choice of behavior for nodes without regulators'),
                                         sg.Combo(boundary_model_list)],
                                        [sg.Text('Genes initial states', text_color='#e4e4e4', background_color='#343434',
                                                 tooltip='In <<specified>> mode, use the gene list panel to highlight the genes to be set to 1 at start,\n all other genes will be set to 0 (0 or 1 in <<random-specified>> mode)'),
                                         sg.Combo(initial_states_choice_list)],
                                        [sg.Text('Number of initial states :', text_color='#e4e4e4',
                                                 background_color='#343434',
                                                 tooltip='Choice of how many initial states the model will run from'),
                                         sg.Radio('All possible states', 'RADIO2', background_color='#343434',
                                                  text_color='#e4e4e4'),
                                         sg.Radio('Given number of states', 'RADIO2', background_color='#343434',
                                                  text_color='#e4e4e4', default=True),
                                         sg.InputText(size=(10, 1), default_text='1')],
                                        [sg.Text('KO mutation type', text_color='#e4e4e4', background_color='#343434',
                                                 tooltip='Choice of Knocked Out mutation to apply to the network \n- single KO will run the model mutating each chosen genes one at a time'
                                                         '\n- multiple KO will run the model once with all chosen genes as KO'
                                                         '\n- combination KO will go through all combinations of double KO for the chosen genes and run the model for each'),
                                         sg.Combo(KO_mutation_type_list)],
                                        [sg.Text('List of KO genes', text_color='#e4e4e4', background_color='#343434',
                                                 tooltip='Enter gene names separated by a colon e.g. ABC, LMN, XYZ'),
                                         sg.InputText()],
                                        [sg.Text('OA mutation type', text_color='#e4e4e4', background_color='#343434',
                                                 tooltip='Choice of Over Activated mutation to apply to the network \n- single OA will run the model mutating each chosen genes one at a time'
                                                         '\n- multiple OA will run the model once with all chosen genes as KO'
                                                         '\n- combination OA will go through all combinations of double KO for the chosen genes and run the model for each'),
                                         sg.Combo(OA_mutation_type_list)],
                                        [sg.Text('List of OA genes', text_color='#e4e4e4', background_color='#343434',
                                                 tooltip='Enter gene names separated by a colon e.g. ABC, LMN, XYZ'),
                                         sg.InputText()],
                                        [sg.Button('Run HELIOS')]
                                        ],
                                title='Predict System Fates', title_color='#e4e4e4', background_color='#343434',
                                relief=sg.RELIEF_GROOVE)],
                      [sg.Text('', background_color='#343434')]]



        layout_visu = [
            [sg.Column(col_visu1, background_color='#343434'), sg.Column(col_visu2, background_color='#343434'),
             sg.Column(col_visu3, background_color='#343434')],
            [sg.Text('Extract core network :', text_color='#e4e4e4', background_color='#343434')],
            [sg.Button('Back to files import'), sg.Exit()]]



        window_visu = sg.Window("Initialize Network", alpha_channel=0.95, layout=layout_visu)

        ###########Execution visu##############

        event_visu, values_visu = window_visu.Read()

        # catch all parameters from the GUI
        genes_selected = values_visu[0]
        interactions_selected = values_visu[1]
        boolean_model_selected = values_visu[2]
        boundary_model_selected = values_visu[3]
        initial_state_selected = values_visu[4]
        all_initial_states_bool_visu = values_visu[5]
        subset_initial_states_bool_visu = values_visu[6]
        number_initial_states = values_visu[7]
        KO_type_selected = values_visu[8]
        if KO_type_selected != 'none':
            KO_genes_selected = values_visu[9]
            KO_genes_selected = KO_genes_selected.rsplit(',')
        else:
            KO_genes_selected = ['foo']
        OA_type_selected = values_visu[10]
        if OA_type_selected != 'none':
            OA_genes_selected = values_visu[11]
            OA_genes_selected = OA_genes_selected.rsplit(',')
        else:
            OA_genes_selected = ['foo']



##########Run Helios Window###########

        # launch modules or exit depending on button pushed
        if event_visu is None or event_visu == 'Exit':
            print("break")

        if event_visu == 'Run HELIOS':
            if all_initial_states_bool_visu:
                flow, stable, start, stable_states, initial_states, genes_names, data, network, time = RunBooleanModelVisu(
                    genes_names=genes_names_list, genes_network=network_dictionary, initial_state_number='all',
                    initial_state_choice=initial_state_selected,
                    model=boolean_model_selected, stimulus=boundary_model_selected, initial_state_genes=genes_selected,
                    KO_genes=KO_genes_selected, OA_genes=OA_genes_selected)

                resMod(genes_names_list, network_dictionary, start, flow, stable, genes_list_filename, network_filename,
                       genes_names, KO_genes_selected, OA_genes_selected, boolean_model_selected,
                       boundary_model_selected, network, initial_states, stable_states, data)

            elif subset_initial_states_bool_visu:
                flow, stable, start, stable_states, initial_states, genes_names, data, network, time = RunBooleanModelVisu(
                    genes_names=genes_names_list, genes_network=network_dictionary,
                    initial_state_number=number_initial_states, initial_state_choice=initial_state_selected,
                    model=boolean_model_selected, stimulus=boundary_model_selected, initial_state_genes=genes_selected,
                    KO_genes=KO_genes_selected, OA_genes=OA_genes_selected)

                # resMod(genes_names_list, network_dictionary, start, flow, stable, genes_list_filename, network_filename,
                #        genes_names, KO_genes_selected, OA_genes_selected, boolean_model_selected,
                #        boundary_model_selected, network, initial_states, stable_states, data)

            else:
                sg.Popup(
                    'Please select value for number of initial states ( <<all>> or enter numerical value ) and run PYTHONIS again',
                    no_titlebar=True, background_color='#343434')



            list_panel = createListPanelGraph(flow)
            layout_graph_drawing = ['bipartite_layout',
           'circular_layout',
           'kamada_kawai_layout',
           'random_layout',
           'rescale_layout',
           'shell_layout',
           'spring_layout',
           'spectral_layout',
           'planar_layout',
           'fruchterman_reingold_layout',
           'spiral_layout']

            window_visu.Close()

            nodewindow_active = False
            interactionwindow_active = False
            importdata_active = False
            savedata_active = False

            # col_visu1 = [[sg.Text('Input : List of genes', text_color='#e4e4e4', background_color='#343434')],
            #              [sg.Listbox(values=genes_names_list,
            #                          select_mode=sg.LISTBOX_SELECT_MODE_MULTIPLE, size=(20, 35),
            #                          tooltip='Highlighted genes will be set to 1 at the start of each simulation run if the choice of genes initial states parameter is << specified >>')],
            #              ]
            #
            # col_visu2 = [
            #     [sg.Text('Input : Interactions in the network', text_color='#e4e4e4', background_color='#343434')],
            #     [sg.Listbox(values=network_as_list,
            #                 select_mode=sg.LISTBOX_SELECT_MODE_MULTIPLE, size=(30, 35),
            #                 tooltip='Format is [Source gene, interaction (-1: repression, 1 : activation), target gene]')],
            #     ]
            #
            # # col_subbatch3 = [[sg.Listbox(values=model_type_list, select_mode=sg.LISTBOX_SELECT_MODE_MULTIPLE, size=(15, 4),
            # #                              tooltip='Highlight the models to run through ARGOS'),
            # #                   sg.Listbox(values=boundary_model_list, select_mode=sg.LISTBOX_SELECT_MODE_MULTIPLE, size=(15, 4),
            # #                              tooltip='Highlight the models to run through ARGOS')]]

            col_graphs1 = [
                         [sg.Frame(layout=[
                             [sg.Button('Launch visu'),sg.Combo(list_panel),sg.Combo(layout_graph_drawing)]],
                             title='Add element to graph', title_color='#e4e4e4',
                             background_color='#343434', relief=sg.RELIEF_GROOVE)],
                         [sg.Frame(layout=[
                             [sg.Button('Add Node'), sg.Button('Add Interaction')]],
                             title='Add element to graph', title_color='#e4e4e4',
                             background_color='#343434', relief=sg.RELIEF_GROOVE)],
                         [sg.Frame(layout=[
                             [sg.Button('Save graph'), sg.Button('Load graph')]],
                             title='Save or load Data', title_color='#e4e4e4',
                             background_color='#343434', relief=sg.RELIEF_GROOVE)],
                         [sg.Text('', background_color='#343434')]]

            layout_graphs = [
                [sg.Column(col_graphs1, background_color='#343434')],
                [sg.Text('Extract core network :', text_color='#e4e4e4', background_color='#343434')],
                [sg.Button('Back to files import'), sg.Exit()]]

            window_graphs = sg.Window("Initialize Network", alpha_channel=0.95, layout=layout_graphs)

            while True:

                event_graphs, values_graph = window_graphs.Read()
                graph_selected = values_graph[0]
                layout_selected = values_graph[1]

                if event_graphs is None or event_graphs == 'Exit':
                    break
                if event_graphs == "Launch visu":


                    layout_graph = [[sg.Canvas(size=(640, 480), key='-CANVAS-'),sg.Exit()]]

                    # define the window layout
                    # layout = [[sg.Text('Plot test', font='Any 18')],
                    #           [sg.Canvas(size=(figure), key='canvas')]]
                    #
                    # # create the form and show it without the plot
                    window_graph = sg.Window('Interaction Graph',
                                             layout_graph, finalize=True)



                    fig, G = drawGraph(genes_non_sort, network_as_list, flow,graph_selected,layout_selected)
                    # add the plot to the window
                    fig_canvas_agg = draw_figure(window_graph['-CANVAS-'].TKCanvas, fig)
                    fig.canvas.callbacks.connect('pick_event', on_pick)

                    canvas_elem = window_graph['-CANVAS-']
                    canvas = canvas_elem.TKCanvas
                    event, values = window_graph.read()
                    window_graph.close()
                    #canvas.delete('all')

                if event_graphs == "Save graph":
                    saveData(network_as_list,flow, genes_names_list)

                if event_graphs == "Load graph":
                    interaction_table, node_table = openData()
                    list_panel_load = createListPanelGraph(node_table)
                    importdata_active = False
                    layout_load=[
                        [sg.Button('Launch visu'), sg.Combo(list_panel_load), sg.Combo(layout_graph_drawing)]]
                    window_load = sg.Window('Load graph',
                                             layout_load, finalize=True)

                    event_load, values_load = window_load.Read()
                    graph_selected = values_load[0]
                    layout_selected = values_load[1]

                    if event_load == "Launch visu":
                        layout_graph_load = [[sg.Canvas(size=(640, 480), key='-CANVAS-')]]

                        # define the window layout
                        # layout = [[sg.Text('Plot test', font='Any 18')],
                        #           [sg.Canvas(size=(figure), key='canvas')]]
                        #
                        # # create the form and show it without the plot
                        window_graph_load= sg.Window('Interaction Graph',
                                                 layout_graph_load, finalize=True)

                        # fig, G = drawGraph(genes_non_sort, network_as_list, flow, graph_selected, layout_selected)
                        fig,G = drawLoadGraph(layout_selected,interaction_table, node_table, graph_selected)
                        # add the plot to the window

                        fig_canvas_agg = draw_figure(window_graph_load['-CANVAS-'].TKCanvas, fig)
                        fig.canvas.callbacks.connect('pick_event', on_pick)

                        canvas_elem = window_graph_load['-CANVAS-']
                        canvas = canvas_elem.TKCanvas
                        event, values = window_graph_load.read()
                        window_graph_load.close()


                if not nodewindow_active and event_graphs == 'Add Node':
                    nodewindow_active = True
                    layout_graph = [[sg.Text('The second window')],
                                    [sg.Input(key='-IN-')],
                                    [sg.Button('Show'), sg.Button('Exit')]]

                    # define the window layout
                    # layout = [[sg.Text('Plot test', font='Any 18')],
                    #           [sg.Canvas(size=(figure), key='canvas')]]
                    #
                    # # create the form and show it without the plot
                    nodewindow = sg.Window('Node/Interaction',
                                             layout_graph, finalize=True)


                if nodewindow_active :
                    event, values = nodewindow.read()

                    if event == 'Exit':
                        # print("Closing window 2", event)
                        nodewindow_active = False
                        nodewindow.close()
                    if event == 'Show':
                        sg.popup('You entered ', values['-IN-'])



                if not interactionwindow_active and event_graphs == 'Add Interaction' :

                    interactionwindow_active = True
                    layout_graph = [[sg.Text('The second window')],
                                    [sg.Input(key='-IN-')],
                                    [sg.Button('Show'), sg.Button('Exit')]]


                    interactionwindow = sg.Window('Add Interaction',
                                             layout_graph, finalize=True)

                if interactionwindow_active :

                    event, values = interactionwindow.read()

                    if event == "Exit" :
                        interactionwindow_active = False
                        interactionwindow.close()
                        if event == 'Show':
                            sg.popup('You entered ', values['-IN-'])





            window_graphs.close()

            # buttons_col = []
            #
            # for i in range(len(flow)):
            #     buttons_col.append([sg.Button('col {}'.format(i))])
            #
            # buttons_row = []
            # for i in range(len(flow)):
            #     buttons_row.append(sg.Button('row {}'.format(i)))
            #
            # layout_init_graph = [
            #     [sg.Text('Your typed chars appear here:'), sg.Text('', key='_OUTPUT_')],
            #     *buttons_col,
            #     [*buttons_row],
            #     [sg.Input(do_not_clear=True, key='_IN_')],
            #     [sg.Button('Show'), sg.Button('Exit')]
            # ]
            #
            # window_graph = sg.Window('Window Title').Layout(layout_init_graph)
            #
            # while True:             # Event Loop
            #     event, values = window_graph.Read()
            #     print(event, values)
            #     if event is None or event == 'Exit':
            #         break
            #     if event == 'Show':
            #         # change the "output" element to be the value of "input" element
            #         window_graph.FindElement('_OUTPUT_').Update(values['_IN_'])
            #
            # window_graph.Close()




            layout_graph = [[sg.Canvas(size=(640, 480), key='-CANVAS-')]]


            # define the window layout
            # layout = [[sg.Text('Plot test', font='Any 18')],
            #           [sg.Canvas(size=(figure), key='canvas')]]
            #
            # # create the form and show it without the plot
            window_graph = sg.Window('Interaction Graph',
                               layout_graph, finalize=True)

            fig = drawGraph(genes_non_sort,network_as_list,flow)
            # add the plot to the window
            fig_canvas_agg = draw_figure(window_graph['-CANVAS-'].TKCanvas, fig)
            fig.canvas.callbacks.connect('pick_event', on_pick)

            canvas_elem = window_graph['-CANVAS-']
            canvas = canvas_elem.TKCanvas
            event, values = window_graph.read()

            window_graph.close()



"""
    if mutant_flag == 'singleKO' :
        for m_g in KO_genes_input:
            for model in model_type_list:
                for boundary in boundary_model_list:
                    a, b, flow, stable, start = fonctionMain(folder.decode("utf-8"), gene_list_filename,
                                                             network_filename, nb_starting_states='1000', KO_genes_param = m_g, model_type= model, boundary_model= boundary)


    elif mutant_flag == 'singleOA' :
        for m_g in OA_genes_input:
            for model in model_type_list:
                for boundary in boundary_model_list:
                    a, b, flow, stable, start = fonctionMain(folder.decode("utf-8"), gene_list_filename,
                                                             network_filename, nb_starting_states='1000', OA_genes_param = m_g, model_type= model, boundary_model= boundary)

    elif mutant_flag == 'doubleKO' :
        KO_pairs = generate_pairs(KO_genes_input)
        for pairs in KO_pairs:
            for model in model_type_list:
                for boundary in boundary_model_list:
                    a, b, flow, stable, start = fonctionMain(folder.decode("utf-8"), gene_list_filename,
                                                             network_filename, nb_starting_states='1000', KO_genes_param = pairs, model_type= model, boundary_model= boundary)

    elif mutant_flag == 'doubleOA':
        KO_pairs = generate_pairs(OA_genes_input)
        for pairs in KO_pairs:
            for model in model_type_list:
                for boundary in boundary_model_list:
                    a, b, flow, stable, start = fonctionMain(folder.decode("utf-8"), gene_list_filename,
                                                             network_filename, nb_starting_states='1000', OA_genes_param=pairs, model_type=model, boundary_model=boundary)


"""

#   a, b, flow, stable, start = fonctionMain(folder.decode("utf-8"), gene_list_filename,
#                                     network_filename, nb_columns_genes = 1, name_index = 1, nb_columns_network = 3, network_headers = 0, KO_genes_param = KO_genes_input, OA_genes_param = OA_genes_input)


"""
if __name__ == "__main__":
    a, b, c = ImportBooleanModel("C:\Pythonis\Files_to_run".decode("utf-8"),
                                 "genes_list.txt",
                                 "network_structure.txt")
    core = ExtractCoreNetwork(b)
    print 'full network keys & values'
    print b
    print len(b.keys())
    print '\n', len(b.values())
    print '\ncore network keys & values'
    print core
    print len(core.keys())
    print '\n', len(core.values())
    flow, stable, start = RunBooleanModel(a, core, 't', 't', initial_state_number=50,
                                          initial_state_choice='random',
                                          initial_state_genes=['foo'],
                                          model='algebraic',
                                          stimulus='transient',
                                          verbose=True)

"""

##########
#
# TODO : when given a set of stable states, extract the list of gene which have similar or different behavior between those stable states
#
# TODO : create a shuffle function for network interactions and/or interactions nature and write a bootstrap ?
#
# TODO : create a clustering function for stable states --> mini project 001 for Awa
#
# TODO : filter the LRP dataset to get potential initial states for the  given genes
#
# TODO : script to return the flow from a given state
#
# TODO : see if we can compatiment / classify the trajectories depending on the initial active genes
#
# TODO : make the script an executable
#
##########
