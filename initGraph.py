from pythonis_filesIO import *
from pythonis_model import *
from random_network import *



if __name__=="__main__":

    #Path des fichiers en Import
    genes_list_filename = "l_gnp-190515.txt"
    network_filename = "TDCor6.32_output_Missinglink_parallel-190515-cytoscape.txt"


    #Values en dur pour le run du modèle booléen
    boolean_model_selected = "logical"
    boundary_model_selected = "transient"
    initial_state_selected = "random"
    number_initial_states = "1"
    KO_genes_selected = ['foo']
    OA_genes_selected = ['foo']
    genes_selected = ['foo']







    genes_names_list, network_dictionary, unused_genes, network_as_list = ImportBooleanModel(genes_list_filename,
                                                                                         network_filename)

    flow, stable, start, stable_states, initial_states, genes_names, data, network, time = RunBooleanModel(
        genes_names=genes_names_list, genes_network=network_dictionary,
        initial_state_number=number_initial_states, initial_state_choice=initial_state_selected,
        model=boolean_model_selected, stimulus=boundary_model_selected, initial_state_genes=genes_selected,
        KO_genes=KO_genes_selected, OA_genes=OA_genes_selected)