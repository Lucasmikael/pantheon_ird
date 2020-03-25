from matplotlib.ticker import NullFormatter  # useful for `logit` scale
import matplotlib.pyplot as plt
import numpy as np
from matplotlib.backends.backend_tkagg import FigureCanvasTkAgg
import PySimpleGUI as sg
import matplotlib
matplotlib.use('TkAgg')
from pythonis_model import *
import networkx as nx



class Node:
    """
    Classe de creation des graphs
    """
    def __init__(self,name,ativity):
        print("creation d'un humain")
        self.nom = name
        seld=f.activite = activity



genes_list_filename = "l_gnp-190515.txt"
network_filename = "jeu_essai_helios.txt"

# Values en dur pour le run du modèle booléen
boolean_model_selected = "logical"
boundary_model_selected = "transient"
initial_state_selected = "random"
number_initial_states = "1"
KO_genes_selected = ['foo']
OA_genes_selected = ['foo']
genes_selected = ['foo']

aretes = []


genes_names_list, network_dictionary, unused_genes, network_as_list = ImportBooleanModel(genes_list_filename,

###########node###########

for i in range (len(genes_names_list))
    nodes.append(Node(genes_names_list[i]))

print(nodes)
