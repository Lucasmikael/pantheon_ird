from matplotlib.ticker import NullFormatter  # useful for `logit` scale
import matplotlib.pyplot as plt
import numpy as np
from matplotlib.backends.backend_tkagg import FigureCanvasTkAgg
import PySimpleGUI as sg
import matplotlib
matplotlib.use('TkAgg')
from pythonis_model import *
import networkx as nx
from grave import plot_network
from grave.style import use_attributes

"""
Demonstrates one way of embedding Matplotlib figures into a PySimpleGUI window.
Basic steps are:
 * Create a Canvas Element
 * Layout form
 * Display form (NON BLOCKING)
 * Draw plots onto convas
 * Display form (BLOCKING)

 Based on information from: https://matplotlib.org/3.1.0/gallery/user_interfaces/embedding_in_tk_sgskip.html
 (Thank you Em-Bo & dirck)
"""


# ------------------------------- PASTE YOUR MATPLOTLIB CODE HERE -------------------------------

# Fixing random state for reproducibility
# np.random.seed(19680801)
#
# # make up some data in the interval ]0, 1[
# y = np.random.normal(loc=0.5, scale=0.4, size=1000)
# y = y[(y > 0) & (y < 1)]
# y.sort()
# x = np.arange(len(y))
#
# # plot with various axes scales
# plt.figure(1)
#
# # linear
# plt.subplot(221)
# plt.plot(x, y)
# plt.yscale('linear')
# plt.title('linear')
# plt.grid(True)
#
#
# # log
# plt.subplot(222)
# plt.plot(x, y)
# plt.yscale('log')
# plt.title('log')
# plt.grid(True)
#
#
# # symmetric log
# plt.subplot(223)
# plt.plot(x, y - y.mean())
# plt.yscale('symlog', linthreshy=0.01)
# plt.title('symlog')
# plt.grid(True)
#
# # logit
# plt.subplot(224)
# plt.plot(x, y)
# plt.yscale('logit')
# plt.title('logit')
# plt.grid(True)

    # Path des fichiers en Import
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


genes_names_list, network_dictionary, unused_genes, network_as_list = ImportBooleanModel(genes_list_filename,network_filename)




def hilighter(event):
    # if we did not hit a node, bail
    if not hasattr(event, 'nodes') or not event.nodes:
        return

    # pull out the graph,
    graph = event.artist.graph

    # clear any non-default color on nodes
    for node, attributes in graph.nodes.data():
        attributes.pop('color', None)

    for u, v, attributes in graph.edges.data():
        attributes.pop('width', None)

    for node in event.nodes:
        graph.nodes[node]['color'] = 'C1'

        for edge_attribute in graph[node].values():
            edge_attribute['width'] = 3

    # update the screen
    event.artist.stale = True
    event.artist.figure.canvas.draw_idle()


def drawGraph(genes_names_list,network_as_list):
    G = nx.MultiDiGraph()


    G.add_nodes_from(genes_names_list)
    val_map = {'ARR23': 1.0,'CRF2': 0.5714285714285714}

    values = [val_map.get(node, 0.85) for node in G.nodes()]

#ajout des aretes et des noeuds
    for i in range (len(network_as_list)):
        gene_source = network_as_list[i][0]
        gene_target = network_as_list[i][2]
        interaction = network_as_list[i][1]
        G.add_edge(gene_source, gene_target, arrowstyle = interaction)

    # if gene_source not in G.nodes :
    #     if gene_target == gene_source :
    #         G.add_node(gene_source, regulation = "auto")
    #     else:
    #         G.add_node(gene_source, regulation = "none")

    print(G.adj)

    pos = nx.circular_layout(G)
    nx.draw_networkx_nodes(G, pos, node_size=1500, node_color= values)
    nx.draw_networkx_labels(G, pos)
    for edge in G.edges(data=True):
        if edge[2]['arrowstyle'] == "-1":
            arrowstyle = "-["
        else :
            arrowstyle = "-|>"
        nx.draw_networkx_edges(G, pos, edgelist=[(edge[0],edge[1])], arrowstyle = arrowstyle, node_size=1900)
# plt.show()

    plt.gca().yaxis.set_minor_formatter(NullFormatter())
    plt.subplots_adjust(top=1, bottom=0, left=0, right=1, hspace=0.25,
                    wspace=0.35)
    fig = plt.gcf()
    art = plot_network(G, node_style=use_attributes(),
                       edge_style=use_attributes())

    art.set_picker(10)     # if using Pyplot then get the figure from the plot

    return fig
# figure= fig.bbox.bounds

# ------------------------------- END OF YOUR MATPLOTLIB CODE -------------------------------

# ------------------------------- Beginning of Matplotlib helper code -----------------------
######Ca c'est bon#########

def draw_figure(canvas, figure, loc=(0, 0)):
    figure_canvas_agg = FigureCanvasTkAgg(figure, canvas)
    figure_canvas_agg.draw()
    figure_canvas_agg.get_tk_widget().pack(side='top', fill='both', expand=1)
    return figure_canvas_agg
# ------------------------------- Beginning of GUI CODE -------------------------------
def on_pick(event):
    artist = event.artist
    xmouse, ymouse = event.mouseevent.xdata, event.mouseevent.ydata
    x, y = artist.get_xdata(), artist.get_ydata()
    ind = event.ind
    print ('Artist picked:', event.artist)

def onclick(event):
    print('%s click: button=%d, x=%d, y=%d, xdata=%f, ydata=%f' %
          ('double' if event.dblclick else 'single', event.button,
           event.x, event.y, event.xdata, event.ydata))



layout = [[sg.Canvas(size=(640, 480), key='-CANVAS-')]]


# define the window layout
# layout = [[sg.Text('Plot test', font='Any 18')],
#           [sg.Canvas(size=(figure), key='canvas')]]
#
# # create the form and show it without the plot
window = sg.Window('Interaction Graph',
                   layout, finalize=True)

fig = drawGraph(genes_names_list,network_as_list)

# add the plot to the window
fig_canvas_agg = draw_figure(window['-CANVAS-'].TKCanvas, fig)
fig.canvas.mpl_connect('pick_event', on_pick)
fig.canvas.mpl_connect('button_press_event', onclick)
fig.canvas.mpl_connect('pick_event', hilighter)

canvas_elem = window['-CANVAS-']
canvas = canvas_elem.TKCanvas
event, values = window.read()

window.close()
