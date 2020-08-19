import matplotlib.pyplot as plt
from Helios_model import *

def drawStateActivationGraph(genes_selected_visu, flow, genes_names_list):

    for i in range (len(genes_selected_visu)):
        for gene in range (len(genes_names_list)):
            if genes_names_list[gene] == genes_selected_visu[i]:
                indexflow = gene
                value_source = getListfromFlow(flow,indexflow)
                list_y = createListPanelGraph(flow)


                fig_plot = plt.figure()
                plt.plot(list_y,value_source)
                fig_plot.suptitle("Activation graph of " +genes_selected_visu[i]+ " gene" , fontsize=20)
                plt.xlabel('State', fontsize=18)
                plt.ylabel('Activation', fontsize=16)
                plt.yticks([0, 1])
                plt.xticks(np.linspace(0,len(list_y)-1,len(list_y),endpoint=True))
                fig_plot.savefig(genes_selected_visu[i]+ " activity.png")
                plt.close()



def getListfromFlow(flow,indexflow):

    value_source = []
    for key, value in flow.items():
        # Cle du dictionnaire
        source = key
        for i in range (len(source)):
            if i == indexflow :
                value_source.append(source[i])
    print(value_source)
    return value_source
