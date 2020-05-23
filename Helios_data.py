import csv

def saveData(G):
    f = open("nodegraph.csv", 'w')
    with f :
        writer = csv.writer(f)
        for node in G.nodes(data=True) :
            writer.writerow(node)
