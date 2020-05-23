import csv


nms = [['caca', 'pipi'], ['cici','pupu']]
f = open("nodegraph.csv", 'w')
with f :
    writer = csv.writer(f)
    for row in nms :
        writer.writerow(row)
