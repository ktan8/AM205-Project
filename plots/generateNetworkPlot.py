import networkx as nx
import matplotlib.pyplot as plt
import re
from networkx.drawing.nx_agraph import graphviz_layout

#data = open("../BIOMD0000000313/jacobian.csv")
data = open("../data/PDL1_pathway/PDL1_pathway.matrix.csv")
dataList = list()


G = nx.DiGraph()

header = data.readline()
headerList = header.strip().split(",")
# Loop through file
for line in data:
	lineArr = line.strip().split(",")
	rowGene = lineArr[0]
	# Loop through the elements in the line
	# and add edges as needed.
	for element in range(1, len(lineArr)):

		#rowGeneColIndex = headerList.index(rowGene)
		colGene = headerList[element]
		statusVal = lineArr[element]
		#print statusVal
		if int(statusVal) == int(1) or int(statusVal) == int(-1):
			rowGene = re.sub('PD-L1', "PD_L1", rowGene)
			colGene = re.sub('PD-L1', "PD_L1", colGene)
			rowGene = re.sub('-', "\n", rowGene)
			colGene = re.sub('-', "\n", colGene)

			#G.add_edge(rowGene, colGene, weight=statusVal)
			G.add_edge(rowGene, colGene, weight=statusVal)
			


		# element = lineArr
		# print 

	# G.add_edges_from

# Need to create a layout when doing
# separate calls to draw nodes and edges
pos = nx.spring_layout(G, k=0.55, iterations=100)
#pos = nx.kamada_kawai_layout(G)
#pos = nx.random_layout(G)
plt.figure()
nx.draw_networkx_nodes(G, pos, cmap=plt.get_cmap('jet'), 
                        node_size = 500)
nx.draw_networkx_labels(G, pos)
nx.draw_networkx_edges(G, pos, edge_color='black', arrows=True)
plt.savefig("PD_L1.pathway.defaultNetworkX.png", format="PNG")
#plt.show()



color_map = []
for node in G.nodes():
	if node == "PD_L1":
		color_map.append("pink")
	elif node in ["JAK2", "STAT5"]:
		color_map.append("lightgreen")
	else:
		color_map.append("lightblue")
plt.figure()
nx.draw(G, pos=graphviz_layout(G), node_size=1200, node_color=color_map,
linewidths=0.25, font_size=10, font_weight='bold', with_labels=True)
fig = plt.gcf()
fig.set_size_inches(10.5, 7.5)
plt.savefig("PD_L1.pathway.graphViz.png", format="PNG")
#plt.show()


