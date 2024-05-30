import networkx as nx



def parseFile(fileName, species_id):
	network = nx.Graph()
	fRead = open(fileName , 'r')
	fRead.readline()
	allowed_capture = ['Two-hybrid', 'Affinity Capture-Luminescence', 'Affinity Capture-MS', 'Affinity Capture-RNA', 'Affinity Capture-Western', 'physical']

	for line in fRead:
		splitted = line.strip().split('\t')
		#print(splitted)
		if splitted[11] in allowed_capture and splitted[12] in allowed_capture and splitted[15] == species_id and splitted[16] == species_id:
			if splitted[1] != splitted[2]:
				network.add_edge(splitted[1], splitted[2])
#				if splitted[1] == '101928094' or splitted[2] == '101928094':
#					print(splitted)

	fRead.close()
	
	print(f"genes: {network.number_of_nodes()}, \n interactions: {network.number_of_edges()}\n, density: {nx.density(network)}\n")
	return network



"""
	Main code starts here
"""

biogridFile = 'input/BIOGRID-ALL-4.4.218.tab3.txt'
species = '9606'
outputFile = 'output/Human_PPI_General.edgelist'

print ("Loading Interactions from BioGrid file\n")

network = parseFile(biogridFile, species)

print ("## Output PPI network\n")
nx.write_edgelist(network, outputFile, data=False)
nodes = [n for n in network.nodes()]

