import networkx as nx

def parseFile(fileName, species_id):
	network = nx.Graph()
	fRead = open(fileName , 'r')
	fRead.readline()
	allowed_capture = ['genetic']

	for line in fRead:
		splitted = line.strip().split('\t')
		#print(splitted)
		if splitted[12] in allowed_capture and splitted[15] == species_id and splitted[16] == species_id:
			if splitted[1] != splitted[2]:
				network.add_edge(splitted[1], splitted[2])

	
	fRead.close()
	
	print(f"genes: {network.number_of_nodes()}, \n interactions: {network.number_of_edges()}\n, density: {nx.density(network)}\n")
	return network




"""
	Main code starts here
"""


biogridFile = 'input/BIOGRID-ALL-4.4.218.tab3.txt'
species = '9606'
outputFile = 'output/Human_GI_General.edgelist'

print ("Loading Interactions from BioGrid file\n")
network = parseFile(biogridFile, species)


print ("## Output GI network\n")
nx.write_edgelist(network, outputFile, data=False)





