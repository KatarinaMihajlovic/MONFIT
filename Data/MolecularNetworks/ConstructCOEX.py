import os
import networkx as nx


def create_gene_coexpress_network(directory_path,threshold_value=3.):
	# Return a networkx graph
	all_files = os.listdir(directory_path)
	G = nx.Graph()
	for e1 in all_files:
		with open(directory_path+e1, 'r') as i_stream:
			for line in i_stream.readlines():
				g1 = e1
				lspt = line.strip().split('\t')
				g2 = lspt[0]
				zscore = float(lspt[1])
				if zscore >= threshold_value:
					G.add_edge(g1, g2)
				else:
					break
	print("density:\t\t", nx.density(G))
	print("nodes:\t\t", G.number_of_nodes())
	print("edges:\t\t", G.number_of_edges())
	return G


coex = create_gene_coexpress_network('input/Hsa-u.v22-05.G16651-S245698.combat_pca.subagging.z.d/')
outputFile = 'output/Human_COEX_General.edgelist'
nx.write_edgelist(coex, outputFile)



