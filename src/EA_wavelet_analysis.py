"""
@author: Yashwanth lagisetty

Main script to run EAWavelet

"""

import networkx as nx 
import numpy as np 
import pandas as pd 
import seaborn as sb 
import scipy.stats as stats 
from sklearn.decomposition import PCA
from sklearn.preprocessing import StandardScaler
import sys 

### Set directory for graphwave repository 
graphwavepath = sys.argv[1]

sys.path.insert(0,graphwavepath)
import graphwave
from graphwave.graphwave import *

dataframe = sys.argv[2]
cases = sys.argv[3]
conts = sys.argv[4]
networkfile = sys.argv[5]

def init(dataframe,cases,conts,network):
	print('Initializing...')
	data = pd.read_csv(dataframe,delimiter='\t',header=0,index_col=0)
	genelist = data.columns.values.tolist()
	caselist = [line.strip('\n') for line in open(cases)]
	contlist = [line.strip('\n') for line in open(conts)]
	shift = None
	case = nx.Graph()
	cont = nx.Graph()

	return data, genelist, caselist,contlist, shift, case,cont 

def create_graphs(data,case,cont,networkfile,caselist,contlist,genelist):
	print('Creating Case and Control Graphs... \n')
	genelistnet = []
	with open(networkfile) as f:
		for line in f:
			node1,node2,_ = line.strip('\n').split(',')
			if node1 not in genelist:
				continue
			if node2 not in genelist:
				continue
			if node2 != 'NOLINK':
				casemet = np.mean(list(np.abs(data.loc[caselist,node1]-data.loc[caselist,node2])))
				contmet = np.mean(list(np.abs(data.loc[contlist,node1]-data.loc[contlist,node2])))
				case.add_node(node1)
				case.add_node(node2)
				cont.add_node(node1)
				cont.add_node(node2)
				case.add_edge(node1,node2,weight=casemet)
				cont.add_edge(node1,node2,weight=contmet)
				genelistnet.append(node1)
				genelistnet.append(node2)
			else:
				continue
	genelistnet = list(set(genelistnet))
	return case,cont,genelistnet

#write out arguments of function 
def run_wavelet_decomp(case_graph,cont_graph, genelist, shift):
	print('Computing embeddings...\n')
	CaseChi,CaseHeat,CaseTau = graphwave_alg(case_graph,np.linspace(0,100,25),taus='auto',verbose=True)
	ContChi,ContHeat,ContTau = graphwave_alg(cont_graph,np.linspace(0,100,25),taus='auto',verbose=True)

	if CaseChi.shape[0] != ContChi.shape[0]:
		raise ValueError('Case and Control embeddings have different number of nodes')
	else:
		shift = CaseChi.shape[0]

	case_node_map = {}
	i1 = 0
	for node in list(case_graph.nodes):
		case_node_map[node] = i1
		i1 +=1
	cont_node_map = {}
	i2 = 0
	for node in list(cont_graph.nodes):
		cont_node_map[node] = i2
		i2 +=1

	TotChi = np.vstack((CaseChi,ContChi))

	#Check if stacking occured correctly 
	randgenes = np.random.choice(genelist,size=5,replace=False)

	for gene in randgenes:
		if np.array_equal(CaseChi[case_node_map[gene]],TotChi[case_node_map[gene]]) != True:
			raise ValueError('CaseChi and TotChi not mapped correctly!!')
		if np.array_equal(ContChi[cont_node_map[gene]],TotChi[cont_node_map[gene]+shift]) != True:
			raise ValueError('ContChi and TotChi not mapped correctly!!')
	print('Passed random mapping check')

	return TotChi,shift, case_node_map, cont_node_map 

#writeout arguments of function
def PCA_distance(embeddings,case_node_map,cont_node_map,shift):
	pca = PCA(n_components=3)
	transdata = pca.fit_transform(StandardScaler().fit_transform(embeddings))

	distances = pd.DataFrame(columns=['distance','gene'])
	gcounter = 0 
	for gene in case_node_map.keys():
		genepidx = case_node_map[gene]
		genenidx = cont_node_map[gene]+shift 
		dist = (transdata[genepidx,0]-transdata[genenidx,0])**2 + (transdata[genepidx,0] - transdata[genenidx,0])**2
		distances.loc[gcounter,'gene'] = gene
		distances.loc[gcounter,'distance'] = np.sqrt(dist)
		gcounter +=1

	distances.sort_values('distance',ascending=False).reset_index(drop=True).to_csv('./PCAdistances.txt',sep='\t',header=True,index=False)
	print('Distance matrix tabulated and saved')
	return distances 

def main():
	data,genelist,caselist,contlist,shift,case,cont = init(dataframe,cases,conts,networkfile)
	case,cont,genelistnet = create_graphs(data,case,cont,networkfile,caselist,contlist,genelist)
	TotChi,shift,case_node_map,cont_node_map = run_wavelet_decomp(case,cont,genelistnet,shift)
	distances = PCA_distance(TotChi,case_node_map,cont_node_map,shift)
	print('EA Wavelet Analysis Complete')

if __name__ == "__main__":
	main()
