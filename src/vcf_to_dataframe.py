#!/usr/bin/env python3

"""
@author: Yashwanth Lagisetty

Main script to generate dataframe for EAWavelet

"""

import pandas as pd 
import numpy as np 
from pysam import VariantFile
import sys
import re
from collections import OrderedDict

vcffile = sys.argv[1]
networkfile = sys.argv[2]
casesfile = sys.argv[3]
controlsfile = sys.argv[4]
genefile = sys.argv[5]
state = sys.argv[6]

savepath = sys.argv[7]

def init(vcffile,networkfile,casesfile,controlsfile,genefile):
	print('Initializing \n')
	networkmap = {}
	repgenelist = []
	with open(networkfile) as f:
		for line in f:
			node1,node2,string = line.strip('\n').split(',')
			if node2 != 'NOLINK':
				edge = node1+'_'+node2
				networkmap[edge] = string
				repgenelist.append(node1)
				repgenelist.append(node2)
	netgenelist = list(set(repgenelist))
	if isinstance(genefile, str):
		vcfgenelist = [line.strip('\n') for line in open(genefile)]
		genelist = list(set(netgenelist) & set(vcfgenelist))
	else:
		genelist = natgenelist
	cases = [line.strip('\n') for line in open(casesfile)]
	controls = [line.strip('\n') for line in open(controlsfile)]
	samples = cases+controls

	vcf=VariantFile(vcffile)
	EA_matrix = np.zeros((len(samples),len(genelist)))
	masterdict = {gene:{sample:[] for sample in samples} for gene in genelist}
	samp_map = OrderedDict((samples,i) for i,samples in enumerate(samples))
	gene_map = OrderedDict((gene,i) for i,gene in enumerate(genelist))

	return genelist, cases, controls, samples, vcf, EA_matrix,masterdict,samp_map,gene_map

def retrieve_EA(recEA):
	EAs = []
	for score in recEA:
		try:
			score = float(score)
		except(ValueError,TypeError):
			if(re.search(r'fs-indel',score) or re.search(r'STOP',score)):
				score=100
			else:
				score=0
		EAs.append(score)
		geneEA = np.max(EAs)
	return geneEA


def retrieve_pph2(rec):
	pph2s = []
	pph2sUnfiltered = rec
	pph2s = [x for x in pph2sUnfiltered if x is not None]
	genePPh2 = np.max(pph2s)
	return genePPh2


def search_vcf(vcf,genelist,samples,masterdict,state):
	print('Creating Matrix \n')
	for rec in vcf.fetch():
		gene = rec.info['gene']
		if isinstance(gene,tuple):
			gl = list(gene)
			gene = gl[0]
		if gene == 'UNKNOWN':
			continue
		if gene == '.':
			continue
		if gene in genelist:
			zyg = 0 
			if state == 'EA':
				geneEA = retrieve_EA(rec.info['EA'])
			elif state == 'PPh2':
				if 'dbNSFP_Polyphen2_HDIV_score' not in rec.info.keys():
					EA = rec.info['EA']
					if EA == 'silent':
						genePPh2 = 0
					elif EA == 'STOP':
						genePPh2 = 1
					elif EA == 'fs-indel':
						genePPh2 = 1
					else:
						continue
				genePPh2 = retrieve_pph2(rec.info['dbNSFP_Polyphen2_HDIV_score'])
			for sample in samples:
				genotype = rec.samples[sample]['GT']
				if genotype == (0,0):
					zyg = 0 
				elif genotype == (1,1):
					zyg = 2
				else:
					zyg = 1
				if state == 'EA':
					masterdict[gene][sample].append((1-(geneEA/100))**zyg)
				elif state == 'PPh2':
					masterdict[gene][sample].append((1-(genePPh2))**zyg) #genepph2 is not /100 because it is already normalize 0-1
	return masterdict

def make_matrix(masterdict,genelist,gene_map,samp_map,EA_matrix,samples):
	missing_genes = []
	for gene in genelist:
		for sample in samples:
			sampidx = samp_map[sample]
			geneidx = gene_map[gene]
			varlist = masterdict[gene][sample]
			if len(varlist) == 0:
				EA_matrix[sampidx,geneidx] = 0 
				missing_genes.append(gene)
			else:
				pEA = 1-np.prod(varlist)
				EA_matrix[sampidx,geneidx] = pEA
	missing_genes_unq = list(set(missing_genes))
	return EA_matrix, missing_genes_unq

def main():
	genelist, cases, controls, samples, vcf, EA_matrix,masterdict,samp_map,gene_map = init(vcffile,networkfile,casesfile,controlsfile,genefile)
	masterdict = search_vcf(vcf,genelist,samples,masterdict,state)
	EA_matrix,missing_genes_unq = make_matrix(masterdict,genelist,gene_map,samp_map,EA_matrix,samples)

	df = pd.DataFrame(data=EA_matrix,index=list(samp_map.keys()),columns=list(gene_map.keys()))
	df.to_csv(savepath+'design_matrix_'+state+'.txt',sep='\t',header=True,index=True)
	with open(savepath+'missing_genes.txt','w') as f:
		for gene in missing_genes_unq:
			f.write('%s\n' %gene)
	print('DataFrame created and saved\n')

if __name__ == "__main__":
	main()









