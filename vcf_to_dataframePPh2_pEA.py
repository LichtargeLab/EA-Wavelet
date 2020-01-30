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
savepath = sys.argv[5]
genefile = sys.argv[6]

save = True

networkmap = {}
repgenelist=[]
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
	genelist = netgenelist
cases = [line.strip('\n') for line in open(casesfile)]
controls = [line.strip('\n') for line in open(controlsfile)]
samples = cases+controls

#Initialize design matrix
print('Initializeing \n')
vcf = VariantFile(vcffile)
# samples = list((vcf.header.samples))
EA_matrix = np.zeros((len(samples),len(genelist)))
masterdict = {gene:{sample:[] for sample in samples} for gene in genelist}
gene_boltz = {gene:1 for gene in genelist}
samp_map = OrderedDict((samples,i) for i,samples in enumerate(samples))
gene_map = OrderedDict((gene,i) for i,gene in enumerate(genelist))

print('Creating matrix \n')
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
		genePPh2 = 0
		zyg = 0
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
		else:
			pph2sUnfiltered = rec.info['dbNSFP_Polyphen2_HDIV_score']
			pph2s = [x for x in pph2sUnfiltered if x is not None]
			genePPh2 = np.max(pph2s)
		for sample in samples:
			# sampidx = samp_map[sample]
			# geneidx = gene_map[gene]
			genotype = rec.samples[sample]['GT']
			if genotype == (0,0):
				zyg = 0
			elif genotype == (1,1):
				zyg = 2
			else:
				zyg = 1
			masterdict[gene][sample].append((1-(genePPh2/100))**zyg)
missing_genes = []
for gene in genelist:
	for sample in samples:
		sampidx = samp_map[sample]
		geneidx = gene_map[gene]
		var_list = masterdict[gene][sample]
		if len(var_list) == 0:
			EA_matrix[sampidx,geneidx] = 0
			missing_genes.append(gene)
		else:
			boltprob = 1 - np.prod(var_list)
			EA_matrix[sampidx,geneidx] = boltprob
missing_genes_unq = list(set(missing_genes)) 

if save == True:
	df = pd.DataFrame(data=EA_matrix,index=list(samp_map.keys()),columns=list(gene_map.keys()))
	df.to_csv(savepath+'PPh2_pEA_design_matrix.txt',sep='\t',header=True,index=True)
	with open(savepath+'missing_genes', 'w') as f:
		for gene in missing_genes_unq:
			f.write('%s\n' %gene)
	print('DataFrame Saved')











