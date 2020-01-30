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
# genelist=[line.strip('\n').split(',')[0] for line in open(genelistfile)]
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
		EAs = []
		zyg = 0
		geneEA = 0
		for score in rec.info['EA']:
			try:
				score = float(score)
			except(ValueError,TypeError):
				if(re.search(r'fs-indel',score) or re.search(r'STOP',score)):
					score = 100
				else:
					score = 0
			EAs.append(score)
		geneEA = np.max(EAs)
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
			masterdict[gene][sample].append((1-(geneEA/100))**zyg)
missing_genes = []
for gene in genelist:
	for sample in samples:
		sampidx = samp_map[sample]
		geneidx = gene_map[gene]
		var_list = masterdict[gene][sample]
		if len(var_list) == 0:
			EA_matrix[sampidx,geneidx] = 1
			# treat missing genes as wild types with EA=0 and only 1 variant in gene
			# if genefile is included, then only genes in both network and genefile 
			# will be used, and this IF loop will not be called
			missing_genes.append(gene)
		else:
			boltprob = 1 - np.prod(var_list)
			EA_matrix[sampidx,geneidx] = boltprob
missing_genes_unq = list(set(missing_genes)) 
# print('Total number of genes in network is: {} number of genes not seen in cohort is: {}'. format(len(genelist),w))
if save == True:
	df = pd.DataFrame(data=EA_matrix,index=list(samp_map.keys()),columns=list(gene_map.keys()))
	df.to_csv(savepath+'bolt_met_design_matrix.txt',sep='\t',header=True,index=True)
	with open(savepath+'missing_genes', 'w') as f:
		for gene in missing_genes_unq:
			f.write('%s\n' %gene)
	print('DataFrame Saved')