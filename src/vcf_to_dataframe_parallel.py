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
from joblib import Parallel, delayed
from tqdm import tqdm

def parallelize(vcf_file, reference, samples, n_jobs, state='EA'):
	gene_dfs = Parallel(n_jobs=n_jobs)(
		delayed(parse_gene)(vcf_file=vcf_file, samples_list=samples, gene=gene, gene_ref=reference.loc[gene],state=state) for gene in tqdm(reference.index.unique()))
	design_matrix = pd.concat(gene_dfs, axis=1)
	design_matrix = 1 - design_matrix
	return design_matrix


def retrieve_EA(rec):
	EAs = []
	for score in rec:
		try:
			score = float(score)
		except(ValueError, TypeError):
			if (re.search(r'fs-indel', score) or re.search(r'STOP', score)):
				score = 100
			else:
				score = 0
		EAs.append(score)
		geneEA = np.max(EAs)
	return geneEA


def retrieve_pph2(rec):
	pph2s = []
	pph2sUnfiltered = rec
	pph2s = [x for x in pph2sUnfiltered if x is not None]
	genePPh2 = np.max(pph2s)
	return genePPh2


def get_zygo(genotype):
    if genotype == (0, 0):
        zyg = 0
    elif genotype == (1, 1):
        zyg = 2
    else:
        zyg = 1
    return zyg


def parse_gene(vcf_file, samples_list, gene, gene_ref, state):
	vcf = VariantFile(vcf_file)
	vcf.subset_samples(samples_list)
	if type(gene_ref) == pd.DataFrame:
		chrom = gene_ref['chrom'].iloc[0].strip('chr')
	else:
		chrom = gene_ref['chrom'].strip('chr')
	start = gene_ref['cdsStart'].min()
	end = gene_ref['cdsEnd'].max()
	temp_df = pd.DataFrame(np.ones((len(samples_list), 1)), index=samples_list, columns=[gene])

	for rec in vcf.fetch(contig=chrom, start=start, end=end):
		rec_gene = rec.info['gene']
		if isinstance(rec_gene, tuple):
			gl = list(rec_gene)
			rec_gene = gl[0]
		if rec_gene != gene:
			continue
		if state == 'EA':
			score = retrieve_EA(rec.info['EA'])
			score = score / 100
		elif state == 'PPh2':
			if 'dbNSFP_Polyphen2_HDIV_score' not in rec.info.keys():
				EA = rec.info['EA']
				if EA == 'silent':
					score = 0
				elif EA == 'STOP':
					score = 1
				elif EA == 'fs-indel':
					score = 1
				else:
					continue
			score = retrieve_pph2(rec.info['dbNSFP_Polyphen2_HDIV_score'])
		for sample in samples_list:
			zyg = get_zygo(rec.samples[sample]['GT'])
			temp_df.loc[sample, rec_gene] *= (1 - score) ** zyg
	return temp_df


def main():
	vcffile = sys.argv[1]
	casesfile = sys.argv[2]
	controlsfile = sys.argv[3]
	genefile = sys.argv[4]
	ref_file = sys.argv[5]
	state = sys.argv[6]
	n_jobs = sys.argv[7]

	n_jobs = int(n_jobs)
	cases = [line.strip('\n') for line in open(casesfile)]
	conts = [line.strip('\n') for line in open(controlsfile)]
	samples = cases + conts

	genelist = [line.strip('\n') for line in open(genefile)]

	ref = pd.read_csv(ref_file, delimiter='\t', header=0, index_col='name2')
	ref = ref.loc[ref.index.intersection(genelist)]

	design_matrix = parallelize(vcf_file=vcffile, reference=ref, samples=samples, n_jobs=n_jobs, state=state)
	design_matrix.to_csv('./design_matrix_' + state + '.txt', sep='\t', header=True, index=True)


if __name__ == "__main__":
	main()
