# EAWavelet 

EAWavelet is a graph learning based approach to identify potential disease-gene relations using molecular network information in conjunction with patient specific mutational data. Node embedding itself is performed using 
the GraphWave algorithm pioneered by Donnat et al (https://github.com/snap-stanford/graphwave). The embedding portion of this repository is a replica of Donnat et al's repository with modifications needed to lift over from python 2 to python 3.

The output is a ranked list of genes with corresponding PCA-distances, z-scores, p-values, and corresponding FDR-values


## Setup
### Requirements

  - Linux (has not been tested on OSx)
  - Anaconda3
  - python 3.7+

### Installation
To install conda environment:
```bash
conda env create -f EAWavelet.yml
```
## Usage 
There are two main scripts for EAWavelet, vcf_to_dataframe_parallel.py and EA_wavelet analysis.py.
The first of the two processes a VCF annotated with EA or PPh2 into a data matrix of patient specific p(VIS) scores.
The second performs the actual EAWavelet analysis.

Before running the pipeline, ensure activation of the virtual enviornment:
```bash
conda activate EAWavelet
```
### VCF input processing 
Required arguments for the vcf_to_dataframe_parallel.py script are as follows:

| argument       | type          | description                                       |
|----------------|---------------|---------------------------------------------------|
| VCF file | \<file\> | Annotated VCF file path                             |
| Case IDs           | \<file\>      | text file of case sample IDs      |
| Control IDs        | \<file\>      | text file of control sample IDs |
| Gene list        | \<file\>      | Text file of genes to assess |
| Genome Reference file        | \<file\>      | genome reference file |
| VIS type        | str      | Type of variant scoring method ('EA' or 'PPh2') |
| number of jobs        | int      | number of jobs for parallelization|

Usage is as follows:
```bash
python vcf_to_dataframe_parallel.py VCFfile CaseIDs ControlIDS
                                    Genelist GenomeReference EA(or PPh2) 20
```

#### Input file Formatting
- Case ID and control ID files should be formatted as one column text files
with only sample IDs
  
- Gene list should be formatted as one column text file with desired gene names as they appear in the reference

- VCF file should follow format as described [here](<https://samtools.github.io/hts-specs/VCFv4.2.pdf>). 
 Additionally, extra information is required in the VCF as follows:
  - 'gene' annotations as necessary in the INFO column
   - Depending on variant scoring method utilized, they should also appear in the INFO column. If using 
   polyphen2, simply follow annotation instructions specified by [SnpSift](<https://pcingola.github.io/SnpEff/ss_dbnsfp/>). If using EA, then 'EA' must appear in INFO field and 
     must be typed as 'String' to allow for recognition of STOP, frameshift-indels, and SILENT mutations.
     
   - All additional INFO fields mus be defined in the header with type information
     
### EA Wavelet Analysis
Required arguments for the EA_wavelet_analysis.py script are as follows:

| argument       | type          | description                                       |
|----------------|---------------|---------------------------------------------------|
| --input | \<file\> | Data matrix output from previous step                             |
| --cases          | \<file\>      | text file of case sample IDs      |
| --conts        | \<file\>      | text file of control sample IDs |
| --network        | \<file\>      | text file for molecular network |
| --thr        | float      | FDR threshold (default 0.01) |
| --full_fdr        |       | pass this flag if full scores for all genes is desired |
| --savepath        | \<path\>      | output file name and path|

Usage is as follows:
```bash
python EA_wavelet_analysis.py --input /path/to/data.txt --cases /path/to/case.txt --conts /path/to/cont.txt
                              --network /path/to/network.csv --thr 0.1 --savepath /path/to/ouput.txt
```
Output from this command will be a text file of 5 columns (distance, gene, z-score, p-value, and fdr value). Only genes whose FDR value is less than the 
requested threshold of 0.1 will be reported. If a FULL output of all genes is desired, then --full_fdr flag should be passed as follows:
```bash
python EA_wavelet_analysis.py --input /path/to/data.txt --cases /path/to/case.txt --conts /path/to/cont.txt
                              --network /path/to/network.csv --full_fdr --savepath /path/to/ouput.txt
```
 #### Input file formats

- CaseID and ControlID files will be the same as from the previou step
- Network file should be a comma separated file of 3 columns. Each row in the file 
should represent an edge. Columns will be node1, node2, edge-weight. If no edge weight exits, simply fill with 1s.
  - For example, edge between gene x and y would appear as:
    
      - x,y,1
   
## Credits
As previously stated this method uses GraphWave node embedding algorithm designed and developed by Donnat et al:

- Donnat, C., Zitnik, M., Hallac, D., & Leskovec, J. (2018). Learning Structural Node Embeddings via Diffusion Wavelets. 1320–1329. https://doi.org/10.1145/3219819.3220025

The file graphwave_py3.py incorporates and uses portions of Donnat et al's github repository for GraphWave. Modifications
were made in order to allow for compatability with python 3.7. 
