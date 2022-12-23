## Introduction:
***
MPASS.pl creates a distance matrix for the construction of metaphylogenomic trees from the metagenomic shotgun sequencing data.<br><br><br><br><br>


## Pre-requisite:
***
[MPASS.pl]
1. Paired-end fastq files (not less than 3 pairs)
2. seqtk (version 1.3-r106)
3. metaSPAdes (SPAdes genome assembler v3.15.2)
4. MetaGeneMark (GeneMark.hmm version 3.25)  
5. Configuration file for MetaGeneMark (MetaGeneMark_v1.mod)
6. blast-legacy-2.2.26
7. BIONJ, phylip-3.697 or other NJ tree construction softwares
<br>

[MPASS_core.pl]
1. Protein FASTA files
2. blast-legacy-2.2.26
3. BIONJ, phylip-3.697 or other NJ tree construction softwares
<br>

[FilteringFastq_for_MPASS.pl] (optional)
1. Trimmomatic (0.39)
2. Adapter FASTA files (optional)
<br><br>
<div style="text-align: left;">
The versions presented for all softwares are recommended, but not required.
</div>
<br><br>
<div style="text-align: left;">


seqtk, metaSPAdes, blast-legacy, phylip and Trimmomatic can be installed using anaconda.
</div>

```vb
conda install -c bioconda seqtk
conda install -c bioconda spades
conda install -c bioconda blast-legacy
conda install -c bioconda phylip
conda install -c bioconda trimmomatic
```  
<br>

MetaGeneMark and its configuration file (MetaGeneMark_v1.mod) can be downloaded from `http://exon.gatech.edu/license_download.cgi`.<br><br>
BIONJ can be downloaded from `http://www.atgc-montpellier.fr/bionj/binaries.php`.



<br><br>

When using both or either MPASS.pl and FilteringFastq_for_MPASS.pl, the paired-end FASTQ files should be named `*_1.fastq, *_2.fastq` and placed `in the current directory`.<br>
Also, the MetaGeneMark configuration file should be placed in the `current directory` (if MPASS.pl is used).
<br><br>

For the MPASS_core.pl, the FASTA files should be named `*_aa` and placed in the `current directory`. No file extension is required.<br>
The header lines of each sequence must include gene No., node No., coverage and length of sequences as follows:

```vb
>gene_x_NODE_y_covaaaa_lengthbbbb

(Example)
>gene_2_NODE_1_cov14.182092_length17667
```  
<br><br><br><br>


## Installation:
***
Download perl script MPASS.pl or MPASS_core.pl.
<br><br><br><br><br>

## Running (Example):
***
MPASS.pl can be run directly from command line.
(E-value threshold = 10, Number of threads = 8)
```vb
perl MPASS.pl -d DISTANCElist.txt -o MPASSmatrix.txt -e 10 -t 8
```
MPASS_core.pl also can be run in basically the same way:
```vb
perl MPASS_core.pl -d DISTANCElist.txt -o MPASSmatrix.txt -e 10 -t 8
```

Created distance matrix can be converted to the Newick format file by BIONJ, Neighbor (in the PHYLIP package), and other software.<br>
Using obtained Newick file, the metaphylogenomic tree can be displayed by the graphical viewer software of phylogenetic trees.<br><br><br>


If trimming and quality filtering of fastq reads are needed, you can use FilteringFastq_for_MPASS.pl to run Trimmomatic automatically. <br>
This script can trim and filter all paired-end fastq files in the current directory and can be run as follows<br>
(Minimum length of filtered sequences = 50bp, Number of threads = 8, Adapter FASTA file = TruSeq3-PE.fa):

```vb
perl FilteringFastq_for_MPASS.pl -m 50 -t 8 -a TruSeq3-PE.fa
```
<br><br>

Help for these scripts can be displayed as follows:
```vb
perl MPASS.pl -h
perl MPASS_core.pl -h
perl FilteringFastq_for_MPASS.pl -h
```

<br><br><br><br>

## Citation:
***
Satoh S, Tanaka R, Yokono M, Endoh D, Yabuki T, and Tanaka A. Phylogeny analysis of whole protein-coding genes in metagenomic data detected an environmental gradient for the microbiota. bioRxiv. https://doi.org/10.1101/2022.07.04.498637
<br><br><br><br><br>

## Recent updates:
***
(12-27-2022)
- MPASS.pl (Full version) is now available.
- Formula to calculate the metagenomic distances is changed from one that uses the simple 16S-substitution rates to one that uses the poisson-corrected 16S-substitution rates (Both MPASS.pl and MPASS_core.pl).
- Requirements, procedures for running scripts, and citation are updated.

(6-20-2022)
- MPASS_core.pl is now available.
<br><br>