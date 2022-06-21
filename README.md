## Introduction:
***
MPASS.pl creates a distance matrix for the construction of meta-phylogenomic trees from protein FASTA files obtained by the metagenomic shotgun sequencing.<br><br>


## Pre-requisite:
***
1. blast-legacy-2.2.26
2. phylip-3.697 or BIONJ
3. protein FASTA files  
<br>
<div style="text-align: left;">
blast-legacy and phylip can be installed using anaconda.
</div>

```vb
conda install -c bioconda blast-legacy
conda install -c bioconda phylip
```  
<br>

BIONJ can be downloaded from `http://www.atgc-montpellier.fr/bionj/binaries.php`.  
<br>

The file name of protein FASTA files should be as `*_aa`. No file extension is required.<br>
FASTA headers must include gene No., node No., coverage and length of sequences as follows:

```vb
>gene_x_NODE_y_covaaaa_lengthbbbb

(Example)
>gene_2_NODE_1_cov14.182092_length17667
```  
<br>

## Installation:
***
Download perl script MPASS.pl.
<br><br>

## Running (Example):
***
After installation, MPASS.pl can be run directly from command line.
The protein FASTA files must be in the current directory.
(E-value threshold = 10, Number of threads = 8)
```vb
perl MPASS.pl -d DISTANCElist.txt -o MPASSmatrix.txt -e 10 -t 8
```
<br><br>

## Citation:
***
<br><br>

## Recent updates:
***
(6-20-2022)
- MPASS.pl is now available.
<br><br>
