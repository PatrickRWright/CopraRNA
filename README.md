# CopraRNA [![GitHub](https://img.shields.io/github/tag/PatrickRWright/CopraRNA.svg)](https://github.com/PatrickRWright/CopraRNA/releases)  [![Bioconda](https://anaconda.org/bioconda/coprarna/badges/version.svg)](https://anaconda.org/bioconda/coprarna) [![Docker Repository on Quay](https://quay.io/repository/biocontainers/coprarna/status "Docker Repository on Quay")](https://quay.io/repository/biocontainers/coprarna)
![CopraRNA](https://raw.githubusercontent.com/PatrickRWright/CopraRNA/master/copra_sRNA.jpg "CopraRNA")

**Phylogenetic target prediction for prokaryotic *trans*-acting small RNAs**

CopraRNA is a tool for sRNA target prediction. It computes whole genome target predictions
by combination of distinct whole genome IntaRNA predictions. As input CopraRNA requires
at least 3 homologous sRNA sequences from 3 distinct organisms in FASTA format.
Furthermore, each organisms' genome has to be part of the NCBI Reference Sequence (RefSeq)
database (i.e. it should have exactly this NZ_* or this NC_XXXXXX format where * stands
for any character and X stands for a digit between 0 and 9). Depending on sequence length
(target and sRNA), amount of input organisms and genome sizes, CopraRNA can take up to 24h
or longer to compute. In most cases it is significantly faster. It is suggested to run CopraRNA
on a machine with at least 8 GB of memory.

CopraRNA produces a lot of file I/O. It is suggested to run CopraRNA in a dedicated
empty directory to avoid unexpected behavior.

For testing or ad hoc use of CopraRNA, you can use its webinterface at the

**==> [Freiburg RNA tools CopraRNA webserver](http://rna.informatik.uni-freiburg.de/CopraRNA/) <==**

## Citation
If you use CopraRNA, please cite our articles
- [Comparative genomics boosts target prediction for bacterial small RNAs](http://dx.doi.org/10.1073/pnas.1303248110)
  Patrick R. Wright, Andreas S. Richter, Kai Papenfort, Martin Mann, Jörg Vogel, Wolfgang R. Hess, Rolf Backofenb, and Jens Georg
  Proceedings of the National Academy of Sciences of the USA, 110, E3487-E3496, 2013, DOI(10.1073/pnas.1303248110).
- [CopraRNA and IntaRNA: predicting small RNA targets, networks and interaction domains](http://dx.doi.org/10.1093/nar/gku359)
  Patrick R. Wright, Jens Georg, Martin Mann, Dragos A. Sorescu, Andreas S. Richter, Steffen Lott, Robert Kleinkauf, Wolfgang R. Hess, and Rolf Backofen
  Nucleic Acids Research, 42, W119-W123, 2014, DOI(10.1093/nar/gku359).

<br /><br />
<a name="doc" />
# Documentation

## Overview

The following topics are covered by this documentation:

- [Installation](#install)
  - [Dependencies](#deps)
  - [CopraRNA via conda](#instconda)
  - [Cloning from github](#instgithub)
- [Usage and Parameters](#usage)
- [Update CopraRNA available organisms](#updateava)

<br /><br />
<a name="install" />
# Installation

In order to use CopraRNA you can either [install it directly via conda](#instconda) or
clone this github repository and install the dependencies individually. 
It is also possible to run CopraRNA [via a provided Docker container](#biocontainer).

<br /><br />
<a name="deps" />
## Dependencies

We provide a [list of dependencies](CopraRNA-deps.yml) within the file [CopraRNA-deps.yml](CopraRNA-deps.yml).

You can easily create a [`conda` environment](https://docs.conda.io/projects/conda/en/latest/user-guide/tasks/manage-environments.html) using the following command (after installing `conda`)
```bash
# within CopraRNA project folder
# once
conda env create --file CopraRNA-deps.yml
chmod a+x CopraRNA.pl
# always once within a shell session
conda activate CopraRNA-deps
# call CopraRNA
./CopraRNA.pl -help
```

<br /><br />
<a name="instconda" />
## CopraRNA via conda (bioconda channel)
The most easy way to locally install CopraRNA is via conda using the 
[bioconda](https://bioconda.github.io/recipes/coprarna/README.html) 
channel (linux and osx only). This way, you will install CopraRNA along
with all dependencies.
Follow
[![install with bioconda](https://img.shields.io/badge/install%20with-bioconda-brightgreen.svg?style=flat-square)](http://bioconda.github.io/recipes/coprarna/README.html)
to get detailed information. We recommend installing into a dedicated environment, to avoid conflicts with
other installed tools (or their versions). Following two commands install CopraRNA into the enviroment `CopraRNA` and activate it:
```bash
conda create -n CopraRNA -c r -c conda-forge -c bioconda coprarna
conda activate CopraRNA
```
<br /><br />
<a name="biocontainer" />

## Usage via biocontainer (docker)

CopraRNA can be retrieved and used as docker container with all dependencies via [docker](https://docs.docker.com/engine/installation/). Once you have docker installed simply type (with changed version):
```bash
       docker run -i -t quay.io/biocontainers/coprarna:3.0.0--0 /bin/bash
```
<br /><br />
<a name="instgithub" />

## Cloning *Source code* from github (or downloading ZIP-file)
```bash
git clone https://github.com/PatrickRWright/CopraRNA
```
If you installed all dependencies you should be able to directly use the source.

<br /><br />
<a name="usage" />
# Usage and parameters

Example call:
```bash
CopraRNA.pl -srnaseq sRNAs.fa -ntup 200 -ntdown 100 -region 5utr -enrich 200 -topcount 200 -cores 4
```

The following options are available:

- `--help` : help
- `--srnaseq` : FASTA file with small RNA sequences (def:input_sRNA.fa)
- `--region` : region to scan in whole genome target prediction (def:5utr)
    - '5utr' for start codon
    - '3utr' for stop codon
    - 'cds' for entire transcript
- `--ntup` : amount of nucleotides upstream of '--region' to parse for targeting (def:200)
- `--ntdown` : amount of nucleotides downstream of '--region' to parse for targeting (def:100)
- `--cores` : amount of cores to use for parallel computation (def:1)
- `--verbose` : switch to print verbose output to terminal during computation (def:off)
- `--websrv` : switch to provide webserver output files (def:off)
- `--noclean` : switch to prevent removal of temporary files (def:off)
- `--enrich` : if entered then DAVID-WS functional enrichment is calculated with given amount of top predictions (def:off)
- `--root` : specifies root function to apply to the weights (def:1)
- `--topcount` : specifies the amount of top predictions to return and use for the extended regions plots (def:200)
- `--genomePath`: path where NCBI genome files (`*.gb`) are to be stored (def:`.` i.e. working directory). Set this path if you (want to) store all your genomes in a dedicated folder to be shared by different CopraRNA calls.
- `--intarnaOptions` : path for IntaRNA parameter file.
- `--CopraRNA_expert_options` : path to parameter file for CopraRNA expert or experimental options.
- `--hybrid_threshold` : interactions are removed from the CopraRNA calculations if the hybrid covers >= hybrid_threshold of the sRNA.

<br /><br />
<a name="updateava" />
# Update CopraRNA available organisms

In the update_kegg2refseq directory you create a new run directory
```bash
mkdir run
```
and change into this directory
```bash
cd run 
```
Here you can execute build_kegg2refseq.pl 
```bash
../build_kegg2refseq.pl 
```
which will download prokaryotes.txt from the
NCBI and process it into the files `CopraRNA_available_organisms.txt` and `kegg2refseqnew.csv`.
These two files must then be copied into `coprarna_aux` where they override their older versions.

Note, some genomes are known to cause problems and might have to be (manually) removed from `CopraRNA_available_organisms.txt` (or commented with `#`):

- `# NC_022528 NZ_FO203526 NC_022543 NZ_FO203527    Vibrio_nigripulchritudo` 

