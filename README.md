# CopraRNA [![GitHub](https://img.shields.io/github/tag/PatrickRWright/CopraRNA.svg)](https://github.com/PatrickRWright/CopraRNA)  [![Bioconda](https://anaconda.org/bioconda/coprarna/badges/version.svg)](https://anaconda.org/bioconda/coprarna) [![Docker Repository on Quay](https://quay.io/repository/biocontainers/coprarna/status "Docker Repository on Quay")](https://quay.io/repository/repository/biocontainers/coprarna)
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

Please note: Version 2.1.1 is currently experimental and changes are actively being pushed.

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

In order to use CopraRNA you can either install it directly via conda or
clone this github repository and install the dependencies individually.

<br /><br />
<a name="deps" />
## Dependencies

The specified versions are tested and functional.

- bzip2 1.0.6 (for the core genome archive)                            // conda install bzip2
- gawk 4.1.3                                                           // conda install gawk
- sed 4.2.2.165-6e76-dirty                                             // conda install sed
- grep 2.14                                                            // conda install grep
- GNU coreutils 8.25                                                   // conda install coreutils 
- IntaRNA 2.1.0                                                        // conda install intarna
- EMBOSS package 6.5.7 - distmat (creates distance matix from msa)    // conda install emboss
- embassy-phylip 3.69.650 - fneighbor (creates tree from dist matrix)  // conda install embassy-phylip
- ncbiblast-2.2.22                                                     // conda install blast-legacy
- DomClust 1.2.8a                                                      // conda install domclust
- MAFFT 7.310                                                          // conda install mafft
- clustalo 1.2.3                                                       // conda install clustalo
- phantomjs 2.1.1-0                                                    // conda install phantomjs

- Perl (5.22.0) Module(s):                                             // conda install perl

    - List::MoreUtils 0.413                                                // conda install perl-list-moreutils
    - Parallel::ForkManager 1.17                                           // conda install perl-parallel-forkmanager
    - Getopt::Long 2.45                                                    // conda install perl-getopt-long
    - Bio::SeqIO (bioperl 1.6.924)                                         // conda install perl-bioperl
    - Bio::DB::EUtilities 1.75                                             // conda install perl-bio-eutilities
    - Cwd 3.56                                                             // included in the conda perl installation       

- R statistics 3.2.2                                                   // conda install r-base==3.2.2

    - seqinr 3.1\_3                                                       // conda install r-seqinr 
    - robustrankaggreg 1.1                                                // conda install r-robustrankaggreg
    - pheatmap 1.0.8                                                      // conda install r-pheatmap

- python                                                              // conda install python

    - sys                                                                  // available from conda python
    - logging                                                              // available from conda python
    - traceback                                                            // available from conda python 
    - suds.metrics (suds-jurko 0.6)                                        // conda install suds-jurko
    - suds         (suds-jurko 0.6)                                        // conda install suds-jurko
    - suds.client  (suds-jurko 0.6)                                        // conda install suds-jurko
    - datetime                                                             // available from conda python

<br /><br />
<a name="instconda" />
## CopraRNA via conda (bioconda channel)
The most easy way to locally install CopraRNA is via conda using the 
[bioconda](https://bioconda.github.io/) 
channel (linux only). This way, you will install CopraRNA along
with all dependencies.
Follow
[![install with bioconda](https://img.shields.io/badge/install%20with-bioconda-brightgreen.svg?style=flat-square)](http://bioconda.github.io/recipes/coprarna/README.html)
to get detailed information. We recommend installing into a dedicated environment, to avoid conflicts with
other installed tools. Following two commands install CopraRNA into the enviroment and activate it:
```bash
conda create -n coprarnaenv -c bioconda -c conda-forge coprarna
source activate coprarnaenv
```
<br /><br />
<a name="biocontainer" />

## Usage via biocontainer (docker)

CopraRNA can be retrieved and used as docker container with all dependencies via [docker](https://docs.docker.com/engine/installation/). Once you have docker installed simply type:
```bash
       docker run -i -t quay.io/biocontainers/coprarna:2.1.0--0 /bin/bash
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
CopraRNA2.pl -srnaseq sRNAs.fa -ntup 200 -ntdown 100 -region 5utr -enrich 200 -topcount 200 -cores 4
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
- `--rcsize` : minimum amount (%) of putative target homologs that need to be available for a target cluster 
               to be considered in the CopraRNA1 part (see --cop1) of the prediction (def:0.5)
- `--winsize`                 IntaRNA target (--tAccW) window size parameter (def:150)
- `--maxbpdist`               IntaRNA target (--tAccL) maximum base pair distance parameter (def:100)
- `--cop1`                    switch for CopraRNA1 prediction (def:off)
- `--cons`                    controls consensus prediction (def:0)
    - '0' for off
    - '1' for organism of interest based consensus
    - '2' for overall consensus based prediction
- `--verbose` : switch to print verbose output to terminal during computation (def:off)
- `--websrv` : switch to provide webserver output files (def:off)
- `--noclean` : switch to prevent removal of temporary files (def:off)
- `--enrich` : if entered then DAVID-WS functional enrichment is calculated with given amount of top predictions (def:off)
- `--nooi` : if set then the CopraRNA2 prediction mode is set not to focus on the organism of interest (def:off)
- `--root` : specifies root function to apply to the weights (def:1)
- `--topcount` : specifies the amount of top predictions to return and use for the extended regions plots (def:200)

<br /><br />
<a name="updateava" />
# Update CopraRNA available organisms

In the update_kegg2refseq directory you create a new run directory and change into this directory.
Here you can execute build_kegg2refseq.pl which will download prokaryotes.txt from the
NCBI and process it into the files CopraRNA_available_organisms.txt and kegg2refseqnew.csv.
These two files must then be copied into coprarna_aux where they override their older versions.

