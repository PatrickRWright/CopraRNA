# CopraRNA version 2.*
![CopraRNA](https://raw.githubusercontent.com/PatrickRWright/CopraRNA/master/copra_sRNA.jpg "CopraRNA")

**Phylogenetic target prediction for prokaryotic *trans*-acting small RNAs**

Please note: Version 2.0.6 is currently experimental and changes are actively being pushed.

For testing or ad hoc use of IntaRNA, you can use its webinterface at the

**==> [Freiburg RNA tools CopraRNA webserver](http://rna.informatik.uni-freiburg.de/CopraRNA/) <==**

## Citation
If you use CopraRNA, please cite our articles
- [Comparative genomics boosts target prediction for bacterial small RNAs](http://dx.doi.org/10.1073/pnas.1303248110)
  Patrick R. Wright, Andreas S. Richter, Kai Papenfort, Martin Mann, Jörg Vogel, Wolfgang R. Hess, Rolf Backofenb, and Jens Georg
  Proceedings of the National Academy of Sciences of the USA, 110, E3487-E3496, 2013, DOI(10.1073/pnas.1303248110).
- [CopraRNA and IntaRNA: predicting small RNA targets, networks and interaction domains](http://dx.doi.org/10.1093/nar/gku359)
  Patrick R. Wright, Jens Georg, Martin Mann, Dragos A. Sorescu, Andreas S. Richter, Steffen Lott, Robert Kleinkauf, Wolfgang R. Hess, and Rolf Backofen
  Nucleic Acids Research, 42, W119-W123, 2014, DOI(10.1093/nar/gku359).

<br /><br /><br /><br />
<a name="doc" />
# Documentation

## Overview

The following topics are covered by this documentation:

- [Installation](#install)
  - [Dependencies](#deps)
  - [IntaRNA via conda](#instconda)
  - [Cloning from github](#instgithub)
- [Usage and Parameters](#usage)

<br /><br /><br /><br />
<a name="install" />
# Installation


<br /><br />
<a name="deps" />
## Dependencies

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
- phantomjs 2.1.1                                                      // conda install phantomjs
- icu 56.1                                                             // conda install icu=56.1

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

- python 2.7.13                                                        // conda install python==2.7.13

    - sys                                                                  // available from conda python (2.7.13)
    - logging                                                              // available from conda python (2.7.13)
    - traceback                                                            // available from conda python (2.7.13) 
    - suds.metrics (suds-jurko 0.6)                                        // conda install suds-jurko
    - suds         (suds-jurko 0.6)                                        // conda install suds-jurko
    - suds.client  (suds-jurko 0.6)                                        // conda install suds-jurko
    - datetime                                                             // available from conda python (2.7.13)

<br /><br />
<a name="instconda" />
## CopraRNA via conda (bioconda channel)

TODO

<br /><br />
<a name="instgithub" />
## Cloning *Source code* from github (or downloading ZIP-file)

TODO

<br /><br /><br /><br />
<a name="usage" />
# Usage and parameters

TODO

<br /><br /><br /><br />
<a name="updateava" />
# Update CopraRNA available organisms

In the update_kegg2refseq directory you create a new run directory and change into this directory.
Here you can execute build_kegg2refseq.pl which will download prokaryotes.txt from the
NCBI and process it into the files CopraRNA_available_organisms.txt and kegg2refseqnew.csv.
These two files must then be copied into coprarna_aux where they override their older versions.

