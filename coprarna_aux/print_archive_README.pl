#!/usr/bin/env perl

use strict;
use warnings;

print "
CopraRNA - v.2.0.5.1  - archive README

License: MIT

When using CopraRNA please cite:
Patrick R. Wright, et al. 
Comparative genomics boosts target prediction for bacterial small RNAs
Proc Natl Acad Sci USA, 2013, 110 (37), E3487–E3496

and/or

Patrick R. Wright, et al.
CopraRNA and IntaRNA: predicting small RNA targets, networks and interaction domains
Nucleic Acids Research, 2014, 42 (W1), W119-W123

If you run into problems or have any questions, please do not hesitate to
contact us at:
rna\@informatik.uni-freiburg.de


-- Archive file explanations

16s_sequences.fa:
This is a FASTA formatted file containing the 16s sequences that
the phylogeny is calculated on, which subsequently leads to the
weights (zscore.weight) the individual organisms attain.

cluster.tab:
This tab delimited file, contains the clusters of homologous genes
throughout the organisms participating in your analysis. The clusters
are computed with DomClust.

compatible.fneighbor:
This file contains the nieghbor-joining tree calculated on basis of
the 16s sequences.

CopraRNA_result.csv:
This is the final result table of your run. It contains the
top 100 predictions.

*regions *pdf *ps *png:
These image files contain the regions plots of your run.

*final.csv:
The final.csv files contain the results of the single IntaRNA whole
genome predictions for each organism participating in the analysis.
One of these files is created per organism.

input_sRNA.fa / ncrna.fa:
These FASTAs contain the sequences you submitted.

16s_sequences.aln:
This is the mafft alignment of the 16s sequences.

CopraRNA_result_all.csv:
This is the whole prediction table. CopraRNA_result.csv was created
from this table. However, usually predictions beyond a CopraRNA p-value
of 0.01 are not considered to be significant.

rhodevelopment.txt:
This file gives you an overview of the development of rho in the
analysis. For more information on this, please resort to the publication.

target_sequences_orgofint.fa:
This FASTA contains the putative target sequences extracted from the
organism of interest's genome. The putative target sequences for the
other organisms are not stored in the archive.

termClusterReport.txt:
This file contains the DAVID functional enrichment for the
candidates with a CopraRNA p-value <= 0.01.

zscore.weight:
The weights for the individual organisms are stored in this file.

copra_heatmap.html:
This html file contains the heatmap for the enriched terms from your prediction.
It can be viewed in your web browser.

";
