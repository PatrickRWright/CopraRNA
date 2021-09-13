#!/usr/bin/env perl

use strict;
use warnings;

print "
CopraRNA - archive README

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

CopraRNA_result.csv:
This is the final result table of your CopraRNA run.
It contains the the amount of predictions specified by the
topcount parameter (def: 200). If less than this amount
of predictions is returned then the the file contains
all possible predictions.

CopraRNA_result_all.csv:
This is the non truncated version of the CopraRNA result.
It contains all possible CopraRNA predictions.

CopraRNA_option_file.txt:
This file contains the set of input options.

16s_sequences.fa/.aln (in Phylogeny):
This is the FASTA file containing the (aligned) 16s sequences that
the phylogeny is calculated on.

cluster.tab:
This tab delimited file contains the clusters of homologous 
protein coding genes throughout the organisms participating 
in the prediction. The clusters are computed with DomClust.

*regions *pdf *png (in Regions_plots):
These image files contain the regions plots of your run.

*.fa.intarna.csv (in IntaRNA):
These files contain the results of the single IntaRNA whole
genome predictions for each organism participating in the analysis.
One of these files is created per organism.

ncrna.fa (in FASTA):
This FASTA contains the sequences you submitted.

*_upfromstartpos_*_down_*.fa (in FASTA):
These FASTA files contain the putative target sequences of the
respective organisms in the prediction.

target_sequences_orgofint.fa: (only if -websrv was set // in FASTA)
This FASTA contains the putative target sequences extracted from the
organism of interest's genome. It is a copy of one of the
*_upfromstartpos_*_down_*.fa files. 

coprarna_websrv_table.csv: (only if -websrv was set)
This csv table contains the data for the Freiburg RNA tools webserver frontend.

CopraRNA_result.map_evo_align: (only if -websrv was set)
The file is used to post process the predictions for the webserver frontend.

termClusterReport.txt: (only if -enrich was specified // in Enrichment)
This file contains the DAVID functional enrichment for the amount of
candidates specified in the -enrich parameter.

aux_table.csv: (in Enrichment)
This is the auxilliary enrichment file for the organism of interest.

enriched_heatmap_big.*: (in Enrichment)
pdf and png files for the functional enrichment heatmap.

sRNA_conservation_heatmap.pdf:
This file shows the interaction conservation of the top predictions.

copra_html.zip:
Zip file containing an offline html file, to view the most important results and
visualizations of the CopraRNA prediction.

Rdata:
Precomputed R data environments. 

README.txt:
This file.

";
