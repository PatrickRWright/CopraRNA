#!/usr/bin/env perl

use strict;
use warnings;

print "
CopraRNA - v.2.0.6  - archive README

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

CopraRNA1_final.csv:
This is the final result table of your CopraRNA1 run.
It contains the the amount of predictions specified by the
topcount parameter (def: 100). If less than this amount
of predictions is returned then the the file contains
all possible predictions.
 
CopraRNA2_final.csv: (only when -cop2 was set)
This is the final result table of your CopraRNA2 run.
t contains the the amount of predictions specified by the
topcount parameter (def: 100). If less than this amount
of predictions is returned then the the file contains
all possible predictions.

CopraRNA1_final_all.csv:
This is the non truncated version of the CopraRNA1 result.
It contains all possible CopraRNA1 predictions.

CopraRNA2_final_all.csv: (only if -cop2 was set)
s is the non truncated version of the CopraRNA2 result.
It contains all possible CopraRNA2 predictions.

CopraRNA_option_file.txt:
This file contains the set of input options.

16s_sequences.fa (in FASTA):
This is the FASTA file containing the 16s sequences that
the phylogeny is calculated on.

compatible.treefile:
This file stores the 16s phylogeny.

zscore.weight/weights.warning:
The unrooted weights for the individual organisms are stored in this file.
If an organism has a weight greater than 0.5 then this is reported in
weights.warning.

cluster.tab:
This tab delimited file contains the clusters of homologous 
protein coding genes throughout the organisms participating 
in the prediction. The clusters are computed with DomClust.

*regions *pdf *ps *png:
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

*tags.clustered*:
These files contain the IntaRNA predictions for for the clusters
of homologous putative targets. _rcsize is the truncated input
for the CopraRNA1 thread and _trunc is the truncated input
for the CopraRNA2 (if -cop2) thread.

target_sequences_orgofint.fa: (only if -websrv was set // in FASTA)
This FASTA contains the putative target sequences extracted from the
organism of interest's genome. It is a copy of one of the
*_upfromstartpos_*_down_*.fa files. 

*_websrv_table.csv: (only if -websrv was set)
These are csv tables for the Freiburg RNA tools webserver frontend.

termClusterReport_cop1.txt: (only if -enrich was specified // in Enrichment)
This file contains the DAVID functional enrichment for the amount of
CopraRNA1 candidates specified in the -enrich parameter.

org_of_interest_aux_enrichment.txt: (in Enrichment)
This is the auxilliary enrichment file for the organism of interest.

copra_heatmap.html/copraRNA.json: (in Enrichment)
This html file contains the heatmap for the enriched terms from your prediction.
It can be viewed in your web browser. The json file is needed for correct
display of copra_heatmap.html.

enriched_heatmap_big_cop1.*: (in Enrichment)
pdf and png files for the functional enrichment heatmap.

README.txt:
This file.

";
