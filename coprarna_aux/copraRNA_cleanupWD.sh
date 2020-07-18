#!/usr/bin/env bash

###################################################################
# cleans a CopraRNA job directly after the job was run, i.e.
# - removes temporary files
# - aggregates input and output files in respective subfolders
###################################################################

rm -f *.gb *.gb.gz;
rm -f compatible.distmat.mapped;
rm -f *anno* padj.csv;
rm -f *ncRNA*;
rm -f *pvalues*;
rm -f *.fa.intarna.sorted.csv *opt.intarna.csv;
rm -f gene_CDS_exception.txt find_gaps.txt;
rm -f merged_refseq_ids.txt;    
rm -f CopraRNA2_prep*;
rm -f fasta_temp_file fasta_temp_file_out;
find -regex ".*fa[0-9]+$" | xargs rm -f;
rm -f ncrna_aligned.fa;
rm -f weights.warning;
rm -f aligned_sRNA.fa;
rm -rf target_alignments;
rm -f cluster_backup.tab;
rm -f opt_tags.clustered;
rm -f zscore.weight;
rm -f weights.txt;
rm -f CopraRNA_result.map_evo_align;
rm -f markdown_final.Rmd;

# blast files
rm -f all.fas*;
rm -f duplicated_CDS.txt;

# enrichment
rm -f enrichment.txt *DAVID* IntaRNA* intarna*;
rm -f *IntaRNA1_ui* *top_targets*;

############### subdir distribution ######################

# local bash function to move files
mv2dir () {
	# check if file exists: if so create folder and move
	[ ! -f $1 ] || (mkdir -p $2; mv $1 $2);
}


# log files for debug  ###########################
for f in \
CopraRNA2_subprocess.* \
err.log \
error.log \
formatdb.log \
; do
	mv2dir $f log;
done

# IntaRNA output files  ###########################
mv2dir *.fa.intarna.csv IntaRNA;

# Phylogeny files  ###########################
for f in \
16s_sequences* \
; do
	mv2dir $f Phylogeny;
done


# FASTA files  ###########################
mv2dir *.fa FASTA;
[ -f input_sRNA.fa ] || mv -f FASTA/input_sRNA.fa .; # move input file back to root folder

# Regions_plots  ###########################
for f in \
*regions* thumbnail_* \
; do
	mv2dir $f Regions_plots;
done

# Enrichment files  ###########################
for f in \
copra_heatmap.html \
copraRNA.json \
enriched_heatmap_big* \
termClusterReport.txt \
aux_table.csv \
; do
	mv2dir $f Enrichment;
done
rm -f org_of_interest_aux_enrichment.txt;

# Rdata files  ###########################
for f in \
16S_tree.Rdata \
int_sites.Rdata \
order_table_all_orgs.Rdata \
peak_list.Rdata \
copra_results_all.Rdata \
cluster.tab \
; do
	mv2dir $f Rdata;
done

