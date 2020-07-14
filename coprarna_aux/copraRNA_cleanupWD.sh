#!/usr/bin/env bash

###################################################################
# cleans a CopraRNA job directy after the job was run, i.e.
# - removes temporary files
# - aggregates input and output files in respective subfolders
###################################################################

rm -f *.gb *.gb.gz;
rm -f compatible.distmat.mapped;
rm -f *anno* padj.csv;
rm -f *ncRNA*;
rm -f *pvalues*;
rm -f *.fa.intarna.sorted.csv *opt.intarna.csv;
rm -f gene_CDS_exception.txt find_gaps.txt distmat.out;
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

# v1 specific
rm -f *pvsample*;
rm -f CopraRNA1_final_all.csv CopraRNA1_final.csv;

############### subdir distribution ######################

# log files for debug
mkdir -p log;
mv CopraRNA2_subprocess.* err.log error.log formatdb.log log/.;

# make subdirs for IntaRNA, FASTA, Enrichment, Phylogeny and regions plots
mkdir -p IntaRNA;
mv *.fa.intarna.csv IntaRNA;

mkdir -p Phylogeny;
mv compatible.* Phylogeny;
mv 16s_sequences* Phylogeny;

mkdir -p FASTA;
mv *.fa FASTA;
mv -f FASTA/input_sRNA.fa .; # move input file back to root folder

mkdir -p Regions_plots;
mv *regions* Regions_plots;
for f in thumbnail_*; do
	mv $f Regions_plots;
done

if [ -s copra_heatmap.html ]; then 
    mkdir -p Enrichment;
    mv copra_heatmap.html Enrichment;
    mv copraRNA.json Enrichment;
    mv enriched_heatmap_big* Enrichment;
    mv termClusterReport.txt Enrichment;
    rm -f org_of_interest_aux_enrichment.txt;
    mv aux_table.csv Enrichment;
fi

# make an archive for the Rdata files
if [ -s int_sites.Rdata ]; then
    mkdir Rdata;
    mv int_sites.Rdata Rdata;
    mv order_table_all_orgs.Rdata Rdata;
    mv peak_list.Rdata Rdata;
	mv 16S_tree.Rdata Rdata;
	mv copra_results_all.Rdata Rdata;
	mv cluster.tab Rdata;	
fi
