
## version 3.0.0

- CopraRNA2.pl :
  - cop1 param and related scripts and calls removed
  - new coprarna_expert_option_file for advanced and experimental parameters
  - new hybrid_coverage option to filter for full hybrids (anti-sense)
  - IntaRNA-specific parameters now set via default intarna_option_file
  

### 191024 Martin Raden

- CopraRNA2-deps.yml
  - gzip removed (obsolete)
- CopraRNA2html.r
  - ensure num to plot does not exceed available rows in output
- copraRNA2_conservation_heatmaps.r
  - fix max_cores 
  - fix genename check
- copraRNA2_find_conserved_sites.r
  - send crazy dialign-tx STDOUT output to hell (/dev/null)
- renamed CopraRNA2_subprocess.err to CopraRNA2_subprocess.oe
  - removed output to CopraRNA2_subprocess.out 
  - redirecting STDOUT to STDERR where not captured explicitely

### 191023 Martin Raden

- new file 
  - copraRNA_cleanupWD.sh: dedicated script to cleanup a CopraRNA job's working directory
- CopraRNA2.pl
  - cleanup removed with calling copraRNA_cleanupWD.sh
- adding 'suppressPackageStartupMessages' to all lib loadings in R scripts
- homology_intaRNA.pl
  - DAVID enrichment visualization done only if enrichment data present
- CopraRNA2-deps.yml
  - r-base >= 3.6.0 (was ==)
- obsolete 'edit' comments removed

### Jens Georg and Martin Raden

- removed files
  - join_pvals_coprarna2.R 
  - evo_heatmap.R 
  - copraRNA2_position_script_for_evo_precalculated_alignments_w_ooi.R 
- script_R_plots_8.R
  - ooi column changed to 4 (was 3)
- prepare_output_for_websrv_new.pl 
  - ooi column changed to 4 (was 3)
- homology_intaRNA.pl 
  - log errors to less files
  - refactoring implementing dedicated functions to reduce code redudancy 
- new files
  - CopraRNA2-deps.yml : holds all dependencies available from conda
  - CopraRNA2html.r : offline html summary
- README.md 
  - docu of new parameters
    - genomePath : optional path where genome files are stored (and loaded from)
    - temperature
  - docu of dep install via conda environment
- STDERR and STDOUT output of subscripts now mainly stored in CopraRNA2_subprocess.err|out
- DAVIDWebService_CopraRNA.py 
  - logging (import) disabled
  - email address changed
- tests added
  - test1.sh* : test based on shorted dummy genomes files
- get_CDS_from_gbk.pl 
  - avoid duplicated CDS output (take first occurrence only, assuming it is the full gene)
- get_refseq_from_refid.pl 
  - unified NCBI access information
- cluster_intarna_csv.pl 
  - reuse of parsed CSV header instead of hard coded new header
- prepare_intarna_out.pl 
  - forward temperature argument to IntaRNA call
  - sort IntaRNA output by energy via call argument rather additional command line call
- CopraRNA2.pl
  - new arguments
    - temperature
    - genomePath
  - tool version via variable
- now gezipped genome files are stored and processed
- blastall output is gzipped

### 191008 Jens Georg

New scripts for
- p-Value combination
- conserved site detection
- conservation heatmap visualization

- homology_intaRNA.pl 
  - call refine_clustertab.r
  - replace call copraRNA2_position_script_for_evo_precalculated_alignments_w_ooi.R with copraRNA2_phylogenetic_sorting.r 
  - special $cop1 figure handling removed
  - replace call evo_heatmap.R with copraRNA2_find_conserved_sites.r and copraRNA2_conservation_heatmaps.r
- new files
  - copraRNA2_conservation_heatmaps.r 
  - copraRNA2_find_conserved_sites.r 
  - copraRNA2_phylogenetic_sorting.r 
  - dialign_conf/BLOSUM.diag_prob_t10 
  - dialign_conf/BLOSUM.scr 
  - dialign_conf/BLOSUM75.diag_prob_t2 
  - dialign_conf/BLOSUM75.scr 
  - dialign_conf/BLOSUM90.scr 
  - dialign_conf/dna_diag_prob_100_exp_110000 
  - dialign_conf/dna_diag_prob_100_exp_220000 
  - dialign_conf/dna_diag_prob_100_exp_330000 
  - dialign_conf/dna_diag_prob_100_exp_550000 
  - dialign_conf/dna_diag_prob_150_exp_110000 
  - dialign_conf/dna_diag_prob_200_exp_110000 
  - dialign_conf/dna_diag_prob_250_exp_110000 
  - dialign_conf/dna_matrix.scr 
  - join_pvals_coprarna_2.r 
  - refine_clustertab.r 
- extract_functional_enriched.R 
  - ooi column changed to 4 (was 3)
- prepare_intarna_out.pl 
  - IntaRNA call extended

## changes before 190101 not documented here


## version 2 (copied from CopraRNA2.pl)


v2.1.3   : input sRNA file is kept, updated organisms-list

v2.1.2   : added R downstream ooi false positive removal

v2.1.1   : added input exceptions
           DAVID python code now py2 and py3 compatible
           changed coloring in evolutionary heatmap
           fixed issue for regions plots for sequences with non ATGC alphabet

v2.1.0   : topcount default 200
           stopped tracking edits manually

v2.0.6   : new p-value combination (no hard cutoff anymore // switching to integrated superior method) 
           standard root function for weights is now 1 instead of 2.5 (no more weight reduction)
           added phylogeny output directory in clean output
           -cop2 option is now -cop1 since CopraRNA2 will be standard
           removed -pvcut option
           added -nooi option
           added gawk, sed, grep, tr, sort as dependencies
           added extended regions plots
           aux enrichment now also done for count specified by $enrich
           adding consensus prediction
           removed phantomjs bin // now installing it via conda

v2.0.5.1 : major restructuring due to changed IntaRNA version (2.0.4)
           added IntaRNA --tAccW and --tAccL as parameters to CopraRNA 
           adjusted update_kegg2refseq for new format of prokaryotes.txt
           added verbose terminal printing option
           added topcount option
           added pvalue cutoff option for CopraRNA 2
           now using Cwd 'abs_path' to make script path locations dynamic
           added warning for run directories that contain many files
           added websrv option to only output websrv files if explicitly called for
           added root option // applies this root function to weights both for CopraRNA1 and CopraRNA2 pvalue combination
           now calculating normalized IntaRNA energy scores internally in IntaRNA // adjusted CopraRNA accordingly
           added enrichment parameter 
           added noclean parameter

v2.0.5   : changed to using IntaRNA2.0 ui
           local mirror for .gbk changed to .gb because file ending in local mirror changed
           removed evir dependency by now sourcing the gev and pgev functions locally
           added MIT license
v2.0.4.2 : now using Hartung.R to combine pvalues // no previous calculation of rho needed
v2.0.4.1 : added new adjust_tags_clustered.R that now makes two files one for pold one for padj // also added rho calculation rho_script.R
           changed scale_clusters.pl to combine_clusters.pl // reimplementation
           weights are no longer subjected to the 2.5th root (only for CopraRNA2)
v2.0.4   : major changes to the code preparing for CopraRNA2 benchmarking and publication (UPGMA tree for example)
v2.0.3.2 : fixed issue with RefSeq IDs longer than 11 that chars caused job fail // changed DAVID parameters to "medium" from the webpage
           mirroring DAVID v6.7 specifically instead of DAVID v6.8
v2.0.3.1 : changed DAVID-WS from perl to python client
v2.0.3   : using local mirrors of old NCBI ID system for compatibility if available
v2.0.2   : support of new NCBI ID system
v2.0.1   : Iterative organism subset analysis enabled. Auxiliary enrichment output added. Minimal relative cluster size parameter added. IntaRNA parameters changed to -w 150 -L 100

## version 1

v1.3.0   : Potential outlier detection; evolutionary tree visualization; minor bugfix in weight calculation.
v1.2.9   : Now using (Benjamini&Hochberg, 1995) for false discovery rate (fdr) estimation. Fixed issue where trees with branch lengths of zero would cause job failures.
v1.2.8   : Fixed the issue where jobs with input organisms with exactly the same 16S sequences would fail
v1.2.7   : Reimplementation of p-value joining (runtime reduction); Minor bugfix for heatmap drawing and regions plots
v1.2.6   : Added heatmap pdf output
v1.2.5   : Added functional enrichment heatmaps
v1.2.4   : Changed DomClust parameters to standard MBGD parameters
v1.2.3   : BLAST speedup
v1.2.2   : Fixed issue with organism: 'sfd'
v1.2.1   : RefSeq files now being downloaded from NCBI FTP
