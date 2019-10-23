
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

### changes before 190101 not documented here

