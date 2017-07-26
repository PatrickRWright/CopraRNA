#!/usr/bin/env perl

use strict;
use warnings;

use Getopt::Long;
use Cwd 'abs_path'; ## edit 2.0.5.1

# CopraRNA 2.0.5.1

 # License: MIT

 # When using CopraRNA please cite:
 # Patrick R. Wright, et al. 
 # Comparative genomics boosts target prediction for bacterial small RNAs
 # Proc Natl Acad Sci USA, 2013, 110 (37), E3487–E3496.

 # and/or

 # Patrick R. Wright, et al.
 # CopraRNA and IntaRNA: predicting small RNA targets, networks and interaction domains
 # Nucleic Acids Research, 2014, 42 (W1), W119-W123.

 # If you run into problems or have questions, please do not hesitate to
 # contact us.
 # rna@informatik.uni-freiburg.de

 #    Comparative prediction algorithm for sRNA targets (CopraRNA)
 #
 #    Contributions by: 
 #
 #    Patrick R. Wright              
 #    Jens Georg
 #    Martin Mann
 #    Rolf Backofen
 #    Steffen C. Lott
 #    Andreas S. Richter 
 #    Robert Kleinkauf
 #    Dragos A. Sorescu
 #    Wolfgang R. Hess
 #    Stephan Klaehn
 #    Joerg Vogel
 #    Kai Papenfort

#####################################################

#### dependencies: (specified versions are tested and functional)

### Bio Software

# IntaRNA 2.0.4                                                        // conda install intarna
# EMOBOSS package 6.5.7 - distmat (creates distance matix from msa)    // conda install emboss
# embassy-phylip 3.69.650 - fneighbor (creates from distance matrix)   // conda install embassy-phylip
# ncbiblast-2.2.22                                                     // conda install blast-legacy
# clustalw 2.1                                                         // TODO remove this dependency in regions plot script jens
# DomClust 1.2.8a                                                      // conda install domclust
# MAFFT 7.310                                                          // conda install mafft

### Perl (5.22.0) Module(s):                                           // perl via conda install perl

# List::MoreUtils 0.413                                                // conda install perl-list-moreutils
# Parallel::ForkManager 1.17                                           // conda install perl-parallel-forkmanager
# Getopt::Long 2.45                                                    // conda install perl-getopt-long
# Bio::SeqIO (bioperl 1.6.924)                                         // conda install perl-bioperl
# Bio::DB::EUtilities 1.75                                             // conda install perl-bio-eutilities
# Cwd 3.56                                                             // included in the conda perl installation	

### R Package(s):

# R statistics 3.2.2                                                   // conda install r-base==3.2.2
# seqinr 3.1_3                                                         // conda install r-seqinr 

### python 2.7.13                                                      // conda install python==2.7.13

## packages:
# sys                                                                  // available from conda python (2.7.13)
# logging                                                              // available from conda python (2.7.13)
# traceback                                                            // available from conda python (2.7.13) 
# suds.metrics (suds-jurko 0.6)                                        // conda install suds-jurko
# suds         (suds-jurko 0.6)                                        // conda install suds-jurko
# suds.client  (suds-jurko 0.6)                                        // conda install suds-jurko
# datetime                                                             // available from conda python (2.7.13)

#####################################################

#### changelog

# v2.0.5.1 : major restructuring due to changed IntaRNA version (2.0.4)
#            added IntaRNA --tAccW and --tAccL as parameters to CopraRNA 
#            adjusted update_kegg2refseq for new format of prokaryotes.txt
#            added verbose terminal printing option
#            added topcount option
#            added pvalue cutoff option for CopraRNA 2
#            now using Cwd 'abs_path' to make script path locations dynamic
#            added warning for run directories that contain many files
#            added websrv option to only output websrv files if explicitly called for
#            added root option // applies this root function to weights both for CopraRNA1 and CopraRNA2 pvalue combination
#            now calculating normalized IntaRNA energy scores internally in IntaRNA // adjusted CopraRNA accordingly
#            added enrichment parameter 
#            added noclean parameter
#
# v2.0.5   : changed to using IntaRNA2.0 ui
#            local mirror for .gbk changed to .gb because file ending in local mirror changed
#            removed evir dependency by now sourcing the gev and pgev functions locally
#            added MIT license
# v2.0.4.2 : now using Hartung.R to combine pvalues // no previous calculation of rho needed
# v2.0.4.1 : added new adjust_tags_clustered.R that now makes two files one for pold one for padj // also added rho calculation rho_script.R
#            changed scale_clusters.pl to combine_clusters.pl // reimplementation
#            weights are no longer subjected to the 2.5th root (only for CopraRNA2)
# v2.0.4   : major changes to the code preparing for CopraRNA2 benchmarking and publication (UPGMA tree for example)
# v2.0.3.2 : fixed issue with RefSeq IDs longer than 11 that chars caused job fail // changed DAVID parameters to "medium" from the webpage
#            mirroring DAVID v6.7 specifically instead of DAVID v6.8
# v2.0.3.1 : changed DAVID-WS from perl to python client
# v2.0.3   : using local mirrors of old NCBI ID system for compatibility if available
# v2.0.2   : support of new NCBI ID system
# v2.0.1   : Iterative organism subset analysis enabled. Auxiliary enrichment output added. Minimal relative cluster size parameter added. IntaRNA parameters changed to -w 150 -L 100
# v1.3.0   : Potential outlier detection; evolutionary tree visualization; minor bugfix in weight calculation.
# v1.2.9   : Now using (Benjamini&Hochberg, 1995) for false discovery rate (fdr) estimation. Fixed issue where trees with branch lengths of zero would cause job failures.
# v1.2.8   : Fixed the issue where jobs with input organisms with exactly the same 16S sequences would fail
# v1.2.7   : Reimplementation of p-value joining (runtime reduction); Minor bugfix for heatmap drawing and regions plots
# v1.2.6   : Added heatmap pdf output
# v1.2.5   : Added functional enrichment heatmaps
# v1.2.4   : Changed DomClust parameters to standard MBGD parameters
# v1.2.3   : BLAST speedup
# v1.2.2   : Fixed issue with organism: 'sfd'
# v1.2.1   : RefSeq files now being downloaded from NCBI FTP

my $help = ""; ## edit 2.0.4
my $sRNAs_fasta = "input_sRNA.fa";
my $upstream = 200;
my $downstream = 100;
my $region = "5utr";
my $RelClusterSize = 0.5;
my $core_count = 1; # how many parallel processes are allowed
my $winsize = 150; # IntaRNA window size
my $maxbpdist = 100; # IntaRNA maximum base pair distance 
my $cop2 = 0;
my $verbose = 0; ## edit 2.0.5.1
my $noclean = 0; ## edit 2.0.5.1
my $websrv = 0; ## edit 2.0.5.1
my $pvalcutoff = 0.15; # p-value cutoff for CopraRNA 2 // ## edit 2.0.5.1
my $topcount = 100; # amount of top predictions // ## edit 2.0.5.1
my $root = 2.5; # root function to apply to the weights // ## edit 2.0.5.1
my $enrich = 0; ## edit 2.0.5.1 // functional enrichment needs to be specifically turned on 
                ##              // this option also allows to specify how many top predictions to use for the enrichment

# get absolute path
my $ABS_PATH = abs_path($0); ## edit 2.0.5.1
# remove script name at the end
# match all non slash characters at the end of the string
$ABS_PATH =~ s|[^/]+$||g; ## edit 2.0.5.1
my $PATH_COPRA = $ABS_PATH; ## edit 2.0.4 // changed $versiondirectory to $PATH_COPRA in all scripts

GetOptions ( ## edit 2.0.4
    'help|?'		=> \$help,
    'srnaseq:s'		=> \$sRNAs_fasta,
    'ntup:i'		=> \$upstream,
    'ntdown:i'		=> \$downstream,
    'cores:i'		=> \$core_count,
    'region:s'		=> \$region, # one of "5utr", "3utr", "cds"
    'rcsize:f'		=> \$RelClusterSize,
    'winsize:i'		=> \$winsize,    ## edit 2.0.5.1 
    'maxbpdist:i'	=> \$maxbpdist,  ## edit 2.0.5.1
    'cop2'		=> \$cop2, # switch for coprarna2, if set then coprarna1 and 2 are run // else only coprarna1
    'verbose'		=> \$verbose, # switch for verbose output during computation
    'websrv'		=> \$websrv, # switch for providing webserver output
    'noclean'		=> \$noclean, # switch to prevent cleaning of files
    'pvcut:f'		=> \$pvalcutoff, # p-value cutoff for CopraRNA 2
    'topcount:i'	=> \$topcount, # amount of top predictions to return ## edit 2.0.5.1
    'enrich:i'		=> \$enrich, # functional enrichment needs to be specifically turned on // also how many top preds to use for enrichment 
    'root:i'		=> \$root, # root function to apply to the weights ## edit 2.0.5.1
);

# TODO:
# - think about enrichment for CopraRNA2 output // also chartreport... second aux enrichment
# - switch nocop1
# - do manual testing
# - make a micro archive of model organisms (E. coli, PCC6803, Bacillus subtilis, Salmonella, Staphylococcus areus, Rhizobia (Agrobacterium and meliloti), Vibrio 
#   supply compressed files // make an option to check that archive
# - replace clustalw in regions plots - also make density plot discrete...

if ($help) { ## edit 2.0.4 // added  help and getopt

print "\nCopraRNA 2.0.5.1\n\n",

"CopraRNA is a tool for sRNA target prediction. It computes whole genome target predictions\n",
"by combination of distinct whole genome IntaRNA predictions. As input, CopraRNA requires\n",
"at least 3 homologous sRNA sequences from 3 distinct organisms in FASTA format.\n", 
"Furthermore, each organisms' genome has to be part of the NCBI Reference Sequence (RefSeq)\n",
"database (i.e. it should have exactly this NZ_* or this NC_XXXXXX format where * stands\n",
"for any character and X stands for a digit between 0 and 9). Depending on sequence length\n",
"(target and sRNA), amount of input organisms and genome sizes, CopraRNA can take up to 24h\n",
"or longer to compute. In most cases it is significantly faster. The central result tables\n",
"are CopraRNA1_final.csv and CopraRNA2_final.csv (if --cop2 is set). Further explanations\n",
"concerning the files in the run directory can be found in README.txt.\n\n",

"It is suggested to run CopraRNA in a dedicated empty directory to avoid unexpected behaviour.\n\n",

"The following options are available:\n\n",
" --help                    this help\n\n",
" --srnaseq                 FASTA file with small RNA sequences (def:input_sRNA.fa)\n",
" --region                  region to scan in whole genome target prediction (def:5utr)\n",
"                           '5utr' for start codon\n",
"                           '3utr' for stop codon\n",
"                           'cds' for entire transcript\n",
" --ntup                    amount of nucleotides upstream of '--region' to parse for targeting (def:200)\n",
" --ntdown                  amount of nucleotides downstream of '--region' to parse for targeting (def:100)\n",
" --cores                   amount of cores to use for parallel computation (def:1)\n",
" --rcsize                  minumum amount (%) of putative target homologs that need to be available \n",
"                           for a target cluster to be considered in the CopraRNA1 part of the prediction (def:0.5)\n",  ## edit 2.0.5.1
" --winsize                 IntaRNA target (--tAccW) window size parameter (def:150)\n",                                 ## edit 2.0.5.1
" --maxbpdist               IntaRNA target (--tAccL) maximum base pair distance parameter (def:100)\n",
" --cop2                    switch for CopraRNA2 prediction (def:off)\n",
" --verbose                 switch to print verbose output to terminal during computation (def:off)\n",  ## edit 2.0.5.1
" --websrv                  switch to provide webserver output files (def:off)\n",  ## edit 2.0.5.1
" --noclean                 switch to prevent removal of temporary files (def:off)\n",  ## edit 2.0.5.1
" --enrich                  if entered then DAVID-WS functional enrichment is calculated with given amount of top predictions (def:off)\n",  ## edit 2.0.5.1
" --pvcut                   specifies the p-values to remove before joined p-value computation (def:0.15)\n",
" --root                    specifies root function to apply to the weights (def:2.5)\n",
" --topcount                specifies the amount of top predictions to return (def:100)\n\n", ## edit 2.0.5.1

"Example call: ./CopraRNA2.pl -srnaseq sRNAs.fa -ntup 200 -ntdown 100 -region 5utr -rcsize 0.5 -winsize 150 -maxbpdist 100 -cop2 -enrich 100 -pvcut 0.15 -topcount 100 -cores 4\n\n",
"License: MIT\n\n",
"References: \n",
"1. Wright PR et al., Comparative genomics boosts target prediction for bacterial small RNAs\n   Proc Natl Acad Sci USA, 2013, 110(37), E3487–E3496\n",
"2. Wright PR et al., CopraRNA and IntaRNA: predicting small RNA targets, networks and interaction domains\n   Nucleic Acids Research, 2014, 42 (W1), W119-W123\n",
"\n";

exit(-1);

}

# input check
unless (-e $sRNAs_fasta) { ## edit 2.0.2
    die("\nError: No input FASTA supplied!\nUse '-h' option for help.\n\n");
}

# create warning for non empty run dir
my @dir_files = <*>; ## edit 2.0.5.1
my $file_count = scalar(@dir_files); ## edit 2.0.5.1
print "\nWarning: your run directory contains many files ($file_count). In general it is suggested to run CopraRNA in an empty directory to avoid unexpected behaviour.\n\n" if ($file_count > 10); ## edit 2.0.5.1

# check core count // if wrong set it to 1 // ## edit 2.0.4
if ($core_count <= 1) {
    $core_count = 1;
}

# check region parameter ## edit 2.0.4
die ("\nError: -region parameter must be one of 5utr, 3utr or cds. You set '$region'.\n\n") unless ($region eq "5utr" or $region eq "3utr" or $region eq "cds");

# check for gaps
system "grep '-' $sRNAs_fasta > find_gaps.txt"; ## edit 2.0.2
if (-s "find_gaps.txt") { die("\nError: Gaps are present in sRNA sequences. Please delete them from the file and restart.\n\n"); } ## edit 2.0.2

# check for correct RefSeq formatted headers and their presence in the availibility table ## edit 2.0.4
my $headerIDs = `grep ">" $sRNAs_fasta | sed 's/>//g' | tr '\n' ';'`;
chop $headerIDs;
my @splitHeaderIDs = split(/;/,$headerIDs);
foreach(@splitHeaderIDs) {
    die("\nError: $_ does not match correct RefSeq ID format (NZ_* or NC_XXXXXX where * stands for any character and X stands for a digit between 0 and 9).\n\n") unless ($_ =~ m/NC_\d{6}|NZ_.*/) ; 
    my $availabilityCheck = `grep '$_' $PATH_COPRA/coprarna_aux/kegg2refseqnew.csv`;
    die("\nError: $_ is not present in the availability list and is thus not compatible with CopraRNA.\n\n") unless (length $availabilityCheck); 
}

# check that maxbpdist ist smaller or equal to windowsize
die("\nError: The maximal basepair distance ($maxbpdist) is larger than the given window size ($winsize) but must be <= to the windows size. Please change the parameters accordingly.\n\n") if ($maxbpdist > $winsize); ## edit 2.0.5.1

# check for ntup + ntdown being >= $winsize because of IntaRNA parameters ## edit 2.0.5.1
my $sumUpDown = $upstream+$downstream;
if ($region eq "5utr" or $region eq "3utr") { ## edit 2.0.4.2
    die("\nError: (-ntup + -ntdown) is $sumUpDown but must be >= $winsize. Please change the parameters accordingly.\n\n") if ( $sumUpDown < $winsize ); ## edit 2.0.5.1
}

# check for correct range of 0.5 <= -rcsize <= 1.0 ## edit 2.0.4
die("\nError: -rcsize can only be specified between 0.5 and 1.0. You set '$RelClusterSize'.\n\n") unless ($RelClusterSize >= 0.5 and $RelClusterSize <= 1.0);

# check for duplicate IDs in FASTA header // not allowed ## edit 2.0.4
my $duplicate_fasta_header = `grep ">" $sRNAs_fasta | sort | uniq -d`;
chomp $duplicate_fasta_header;
die("\nError: Duplicate organisms ($duplicate_fasta_header) are present in $sRNAs_fasta.\n\n") if($duplicate_fasta_header);

# write input options to file ## edit 2.0.4
open WRITETOOPTIONS, ">", "CopraRNA_option_file.txt";
    print WRITETOOPTIONS "sRNA FASTA:" . $sRNAs_fasta . "\n";
    print WRITETOOPTIONS "nt upstream:" . $upstream . "\n";
    print WRITETOOPTIONS "nt downstream:" . $downstream . "\n";
    print WRITETOOPTIONS "region:" . $region . "\n";
    print WRITETOOPTIONS "relative clustersize:" . $RelClusterSize . "\n"; 
    print WRITETOOPTIONS "core count:" . $core_count . "\n";
    print WRITETOOPTIONS "win size:" . $winsize . "\n";
    print WRITETOOPTIONS "max bp dist:" . $maxbpdist . "\n";
    print WRITETOOPTIONS "CopraRNA2:" . $cop2 . "\n";
    print WRITETOOPTIONS "verbose:" . $verbose . "\n";
    print WRITETOOPTIONS "websrv:" . $websrv . "\n";
    print WRITETOOPTIONS "p-value cutoff:" . $pvalcutoff . "\n";
    print WRITETOOPTIONS "top count:" . $topcount . "\n";
    print WRITETOOPTIONS "root:" . $root . "\n";
    print WRITETOOPTIONS "enrich:" . $enrich . "\n";
    print WRITETOOPTIONS "noclean:" . $noclean . "\n";
    print WRITETOOPTIONS "version:CopraRNA 2.0.5.1\n";  ## edit 2.0.4.2
close WRITETOOPTIONS;
# end write options

# CopraRNA error log
system "touch err.log";

system "cp $sRNAs_fasta input_sRNA.fa" unless ($sRNAs_fasta eq "input_sRNA.fa");
$sRNAs_fasta = "input_sRNA.fa";

# format sRNA fasta - put sequence in 1 line
system $PATH_COPRA . "coprarna_aux/format_fasta.pl $sRNAs_fasta" . " > $sRNAs_fasta.temp";
system "mv $sRNAs_fasta.temp $sRNAs_fasta";

# build RefSeq input based on the sRNA input fasta (can only contain refseq IDs in header)
my $RefSeqIds = `grep '>' $sRNAs_fasta | sed 's/>//g' | tr '\n' ' '`;

# run homology_intaRNA.pl 
print $PATH_COPRA . "coprarna_aux/homology_intaRNA.pl $sRNAs_fasta $upstream $downstream $region $RefSeqIds\n" if ($verbose);
system $PATH_COPRA . "coprarna_aux/homology_intaRNA.pl $sRNAs_fasta $upstream $downstream $region $RefSeqIds";

# get organism of interest
my $ncrnaRIDs = `grep ">" ncrna.fa | sed 's/>ncRNA_//g' | tr '\n' ' '`;
my @splitRID = split(/\s/, $ncrnaRIDs);
my $organismOfInterest = $splitRID[0];
chomp $organismOfInterest;
# final optimal IntaRNA result for organism of interest
my $MainFinalCSV = $organismOfInterest . "_upfromstartpos_" . $upstream . "_down_" . $downstream . "_opt.intarna.csv"; ## edit 2.0.5.1

if ($enrich) { ## edit 2.0.5.1
    # add IntaRNA single organisms chart reports for aux enrichment // sort by p-value
    system "env LC_ALL=C sort -t';' -g -k36 $MainFinalCSV -o intarna_websrv_table.csv";
    system $PATH_COPRA . "coprarna_aux/add_GI_genename_annotation_intarna.pl";
    system $PATH_COPRA . "coprarna_aux/DAVIDWebService_IntaRNA_chartReport.py intarna_websrv_table_ncbi.csv > IntaRNA_chartReport.txt"; ## edit 2.0.3.1
    system "grep -P 'geneIds\\s=|termName\\s=' IntaRNA_chartReport.txt | sed 's/\\s//g' | sed 's/\"//g' > IntaRNA_chartReport_grepped.txt"; ## edit 2.0.3.1
    system $PATH_COPRA . "coprarna_aux/find_single_specific_targets_in_termCluster.pl > org_of_interest_aux_enrichment.txt";
}

# output warnings
system "awk -F ';' '{if (\$2 > 0.5) { print toupper(\$1) \" may be overweighted. It has weight\"; printf(\"\%.2f\", \$2); print \". You should consider checking the 16S rDNA tree. We suggest removal of outliers from your input and restarting.\";} }' zscore.weight | tr '\n' ' ' > weights.warning"; ## edit 2.0.2

# check err.log // should be empty // err.log contains information on 
# 1. not correctly downloaded RefSeq files 
# 2. gene no CDS issue 
# 3. wrong 16S counts
# 4. empty CopraRNA1_anno_addhomologs_padj_amountsamp.csv
# 5. the exception in add_pval_to_csv_evdfit.R
# -s  File has nonzero size (returns size in bytes).
if (-s "err.log") { die("\nError: CopraRNA failed. Check err.log for details.\n\n"); } ## edit 2.0.4 // added another check here at the bottom maybe we need some more check hooks in the new scripts

# move full result files 
system "mv CopraRNA1_anno_addhomologs_padj_amountsamp.csv CopraRNA1_final_all.csv";
system "mv CopraRNA2_anno_addhomologs_padj_amountsamp.csv CopraRNA2_final_all.csv" if ($cop2);

# clean up
unless ($noclean) {

    system "rm *.gb 16s_sequences.aln *pvsample* enrichment_cop1.txt";
    system "rm compatible.fneighbor compatible.distmat.mapped compatible.distmat";
    system "rm err.log *IntaRNA1_ui* *anno* padj.csv";
    system "rm *top_targets* *pvalues* ncRNA_* rhodevelopment.txt";
    system "rm *.fa.intarna.sorted.csv *opt.intarna.csv";
    system "rm gene_CDS_exception.txt find_gaps.txt distmat.out";
    system "rm -r weight_permutations";
    system "rm *DAVID* input_sRNA.fa IntaRNA* intarna* merged_refseq_ids.txt";    

    # fix warning "rm: missing operand Try 'rm --help' for more information." ## edit 2.0.1
    my $temp_fasta_check = `find -regex ".*fa[0-9]+\$"`;
    if ($temp_fasta_check) {
        system 'find -regex ".*fa[0-9]+$" | xargs rm';
    }

    # blast files
    system "rm all.fas" if (-e "all.fas");
    system "rm all.fas.blast" if (-e "all.fas.blast");
    system "rm all.fas.hom" if (-e "all.fas.hom");
    system "rm all.fas.pin" if (-e "all.fas.pin");
    system "rm all.fas.phr" if (-e "all.fas.phr");
    system "rm all.fas.tit" if (-e "all.fas.tit");
    system "rm all.fas.psq" if (-e "all.fas.psq");
    system "rm all.fas.gene" if (-e "all.fas.gene");
    system "rm error.log" if (-e "error.log");
    system "rm formatdb.log" if (-e "formatdb.log");
    system "rm N_chars_in_CDS.txt" if (-e "N_chars_in_CDS.txt");

}


