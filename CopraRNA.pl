#!/usr/bin/env perl
use strict;
use warnings;
use Getopt::Long;
use Cwd 'abs_path';

my $COPRARNA_VERSION="3.0.0";

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
 #    Contributors (in lexicographic order of first name): 
 #
 #    Andreas S. Richter 
 #    Dragos A. Sorescu
 #    Fayyaz Hussain
 #    Florian Eggenhofer
 #    Jens Georg
 #    Joerg Vogel
 #    Kai Papenfort
 #    Martin Raden, nee Mann
 #    Patrick R. Wright              
 #    Robert Kleinkauf
 #    Rolf Backofen
 #    Steffen C. Lott
 #    Stephan Klaehn
 #    Wolfgang R. Hess

#####################################################

#### dependencies: see file CopraRNA2-deps.yml

#####################################################

# get absolute path
my $ABS_PATH = abs_path($0);
# remove script name at the end
# match all non slash characters at the end of the string
$ABS_PATH =~ s|[^/]+$||g;
my $PATH_COPRA = $ABS_PATH; 

my $help = "";
my $sRNAs_fasta = "input_sRNA.fa";
my $upstream = 200;
my $downstream = 100;
my $region = "5utr";
my $core_count = 1; # how many parallel processes are allowed
my $verbose = 0;
my $noclean = 0;
my $websrv = 0;
my $topcount = 100; # amount of top predictions //
my $root = 1; # root function to apply to the weights //
my $enrich = 0; ## functional enrichment needs to be specifically turned on 
                ## this option also allows to specify how many top predictions to use for the enrichment
my $genomePath = "."; # where to look for and store genome files
my $intarnaParamFile = $PATH_COPRA . "coprarna_aux/intarna_options.cfg";
my $CopraRNA_expert_options = $PATH_COPRA . "coprarna_aux/coprarna_options.cfg";
my $hybrid_threshold = 0.6; # interactions are removed from the CopraRNA calculations if a continuous hybrid covers >= "hybrid_threshold" of the sRNA
my $jalview_cores = 5; # number of cores used for jalview during post-processing. 0 indicates as many cores as defined in core_count



GetOptions (
    'help|?'			=> \$help,
    'srnaseq:s'			=> \$sRNAs_fasta,
    'ntup:i'			=> \$upstream,
    'ntdown:i'			=> \$downstream,
    'cores:i'			=> \$core_count,
	'jalview_cores:i'	=> \$jalview_cores,
    'region:s'			=> \$region, # one of "5utr", "3utr", "cds"
    'verbose'			=> \$verbose, # switch for verbose output during computation
    'websrv'			=> \$websrv, # switch for providing webserver output
    'noclean'		=> \$noclean, # switch to prevent cleaning of files
    'topcount:i'	=> \$topcount, # amount of top predictions to return
    'enrich:i'		=> \$enrich, # functional enrichment needs to be specifically turned on // also how many top preds to use for enrichment 
    'root:i'		=> \$root, # root function to apply to the weights
    'genomePath:s'		=> \$genomePath,
    'intarnaOptions:s'		=> \$intarnaParamFile,
    'CopraRNA_expert_options:s'		=> \$CopraRNA_expert_options,
	'hybrid_threshold:f'		=> \$hybrid_threshold
);

if ($help) {

print "\nCopraRNA ".$COPRARNA_VERSION."\n\n",

"CopraRNA is a tool for sRNA target prediction. It computes whole genome target predictions\n",
"by combination of distinct whole genome IntaRNA predictions. As input CopraRNA requires\n",
"at least 3 homologous sRNA sequences from 3 distinct organisms in FASTA format.\n", 
"Furthermore, each organisms' genome has to be part of the NCBI Reference Sequence (RefSeq)\n",
"database (i.e. it should have exactly this NZ_* or this NC_XXXXXX format where * stands\n",
"for any character and X stands for a digit between 0 and 9). Depending on sequence length\n",
"(target and sRNA), amount of input organisms and genome sizes, CopraRNA can take up to 24h\n",
"or longer to compute. In most cases it is significantly faster. It is suggested to run CopraRNA\n", 
"on a machine with at least 8 GB of memory.\n\n",
 
"CopraRNA produces a lot of file I/O. It is suggested to run CopraRNA in a dedicated\n",
"empty directory to avoid unexpected behavior.\n\n",

"The central result table is CopraRNA_result.csv. Further explanations concerning the files\n",
"in the run directory can be found in README.txt.\n\n",

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
" --jalview_cores           number of cores used for jalview during post-processing. 0 indicates as many cores as defined in cores (def:5)\n",
" --verbose                 switch to print verbose output to terminal during computation (def:off)\n",
" --websrv                  switch to provide webserver output files (def:off)\n",
" --noclean                 switch to prevent removal of temporary files (def:off)\n",
" --enrich                  if entered then DAVID-WS functional enrichment is calculated with given amount of top predictions (def:off)\n",
" --root                    specifies root function to apply to the weights (def:1)\n",
" --topcount                specifies the amount of top predictions to return and use for the extended regions plots (def:200)\n",
" --genomePath              path where NCBI genome files (*.gb) are to be stored (def:"." == working directory)\n\n",
" --intarnaOptions			path for IntaRNA parameter file\n\n",
" --CopraRNA_expert_options	path to parameter file for CopraRNA expert or experimental options\n\n",
" --hybrid_threshold		interactions are removed from the CopraRNA calculations if the hybrid covers >= hybrid_threshold of the sRNA\n\n",
"\n",
"Example call: CopraRNA2.pl -srnaseq sRNAs.fa -ntup 200 -ntdown 100 -region 5utr -enrich 200 -topcount 200 -cores 4\n\n",
"License: MIT\n\n",
"References: \n",
"1. Wright PR et al., Comparative genomics boosts target prediction for bacterial small RNAs\n   Proc Natl Acad Sci USA, 2013, 110 (37), E3487–E3496\n",
"2. Wright PR et al., CopraRNA and IntaRNA: predicting small RNA targets, networks and interaction domains\n   Nucleic Acids Research, 2014, 42 (W1), W119-W123\n",
"\n";

exit(0);

}

# input check
unless (-e $sRNAs_fasta) {
    die("\nError: No input FASTA supplied!\nUse '-h' option for help.\n\n");
}

# create genome path if necessary
(system("mkdir -p $genomePath") == 0) or die("\nError: could not create genome path '$genomePath'.\n\n");

# rudimentary check for fasta
my $check_fa = `grep '>' $sRNAs_fasta`;
chomp $check_fa;
die("\nError: The input file ($sRNAs_fasta) supplied does not appear to be a FASTA file!\n\n") unless($check_fa);

# check for sequence count in input being more than 2
my $count_fa = `grep -c '>' $sRNAs_fasta`;
chomp $count_fa;
die("\nError: The input file ($sRNAs_fasta) seems to contain less than 3 sequences!\n\n") unless($count_fa>2);

# create warning for non empty run dir
my @dir_files = <*>;
my $file_count = scalar(@dir_files);
print "\nWarning: your run directory contains many files ($file_count). In general it is suggested to run CopraRNA in an empty directory to avoid unexpected behaviour.\n\n" if ($file_count > 10);

# check core count // if wrong set it to 1 //
if ($core_count <= 1) {
    $core_count = 1;
}

# check region parameter
die ("\nError: -region parameter must be one of 5utr, 3utr or cds. You set '$region'.\n\n") unless ($region eq "5utr" or $region eq "3utr" or $region eq "cds");

# check for gaps
system "grep '-' $sRNAs_fasta > find_gaps.txt";
if (-s "find_gaps.txt") { die("\nError: Gaps are present in sRNA sequences. Please delete them from the file and restart.\n\n"); }

# check for correct RefSeq formatted headers and their presence in the availibility table
my $headerIDs = `grep ">" $sRNAs_fasta | sed 's/>//g' | tr '\n' ';'`;
chop $headerIDs;
my @splitHeaderIDs = split(/;/,$headerIDs);
foreach(@splitHeaderIDs) {
    die("\nError: $_ does not match correct RefSeq ID format (NZ_* or NC_XXXXXX where * stands for any character and X stands for a digit between 0 and 9).\n\n") unless ($_ =~ m/NC_\d{6}|NZ_.*/);
    $_ =~ s/^\s+|\s+$//g; 
    my $availabilityCheck = `grep '$_' $PATH_COPRA/coprarna_aux/kegg2refseqnew.csv`;
    die("\nError: '$_' is not present in the availability list and is thus not compatible with CopraRNA.\n\n") unless (length $availabilityCheck); 
}

# check that maxbpdist ist smaller or equal to windowsize

#die("\nError: The maximal basepair distance ($maxbpdist) is larger than the given window size ($winsize) but must be <= to the windows size. Please change the parameters accordingly.\n\n") if ($maxbpdist > $winsize);

# check for ntup + ntdown being >= $winsize because of IntaRNA parameters
# my $sumUpDown = $upstream+$downstream;
# if ($region eq "5utr" or $region eq "3utr") {
    # die("\nError: (-ntup + -ntdown) is $sumUpDown but must be >= $winsize (--winsize). Please change the parameters accordingly.\n\n") if ( $sumUpDown < $winsize );
# }

# check for duplicate IDs in FASTA header // not allowed
my $duplicate_fasta_header = `grep ">" $sRNAs_fasta | sort | uniq -d`;
chomp $duplicate_fasta_header;
die("\nError: Duplicate organisms ($duplicate_fasta_header) are present in $sRNAs_fasta.\n\n") if($duplicate_fasta_header);

# write input options to file
open WRITETOOPTIONS, ">", "CopraRNA_option_file.txt";
    print WRITETOOPTIONS "sRNA FASTA:" . $sRNAs_fasta . "\n";
    print WRITETOOPTIONS "nt upstream:" . $upstream . "\n";
    print WRITETOOPTIONS "nt downstream:" . $downstream . "\n";
    print WRITETOOPTIONS "region:" . $region . "\n";
    print WRITETOOPTIONS "core count:" . $core_count . "\n";
	print WRITETOOPTIONS "jalview_cores:" . $jalview_cores . "\n";
    print WRITETOOPTIONS "verbose:" . $verbose . "\n";
    print WRITETOOPTIONS "websrv:" . $websrv . "\n";
    print WRITETOOPTIONS "top count:" . $topcount . "\n";
    print WRITETOOPTIONS "root:" . $root . "\n";
    print WRITETOOPTIONS "enrich:" . $enrich . "\n";
    print WRITETOOPTIONS "noclean:" . $noclean . "\n";
    print WRITETOOPTIONS "version:CopraRNA ".$COPRARNA_VERSION."\n";
    print WRITETOOPTIONS "genomePath:$genomePath\n";
    print WRITETOOPTIONS "intarnaOptions:$intarnaParamFile\n";
	print WRITETOOPTIONS "CopraRNA_expert_options:$CopraRNA_expert_options\n";	
	print WRITETOOPTIONS "hybrid_threshold:$hybrid_threshold\n";
	
	
	
close WRITETOOPTIONS;
# end write options

# CopraRNA error log
system "touch err.log";

system "cp $sRNAs_fasta input_sRNA.fa" unless ($sRNAs_fasta eq "input_sRNA.fa");
$sRNAs_fasta = "input_sRNA.fa";

# format sRNA fasta - put sequence in 1 line
system $PATH_COPRA . "coprarna_aux/format_fasta.pl $sRNAs_fasta" . " > $sRNAs_fasta.temp";
system "mv $sRNAs_fasta.temp $sRNAs_fasta";

# defining intarna options
if($intarnaParamFile ne 'NA') {
    $intarnaParamFile = `readlink -f $intarnaParamFile`;
} else {
    $intarnaParamFile = $PATH_COPRA . "coprarna_aux/intarna_options.cfg";
}
chomp $intarnaParamFile;

# build RefSeq input based on the sRNA input fasta (can only contain refseq IDs in header)
my $RefSeqIds = `grep '>' $sRNAs_fasta | sed 's/>//g' | tr '\n' ' '`;

# print org of interest
if ($verbose) {
    my $ooi_rid = `grep ">" input_sRNA.fa | head -n1 | sed 's/>//g'`;
    chomp $ooi_rid;
    my $full_ooi = `grep '$ooi_rid' "$PATH_COPRA/coprarna_aux/CopraRNA_available_organisms.txt"`;
    chomp $full_ooi;
    print "\nOrganism of interest: $full_ooi\n\n";
}

# my $homology_intaRNA_call=$PATH_COPRA . "coprarna_aux/homology_intaRNA.pl $sRNAs_fasta $upstream $downstream $region $RefSeqIds";
my $homology_intaRNA_call=$PATH_COPRA . "coprarna_aux/homology_intaRNA.pl $sRNAs_fasta $upstream $downstream $region $core_count $intarnaParamFile $RefSeqIds";
print $homology_intaRNA_call . "\n" if ($verbose);

my $homology_intaRNA_exitStatus = system $homology_intaRNA_call;
$homology_intaRNA_exitStatus /= 256; # get original exit value
# check exit status
if ($homology_intaRNA_exitStatus != 0) { 
	die ("\nERROR: homology_intaRNA.pl returned with exit code $homology_intaRNA_exitStatus. Something went wrong, so please check the error files!\n\n");
}

# get organism of interest
my $ncrnaRIDs = `grep ">" ncrna.fa | sed 's/>ncRNA_//g' | tr '\n' ' '`;
my @splitRID = split(/\s/, $ncrnaRIDs);
my $organismOfInterest = $splitRID[0];
chomp $organismOfInterest;
# final optimal IntaRNA result for organism of interest
my $MainFinalCSV = $organismOfInterest . "_upfromstartpos_" . $upstream . "_down_" . $downstream . "_opt.intarna.csv";

if ($enrich) {
	############################################################
    print "Performing auxiliary enrichment\n" if ($verbose);
	############################################################
    # add IntaRNA single organisms chart reports for aux enrichment // sort by p-value
    system "env LC_ALL=C sort -t';' -g -k36 $MainFinalCSV -o intarna_websrv_table.csv";
    system $PATH_COPRA . "coprarna_aux/add_GI_genename_annotation_intarna.pl";
    system $PATH_COPRA . "coprarna_aux/DAVIDWebService_IntaRNA_chartReport.py intarna_websrv_table_ncbi.csv $enrich > IntaRNA_chartReport.txt"; ## aux einrich for same amout as regular enrichment
    system "grep -P 'geneIds\\s=|termName\\s=' IntaRNA_chartReport.txt | sed 's/^[ ]*//g' | sed 's/ = /=/g' | sed 's/, /,/g' | sed 's/\"//g' > IntaRNA_chartReport_grepped.txt";
    system $PATH_COPRA . "coprarna_aux/find_single_specific_targets_in_termCluster.pl > org_of_interest_aux_enrichment.txt";
    system "echo 'locus_tag,start_tar,stop_tar,start_query,stop_query,energy,p-value,gene_name,gene_id,annotation,functional_terms' > aux_table.csv";
    system "awk -F';' '{ print \$1\",\"\$9\",\"\$10\",\"\$11\",\"\$12\",\"\$15\",\"\$36\",\"\$37\",\"\$38\",\"\$39 }' intarna_websrv_table_ncbi.csv > intarna_websrv_table_ncbi_awk.csv";
    my $aux_gids = `grep -oP ';\\d+\\(' org_of_interest_aux_enrichment.txt | sed 's/[;(]//g' | sort -u | tr '\n' ';'`;
    chop $aux_gids; # remove trailing ';' 
    my @split_aux_gids = split(/;/, $aux_gids);
    foreach(@split_aux_gids) {
        system "grep -P ',$_,' intarna_websrv_table_ncbi_awk.csv | tr '\n' ',' >> aux_table.csv";
        system "grep -P '$_' org_of_interest_aux_enrichment.txt | awk -F';' '{ print \$1 }' | tr '\n' ';' | sed 's/,/ /' | sed 's/.\$//' >> aux_table.csv";
        system "echo >> aux_table.csv";
    }
}

# output warnings
#system "awk -F ';' '{if (\$2 > 0.5) { print toupper(\$1) \" may be overweighted. It has weight\"; printf(\"\%.2f\", \$2); print \". You should consider checking the 16S rDNA tree. We suggest removal of outliers from your input and restarting.\";} }' zscore.weight | tr '\n' ' ' > weights.warning";

# check err.log // should be empty // err.log contains information on 
# 1. not correctly downloaded RefSeq files 
# 2. gene no CDS issue 
# 3. wrong 16S counts
# 4. empty CopraRNA_result.csv 
# 5. the exception in add_pval_to_csv_evdfit.R
# -s  File has nonzero size (returns size in bytes).
if (-s "err.log") { die("\nError: CopraRNA failed. Check err.log for details.\n\n"); } ## added another check here at the bottom maybe we need some more check hooks in the new scripts

# create regions plots
print "Preparing interaction plots\n" if ($verbose);
system "R --slave -f " . $PATH_COPRA . "coprarna_aux/script_R_plots_8.R --args CopraRNA_result_all.csv $topcount 2> /dev/null > /dev/null"; ## changed input file and piping command line output to /dev/null for silencing // 


# thumbnails png
if ($websrv) {
    system "convert -size 170x170 -resize 170x170 sRNA_regions_with_histogram.png thumbnail_sRNA.png";
    system "convert -size 170x170 -resize 170x170 mRNA_regions_with_histogram.png thumbnail_mRNA.png";
}



#######################################################
print "prepare html output\n" if ($verbose);
#######################################################
print "CopraRNA2html.r\n" if ($verbose);
system "R --slave -f " . $PATH_COPRA . "coprarna_aux/CopraRNA2html.r 2>> CopraRNA2_subprocess.oe 1>&2";



# clean up
unless ($noclean) {

    print "Cleaning run directory\n" if ($verbose);
	system "bash $PATH_COPRA/coprarna_aux/copraRNA_cleanupWD.sh";

}
