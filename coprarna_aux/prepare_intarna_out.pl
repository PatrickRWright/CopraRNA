#!/usr/bin/env perl

use strict;
use warnings;

use Parallel::ForkManager;
use List::MoreUtils qw(uniq);
use Cwd 'abs_path'; ## edit 2.0.5.1

my $ncrnas = $ARGV[0];
my $upfromstartpos = $ARGV[1]; 
my $down = $ARGV[2]; 
my $mrnapart = $ARGV[3]; 
my $refseqid = '';

my $orgcnt = (scalar(@ARGV) - 4);

# get absolute path
my $ABS_PATH = abs_path($0); ## edit 2.0.5.1
# remove script name at the end
# match all non slash characters at the end of the string
$ABS_PATH =~ s|[^/]+$||g; ## edit 2.0.5.1
my $PATH_COPRA_SUBSCRIPTS = $ABS_PATH;

# get core count from option file
my $cores = `grep 'core count:' CopraRNA_option_file.txt | grep -oP '\\d+'`; ## edit 2.0.4
chomp $cores; ## edit 2.0.4

# check if CopraRNA2 prediction should be made
my $cop2 = `grep 'CopraRNA2:' CopraRNA_option_file.txt | sed 's/CopraRNA2://g'`; ## edit 2.0.5.1
chomp $cop2;

# check for verbose printing
my $verbose = `grep 'verbose:' CopraRNA_option_file.txt | sed 's/verbose://g'`; ## edit 2.0.5.1
chomp $verbose;

# get window size option
my $winsize = `grep 'win size:' CopraRNA_option_file.txt | sed 's/win size://g'`; ## edit 2.0.5.1
chomp $winsize;

# get maximum base pair distance 
my $maxbpdist = `grep 'max bp dist:' CopraRNA_option_file.txt | sed 's/max bp dist://g'`; ## edit 2.0.5.1
chomp $maxbpdist;

for(my $i=4;$i<=scalar(@ARGV) - 1;$i++) {

    my @splitarg = split(/,/, $ARGV[$i]);

    if ($splitarg[0] =~ m/(N[ZC]_.+?)\.gb/) { ## edit 2.0.2
        $refseqid = $1;
    }

    my $outfile = $refseqid . '_upfromstartpos_' . $upfromstartpos . '_down_' . $down . '.fa';

    my $splitargcount = 1;

    foreach (@splitarg) {
        my $tempOutfile = $outfile; ## edit 1.2.2
        $tempOutfile = $tempOutfile . $splitargcount; ## edit 1.2.2
        system $PATH_COPRA_SUBSCRIPTS . "parse_region_from_genome.pl $_ $upfromstartpos $down $mrnapart > $tempOutfile"; ## edit 1.2.2
        $splitargcount++;   
    }
}

my @files = ();
@files = <*>;
my @rfids = ();

foreach (@files) {
    if ($_ =~ m/(N[ZC]_.+?)\.gb/) { # edit 2.0.2
        push (@rfids, $1)
    }
}

@rfids = uniq(@rfids);

foreach my $id (@rfids) {
    foreach my $file (@files) {
        if ($file =~ m/($id\S+\.fa)\d+/) {
             my $tempfile = $1;
             system "cat $file >> $tempfile";
        }
    }
}

## pairwise whole genome IntaRNA 
@files = <*>;
my $switch = 0;
my @ncrnaarray = (); 

# preparing individual sRNA files for IntaRNA whole genome predictions
open(MYDATA, $ncrnas) or die("\nError: cannot open file $ncrnas in prepare_intarna_out.pl\n\n"); ## edit 2.0.5.1 // added script name
    my @ncrnalines = <MYDATA>;
close MYDATA;

foreach (@ncrnalines) {
    if ($_ =~ m/>/) {
        $_ =~ s/\r|\n|\s|\t//g;
        $_ = reverse $_;
        chop $_;
        $_ = reverse $_;
        $_ = $_ . ".fa";
        open (FILE, ">$_");
        chop $_;
        chop $_;
        chop $_;
        print FILE (">" . $_ . "\n");
    } else {
        $_ =~ s/\r|\n|\s|\t|-//g;
        $_ = lc($_);
        print FILE ($_);
    }
}
close FILE;

my $suffix = '_upfromstartpos_' . $upfromstartpos . '_down_' . $down . '.fa';

my $pm = new Parallel::ForkManager($cores);   
foreach (@files) {
    if ($_ =~ m/(N[ZC]_.+)$suffix$/) { ## edit 2.0.2
        my $refid = $1;  # get refseq id ## edit 2.0.2
        foreach my $line (@ncrnalines) {
            if ($switch) {                  
                push(@ncrnaarray, $line);   
                $switch = 0;                
            }                               
            if($line =~ m/$refid/) { ## edit 2.0.2
                $switch = 1;
                push(@ncrnaarray, $line);
            }
        }
        my $ncrnafilename = $ncrnaarray[0]; 
        $ncrnafilename = $ncrnafilename . ".fa";
        @ncrnaarray = ();
        my $intarnaout = $_ . ".intarna.csv"; ## edit 2.0.5.1 // added .csv
            $pm->start and next;            
            system("IntaRNA --tAccW $winsize --tAccL $maxbpdist --outNumber 2 --target $_ --query $ncrnafilename --outCsvCols 'id1,id2,seq1,seq2,subseq1,subseq2,subseqDP,subseqDB,start1,end1,start2,end2,hybridDP,hybridDB,E,ED1,ED2,Pu1,Pu2,E_init,E_loops,E_dangleL,E_dangleR,E_endL,E_endR,seedStart1,seedEnd1,seedStart2,seedEnd2,seedE,seedED1,seedED2,seedPu1,seedPu2,E_norm' --outMode=C --out $intarnaout") unless (-e $intarnaout); ## edit 2.0.4 added -s 1 for consensus prediction // ## edit 2.0.5 changed to IntaRNA2 // ## edit 2.0.5.1 changed to IntaRNA 2.0 csv output
            print("IntaRNA --tAccW $winsize --tAccL $maxbpdist --outNumber 2 --target $_ --query $ncrnafilename --outCsvCols 'id1,id2,seq1,seq2,subseq1,subseq2,subseqDP,subseqDB,start1,end1,start2,end2,hybridDP,hybridDB,E,ED1,ED2,Pu1,Pu2,E_init,E_loops,E_dangleL,E_dangleR,E_endL,E_endR,seedStart1,seedEnd1,seedStart2,seedEnd2,seedE,seedED1,seedED2,seedPu1,seedPu2,E_norm' --outMode=C --out $intarnaout\n") if ($verbose);
            $pm->finish;
     }                                                                                   
}
                                                                                        
$pm->wait_all_children;  

## sort IntaRNA output by energy
@files = <*intarna.csv>; ## edit 2.0.4 // added *csv // ## edit 2.0.5.1 // changed to *intarna.csv

foreach (@files) {
    ## edit 2.0.4 // removed pattern match for files
    my $temp = $_;
    chomp $temp;
    chop $temp;
    chop $temp;
    chop $temp;
    my $sortedcsv = $temp . "sorted.csv"; 
    system "env LC_ALL=C sort -t';' -g -k15 $_ -o $sortedcsv"; ## edit 2.0.5.1 // removed sort_intarna_csv_results.pl and switched to linux sort
}

@files = ();
@files = <*>;

my %ncrnalengthhash = (); 
my @lines = (); 
my @datalines = (); 

## disentangle *.fa.intarna.sorted.csv

# this needs to be here so we can calculate
# the pvalues on the distribution returned
# from the optimal results and not from the
# mixture of optimal and suboptimals

## creates *_opt.intarna.csv and *_subopt.intarna.csv files
system $PATH_COPRA_SUBSCRIPTS . "disentangle_sorted_intarna_CSV.pl"; ## edit 2.0.4

## adds pvalues to *_opt.intarna.csv and *_subopt.intarna.csv files // file names stay the same
system "R --slave -f " . $PATH_COPRA_SUBSCRIPTS . "add_pval_to_csv_evdfit.R"; ## edit 2.0.4 // changed this to an R script
                                                                              ## makes add_pval_to_csv_evdfit.pl obsolete
## create opt_tags.clustered
system $PATH_COPRA_SUBSCRIPTS . "cluster_intarna_csv.pl > opt_tags.clustered"; ## edit 2.0.4.1 // reimplemented cluster_intarna_csv.pl from hash_clusters

## create opt_tags.clustered_rcsize
system "R --slave -f " . $PATH_COPRA_SUBSCRIPTS . "remove_clusters_under_percantage.R"; ## edit 2.0.1

## need to do this now -> compatible.distmat needed in prev. weight calc 
system "mafft --localpair --quiet 16s_sequences.fa > 16s_sequences.aln";
system "distmat -sequence 16s_sequences.aln -nucmethod 1 -outfile distmat.out 2> /dev/null"; ## edit 2.0.5.1 // added 2> /dev/null to prevent output to the terminal
system $PATH_COPRA_SUBSCRIPTS . "transform_distmat.pl distmat.out > compatible.distmat";

## prepare opt_tags.clustered_trunc for cop2 prediction // ## edit 2.0.5.1

if ($cop2) {
    ## get p-value cutoff
    my $pvcut = `grep 'p-value cutoff:' CopraRNA_option_file.txt | grep -oP '\\d+\\.\\d+'`; ## edit 2.0.5.1
    chomp $pvcut; ## edit 2.0.5.1

    ## truncate the tags.clustered file according to the specified pvalue cutoff
    ## 36 is the field with the p-value and 0 is the line
    system "awk -F';' '{ if (\$36 <= $pvcut) print \$0 }' opt_tags.clustered > opt_tags.clustered_trunc_no_head";
    system "head -n 1 opt_tags.clustered > temp_head";
    system "cat temp_head opt_tags.clustered_trunc_no_head > opt_tags.clustered_trunc";
    system "rm opt_tags.clustered_trunc_no_head temp_head";

    # cop1: opt_tags.clustered_rcsize
    # cop2: opt_tags.clustered_trunc

    # calculate all weights for coprarna2 prediction
    system $PATH_COPRA_SUBSCRIPTS . "parallelize_mafft_for_weights.pl opt_tags.clustered_trunc";
}


