#!/usr/bin/env perl

use strict;
use warnings;

use Parallel::ForkManager;
use List::MoreUtils qw(uniq);
use Cwd 'abs_path'; 

my $ncrnas = $ARGV[0];
my $upfromstartpos = $ARGV[1]; 
my $down = $ARGV[2]; 
my $mrnapart = $ARGV[3];
my $core_count = $ARGV[4]; 
my $intarnaParamFile = $ARGV[5];
my $refseqid = '';

# files dedicated to capture output of subcalls for debugging
my $OUT_ERR = "CopraRNA2_subprocess.oe";

my $orgcnt = (scalar(@ARGV) - 6);

# get absolute path
my $ABS_PATH = abs_path($0); 
# remove script name at the end
# match all non slash characters at the end of the string
$ABS_PATH =~ s|[^/]+$||g; 
my $PATH_COPRA_SUBSCRIPTS = $ABS_PATH;

# get core count from option file
my $cores = `grep 'core count:' CopraRNA_option_file.txt | grep -oP '\\d+'`; 
chomp $cores;

# check for verbose printing
my $verbose = `grep 'verbose:' CopraRNA_option_file.txt | sed 's/verbose://g'`; 
chomp $verbose;


my $pm = new Parallel::ForkManager($cores);

for(my $i=6;$i<=scalar(@ARGV) - 1;$i++) {

    my @splitarg = split(/,/, $ARGV[$i]);

    if ($splitarg[0] =~ m/(N[ZC]_.+?)\.gb(\.gz)?/) { 
        $refseqid = $1;
    }
    $pm->start and next; 
    my $outfile = $refseqid . '_upfromstartpos_' . $upfromstartpos . '_down_' . $down . '.fa';

    my $splitargcount = 1;
    foreach (@splitarg) {
        my $tempOutfile = $outfile; 
        $tempOutfile = $tempOutfile . $splitargcount; 
        system $PATH_COPRA_SUBSCRIPTS . "parse_region_from_genome.pl $_ $upfromstartpos $down $mrnapart > $tempOutfile"; 
        $splitargcount++;   
    }
    $pm->finish;
}
$pm->wait_all_children;

my @files = ();
@files = <*>;
my @rfids = ();

foreach (@files) {
    if ($_ =~ m/(N[ZC]_.+?)\.gb(\.gz)?/) { 
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
open(MYDATA, $ncrnas) or die("\nError: cannot open file $ncrnas in prepare_intarna_out.pl\n\n"); 
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

# my $pm = new Parallel::ForkManager($cores);   
foreach (@files) {
    if ($_ =~ m/(N[ZC]_.+)$suffix$/) { 
        my $refid = $1;  # get refseq id
        foreach my $line (@ncrnalines) {
            if ($switch) {                  
                push(@ncrnaarray, $line);   
                $switch = 0;                
            }                               
            if($line =~ m/$refid/) {
                $switch = 1;
                push(@ncrnaarray, $line);
            }
        }
        my $ncrnafilename = $ncrnaarray[0]; 
		my $ncrnafilename1 = $ncrnaarray[0];
			$ncrnafilename = $ncrnafilename . ".fa";
			@ncrnaarray = ();
	        my $intarnaout = $_ . ".intarna.csv"; 
			my $intarnasortedout = $_ . ".intarna.sorted.csv"; 
			# create temporary "sorted" file to support old pipeline
			my $intarna_call = 
					"IntaRNA"
					." --outOverlap=Q"
					." --target $_"
					." --query $ncrnafilename"
					." --outNumber=2"
                    ." --parameterFile $intarnaParamFile"
                    ." --threads $cores"
					." --outMode=C --outCsvCols 'id1,id2,seq1,seq2,subseq1,subseq2,subseqDP,subseqDB,start1,end1,start2,end2,hybridDP,hybridDB,E,ED1,ED2,Pu1,Pu2,E_init,E_loops,E_dangleL,E_dangleR,E_endL,E_endR,seedStart1,seedEnd1,seedStart2,seedEnd2,seedE,seedED1,seedED2,seedPu1,seedPu2,E_norm'"
					." --out $intarnaout"
					." --outCsvSort E"
					."; "
					."ln -s $intarnaout $intarnasortedout";
			print($intarna_call . "\n") if ($verbose);
            system($intarna_call) unless (-e $intarnaout);
     }                                                                                   
}
                                                                                        


my %ncrnalengthhash = (); 
my @lines = (); 
my @datalines = (); 

## disentangle *.fa.intarna.sorted.csv

# this needs to be here so we can calculate
# the pvalues on the distribution returned
# from the optimal results and not from the
# mixture of optimal and suboptimals

## creates *_opt.intarna.csv and *_subopt.intarna.csv files
system $PATH_COPRA_SUBSCRIPTS . "disentangle_sorted_intarna_CSV.pl"; 

## adds pvalues to *_opt.intarna.csv and *_subopt.intarna.csv files // file names stay the same
system "R --slave -f " . $PATH_COPRA_SUBSCRIPTS . "add_pval_to_csv_evdfit.R"; 
                                                                              ## makes add_pval_to_csv_evdfit.pl obsolete
																			  
																			  
# re-cluster based on 5'UTRs
my $refineClusterCall = $PATH_COPRA_SUBSCRIPTS . "refine_clustertab.r"; 
print "$refineClusterCall\n"  if ($verbose);
system "Rscript --slave $refineClusterCall"; 																			  


# remove interactions with full hybrids
print "remove_full_hybrids\n"  if ($verbose);
system "R --slave -f " . $PATH_COPRA_SUBSCRIPTS . "remove_full_hybrids.r>> $OUT_ERR 1>&2"; 
 
																			  
																			  
## create opt_tags.clustered
system $PATH_COPRA_SUBSCRIPTS . "cluster_intarna_csv.pl > opt_tags.clustered"; 


system "mafft --localpair --quiet 16s_sequences.fa > 16s_sequences.aln";
#system "distmat -sequence 16s_sequences.aln -nucmethod 1 -outfile distmat.out 2>> $OUT_ERR 1>&2"; 
#system $PATH_COPRA_SUBSCRIPTS . "transform_distmat.pl distmat.out > compatible.distmat";

