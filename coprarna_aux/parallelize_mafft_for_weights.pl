#!/usr/bin/env perl

use warnings;
use strict;

# ## edit 2.0.4 // wrote this script for 2.0.4
# run with $tags_clustered as argument in CopraRNA directory
# creates a directory "weight_permutations" and aligns
# 16s sequences for all combinations present in $tags_clustered
# without recomputing already used combinations

# this script needs to run after the compatible.distmat for all
# organisms has been created (because of the long ID issue)

# tested this with an individual example of 26 sequences out
# of a set of 32 and the results in the weight file were 
# exactly the same

use Parallel::ForkManager;
use Cwd 'abs_path'; ## edit 2.0.5.1

# get absolute path
my $ABS_PATH = abs_path($0); ## edit 2.0.5.1
# remove script name at the end
# match all non slash characters at the end of the string
$ABS_PATH =~ s|[^/]+$||g; ## edit 2.0.5.1
my $PATH_COPRA_SUBSCRIPTS = $ABS_PATH;

my $tags_clustered = $ARGV[0];

my $cores = `grep 'core count:' CopraRNA_option_file.txt | grep -oP '\\d+'`; ## edit 2.0.4
chomp $cores; ## edit 2.0.4

system "mkdir weight_permutations";

my %org_combination_hash = (); # RefSeqID_RefSeqID_RefSeqID -> 1

open (MYDATA, $tags_clustered) or die ("\nError: cannot open file $tags_clustered at parallelize_mafft_for_weights.pl\n\n"); ## edit 2.0.5.1
    my @tags_clustered = <MYDATA>;
close (MYDATA);

my $curr_cluster_num = 0;
my @org_array = ();

for (my $i=1;$i<scalar(@tags_clustered);$i++) { 

    my $curr_line = $tags_clustered[$i];
    chomp $curr_line;
    my @split_curr_line = split(/;/, $curr_line);
    $curr_cluster_num = $split_curr_line[-1];
    my $curr_org = uc($split_curr_line[1]); ## edit 2.0.5.1 // org is now in field 1
    $curr_org =~ s/NCRNA_//g;
    push(@org_array, $curr_org);

    my $next_line = "";
    my $next_line_cl_num = "";
    
    if(exists($tags_clustered[($i+1)])) {
        $next_line = $tags_clustered[($i+1)];
        chomp $next_line;
        my @split_next_line = split(/;/, $next_line);
        $next_line_cl_num = $split_next_line[-1];
    } else { # exception for the last line in the file -> create hash entry
        @org_array = sort(@org_array);
        my $org_hash_key = join("-",@org_array);
        if(exists $org_combination_hash{$org_hash_key}) {
            # nothing
        } else {
            $org_combination_hash{$org_hash_key} = 1 unless (scalar(@org_array)<3); # init
        }
        last;
    } 
    
    if($curr_cluster_num ne $next_line_cl_num) { # next line is different cluster -> create hash entry
        @org_array = sort(@org_array); 
        my $org_hash_key = join("-",@org_array);
        if(exists $org_combination_hash{$org_hash_key}) {
            # nothing
        } else {
            $org_combination_hash{$org_hash_key} = 1 unless (scalar(@org_array)<3); # init
        }
        @org_array = ();
    }
}

my $file_number = 0;

my $pm = new Parallel::ForkManager($cores);

foreach my $key (keys %org_combination_hash) {
    # file specs need to be done before parallel processes
    $file_number++;

    $pm->start and next;
    my @split_orgs = split("-", $key);
    my $file_name = $file_number . "_16s.fa";
    foreach(@split_orgs) {
        system "grep -A 1 '$_' 16s_sequences.fa >> $file_name"; # sequence is only in one line
    } 
    system "mv $file_name weight_permutations";
    $pm->finish;
}
$pm->wait_all_children;

my @sixteens_dirs = <weight_permutations/*fa>;

## run mafft in parallel
## takes about 60 mins on 4 cores with 32 orgs
## mafft is the most intensive part
foreach(@sixteens_dirs) {
    
    # file specs need to be done before parallel processes
    my @split = split(/\./, $_);
    my $mafft_out = $split[0] . ".aln";
    my $distmat_out = $split[0] . "_distmat.out";
    my $compatible_distmat_out = $split[0] . "_compatible.distmat";   
    
    $pm->start and next;
    system "mafft --localpair --quiet $_ > $mafft_out" unless (-e $mafft_out);
    system "distmat -sequence $mafft_out -nucmethod 1 -outfile $distmat_out 2> /dev/null" unless (-e $distmat_out); 
    system "$PATH_COPRA_SUBSCRIPTS/transform_distmat.pl $distmat_out > $compatible_distmat_out" unless (-e $compatible_distmat_out);
    $pm->finish;
}
$pm->wait_all_children;

# do the mapping of the ID issue
open (DISTMAT, "compatible.distmat") or die ("\nError: cannot open compatible.distmat in parallelize_mafft_for_weights.pl\n\n"); ## edit 2.0.5.1 // adjusted Error output
    my @distmat_lines = <DISTMAT>;
close (DISTMAT);

my $c = 1000; ## edit 2.0.5.1 // changed this to 1000 because it was causing  wrong matches when more than 10 long IDs were present
my %ID_to_ID_hash = (); # new_id -> old_id

for (my $i=1;$i<scalar(@distmat_lines);$i++) {
    my $curr_line = $distmat_lines[$i];
    chomp $curr_line;
    my @split = split(/\t+/,$curr_line);
    my $RID = $split[0];
    if( length($RID)>10 ) { # remap ID
        my $new_id = "na_" . $c;
        $ID_to_ID_hash{$new_id} = $RID;
        $c++;
    }
}

# map replacement IDs
foreach my $key (keys %ID_to_ID_hash) {

    my $new_id = $key;
    my $refseq_id = $ID_to_ID_hash{$key};

    system "sed -i 's/$refseq_id/$new_id/g' weight_permutations/*_compatible.distmat";
}

# get trees
my @compatible_distmat = <weight_permutations/*_compatible.distmat>;

foreach(@compatible_distmat) {

    my @split = split(/\./, $_);
    my $fneighbor_out = $split[0] . ".fneighbor";
    my $treefile_out = $split[0] . ".treefile";

    $pm->start and next;
    system "fneighbor -datafile $_ -outfile $fneighbor_out -outtreefile $treefile_out -treetype u > /dev/null 2> /dev/null";
    system "sed -i 's/0.00000/0.00001/g' $fneighbor_out";
    system "sed -i 's/0.00000/0.00001/g' $treefile_out";
    $pm->finish;
}
$pm->wait_all_children;

# map original IDs back
foreach my $key (keys %ID_to_ID_hash) {

    my $new_id = $key;
    my $refseq_id = $ID_to_ID_hash{$key};

    system "sed -i 's/$new_id/$refseq_id/g' weight_permutations/*_compatible.fneighbor";
    system "sed -i 's/$new_id/$refseq_id/g' weight_permutations/*_compatible.treefile";
}

foreach(@compatible_distmat) {
    my @split = split(/\./, $_);
    my $fneighbor_out = $split[0] . ".fneighbor";
    my $treefile_out = $split[0] . ".treefile";
    my $zscore_weight_out = $split[0] . "_zscore.weight";
    system "$PATH_COPRA_SUBSCRIPTS/read_fneighfile.pl $treefile_out $fneighbor_out > $zscore_weight_out";
}

