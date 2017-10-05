#!/usr/bin/env perl

use strict;
use warnings;

use Parallel::ForkManager;

# computes mafft alignments for putative targets in parallel

my $input_table = $ARGV[0]; # CopraRNA2_prep_anno_addhomologs_padj_amountsamp.csv

# get core count from option file
my $cores = `grep 'core count:' CopraRNA_option_file.txt | grep -oP '\\d+'`;
chomp $cores;
 
open(MYDATA, $input_table) or die("\nError: cannot open file $input_table at parallelize_target_alignments.pl\n\n");
    my @input_table_lines = ();
    @input_table_lines = <MYDATA>;
close MYDATA;

# read sequences into hash to speed up searching
my %ltag_to_seq_hash = (); # ltag -> sequence
my @target_fastas = <*upfromstartpos*down*fa>;

foreach (@target_fastas) {

    open(FASTA, $_) or die("\nError: cannot open file $_ at parallelize_target_alignments.pl\n\n");
        my @fasta_lines = <FASTA>;
    close FASTA;

    for (my $l=0;$l<scalar(@fasta_lines);$l=$l+2) {
        if ($fasta_lines[$l] =~ m/>(.+)/) {
            $ltag_to_seq_hash{lc($1)} = $fasta_lines[($l+1)];
        } else {
            print "\nError unexpected behaviour!\n\n";
        }
    }
}

system "mkdir target_alignments";

# build fasta files
for (my $i=1;$i<scalar(@input_table_lines);$i++) {

    my $curr_line = $input_table_lines[$i];
    chomp $curr_line;

    system "touch target_alignments/$i.fa"; 
 
    open(WRITEFASTA, ">>target_alignments/$i.fa");    
 
    my @split = split(/,/, $curr_line);
    pop @split; # remove amount sampled
    pop @split; # remove add homologs
    pop @split; # remove annotation
    shift @split; # remove fdr
    shift @split; # remove p-value
    foreach (@split) {
        if ($_ =~ m/(.+)\(.+/) { # get locus tag
            my $ltag = lc($1);
            print WRITEFASTA ">$ltag\n$ltag_to_seq_hash{lc($ltag)}";
        }
    }
    close WRITEFASTA; 
}

# make alignments
my @target_fasta_files = <target_alignments/*fa>;

my $pm = new Parallel::ForkManager($cores);

foreach(@target_fasta_files) {
    print $_ . "\n";
    my @split = split(/\./, $_);
    my $mafft_out = $split[0] . ".aln";

    $pm->start and next;
    system "mafft --localpair --quiet $_ > $mafft_out";
    $pm->finish;
}

$pm->wait_all_children;

