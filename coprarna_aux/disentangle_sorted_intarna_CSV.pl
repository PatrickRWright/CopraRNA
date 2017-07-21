#!/usr/bin/env perl

use strict;
use warnings;

# written for ## edit 2.0.4
# run this in a CopraRNA output dir without arguments
# the directory has to have run a CopraRNA job that returned
# one suboptimal interaction within the IntaRNA part 
# IntaRNA --outNumber 2

my @CSVs = <*.fa.intarna.sorted.csv>;

foreach (@CSVs) {

    my @split_file_name = split(/\./, $_);

    my $opt_print = $split_file_name[0] . "_opt.intarna.csv";
    my $subopt_print = $split_file_name[0] . "_subopt.intarna.csv";

    open WRITETOOPT, ">", $opt_print;
    open WRITETOSUBOPT, ">", $subopt_print;

    # headers   
    print WRITETOOPT "id1;id2;seq1;seq2;subseq1;subseq2;subseqDP;subseqDB;start1;end1;start2;end2;hybridDP;hybridDB;E;ED1;ED2;Pu1;Pu2;E_init;E_loops;E_dangleL;E_dangleR;E_endL;E_endR;seedStart1;seedEnd1;seedStart2;seedEnd2;seedE;seedED1;seedED2;seedPu1;seedPu2\n"; 
    print WRITETOSUBOPT "id1;id2;seq1;seq2;subseq1;subseq2;subseqDP;subseqDB;start1;end1;start2;end2;hybridDP;hybridDB;E;ED1;ED2;Pu1;Pu2;E_init;E_loops;E_dangleL;E_dangleR;E_endL;E_endR;seedStart1;seedEnd1;seedStart2;seedEnd2;seedE;seedED1;seedED2;seedPu1;seedPu2\n"; 

    my %opt_hash = ();
    my %subopt_hash = ();

    open MYDATA, $_ or die("\nError: Cannot open $_ !\n\n");
        my @lines = <MYDATA>;
    close MYDATA;


    for (my $i=1;$i<scalar(@lines);$i++) { # skip header

        my $curr_line = $lines[$i];
        chomp $curr_line;

        my @split = split(/;/, $curr_line);
        my $ltag = $split[0];

        # because the *.fa.intarna.sorted.csv is sorted the optimal
        # result always appears first in the table
        if (exists $opt_hash{$ltag}) {
            $subopt_hash{$ltag} = $curr_line;
        } else {
            $opt_hash{$ltag} = $curr_line;
        }
    }

    foreach my $key (keys %opt_hash) {
        print WRITETOOPT $opt_hash{$key} . "\n";
    }

    foreach my $key (keys %subopt_hash) {
        print WRITETOSUBOPT $subopt_hash{$key} . "\n";
    }

    close WRITETOOPT;
    close WRITETOSUBOPT;

}

