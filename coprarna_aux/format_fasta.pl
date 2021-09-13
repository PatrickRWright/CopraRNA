#!/usr/bin/env perl

use strict;
use warnings;
use Bio::SeqIO;

# this script makes sure that the sequence is only in one line

my $fasta = $ARGV[0];

my $seqio = Bio::SeqIO->new(-file => $fasta, -format => "fasta");

while (my $seq = $seqio->next_seq) {
    my $curr_seq = $seq->seq();
    my $curr_id = $seq->id();
    print ">$curr_id\n$curr_seq\n";
}

