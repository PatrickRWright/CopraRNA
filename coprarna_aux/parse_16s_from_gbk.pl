#!/usr/bin/env perl

use strict;
use warnings;

use Bio::SeqIO;


foreach(@ARGV) {

my @splitargv = split(/,/, $_);
# save split[0] and put this in as header. remove the .gb !!

my $MainRefID = $splitargv[0];
chomp $MainRefID;
chop $MainRefID;
chop $MainRefID;
chop $MainRefID;


for (my $i=0; $i<scalar(@splitargv); $i++) {
my $seqin = Bio::SeqIO->new( -format => 'genbank', -file => $splitargv[$i]);


while( (my $seq = $seqin->next_seq()) ) {
    foreach my $sf ( $seq->get_SeqFeatures() ) {
        if( $sf->primary_tag eq 'rRNA' ) {
            my $product = "";
            if ($sf->has_tag("product")) {
                my @productlist = $sf->get_tag_values("product");
                $product = $productlist[0];
            }
            my $id = $seq->display_id;
            if ($product =~ m/16S/i or $product =~ m/Small subunit ribosomal RNA/i) {
                print ">$MainRefID\n";
                print $sf->spliced_seq->seq, "\n";
                $i = 10000;
                last;
        }
    }
  }
if ($i == 10000) { last; }
}
}


}
