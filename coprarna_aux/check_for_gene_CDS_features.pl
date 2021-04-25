#!/usr/bin/env perl

use strict;
use warnings;
# file handles
use IO::File;
use IO::Uncompress::Gunzip qw(gunzip $GunzipError);
# genbank file parsing
use Bio::SeqIO;

# this script checks if a gene feature is available but no CDS feature
# if yes, it throws an exception, causing CopraRNA to terminate
# this is a more sophisticated way of checking than checking solely for
# empty fasta files. issue came from sfd NC_018190 (no annotated orf)


# ./check_for_gene_CDS_features.pl NC_003197.gb > gene_CDS_exception.txt

my $refid = $ARGV[0];

my $CDScount = 0;
my $GENEcount = 0;

# need the full path of gunzip for Bio::SeqIO call

my $fileHandle = undef;
if ($refid =~ m/.+\.gz$/) {
	$fileHandle = new IO::Uncompress::Gunzip $refid or die "IO::Uncompress::Gunzip failed: $GunzipError\n";
} else {
	$fileHandle = IO::File->new($refid, "r");
}
my $seqin = Bio::SeqIO->new( -format => 'genbank', -fh => $fileHandle );

while( (my $seq = $seqin->next_seq()) ) {
    foreach my $sf ( $seq->get_SeqFeatures() ) {
        if( $sf->primary_tag eq 'CDS' ) {
            $CDScount++;
        }
        if( $sf->primary_tag eq 'gene' ) {
            $GENEcount++;
        }
    }
}

if ($CDScount == 0 and $GENEcount != 0) {
    print "CDS features are $CDScount while gene features are $GENEcount in $refid\n";
}



