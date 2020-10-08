#!/usr/bin/env perl

use strict;
use warnings;
# file handles
use IO::File;
use IO::Uncompress::Gunzip qw(gunzip $GunzipError);
# genbank file parsing
use Bio::SeqIO;
use Parallel::ForkManager;

# parses 16s rRNA sequences from genbank files
my $cores = $ARGV[0];

my $pm = new Parallel::ForkManager($cores);

for (my $i=1;$i<scalar(@ARGV);$i++) {
    my @splitargv = split(/,/, $ARGV[$i]);

    # save splitiargv[0] and put this in as header. remove the .gb !!
    my $MainRefID = $splitargv[0];
    chomp $MainRefID;
    if ($MainRefID =~ m/^(.+)\.gb(\.gz)?$/) {
       $MainRefID = $1;
    } else {
       die("\n parse_16s_from_gbk.pl : given genome file does not end in '.gb' or '.gb.gz'\n\n");
    }

    $pm->start and next;    
    for (my $i=0; $i<scalar(@splitargv); $i++) {

	my $genomeFile = $splitargv[$i];
    my $fileHandle = undef;
    if ($genomeFile =~ m/.+\.gz$/) {
    	$fileHandle = new IO::Uncompress::Gunzip $genomeFile or die "IO::Uncompress::Gunzip failed: $GunzipError\n";
    } else {
    	$fileHandle = IO::File->new($genomeFile, "r");
    }
        my $seqin = Bio::SeqIO->new( -format => 'genbank', -fh => $fileHandle);

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
    $pm->finish;
}
$pm->wait_all_children;
