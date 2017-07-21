#!/usr/bin/env perl

use strict;
use warnings;

use Bio::SeqIO;
use Cwd 'abs_path'; ## edit 2.0.5.1

# ./get_CDS_from_gbk.pl NC_003197.gb > stm.fas

# gets a refseq file and creates a fasta
# with all CDS AA sequences in the following format:
# >eco:b0001
# ISLLGQCPDSINAALICRGEKMSIAIMAGVLEARGHNVTVIDPVEKLL

# >KEGGID:locus_tag
# Aminoacid sequence

# get absolute path
my $ABS_PATH = abs_path($0); ## edit 2.0.5.1
# remove script name at the end
# match all non slash characters at the end of the string
$ABS_PATH =~ s|[^/]+$||g; ## edit 2.0.5.1
my $PATH_COPRA_SUBSCRIPTS = $ABS_PATH;

my $kegg2refseq = $PATH_COPRA_SUBSCRIPTS . "kegg2refseqnew.csv";

open(MYDATA, $kegg2refseq) or die("\nError: cannot open file $kegg2refseq in get_CDS_from_gbk.pl\n\n");
    my @kegg2refseq = <MYDATA>;
close MYDATA;

my $refid = $ARGV[0];
my $refidchopped = $refid;
chop $refidchopped;
chop $refidchopped;
chop $refidchopped;

my $keggcode = "";

foreach(@kegg2refseq) {
    if($_ =~ m/$refidchopped/) {
        my @split = split(/\s/,$_);
        $keggcode = $split[0];
        last;
    }
}

my $seqin = Bio::SeqIO->new( -format => 'genbank', -file => $refid);

while( (my $seq = $seqin->next_seq()) ) {
    foreach my $sf ( $seq->get_SeqFeatures() ) {
        if( $sf->primary_tag eq 'CDS' ) {
            my $ltag = "";
            my $protein = "";
            if ($sf->has_tag("locus_tag") and $sf->has_tag("translation")) {
                my @ltaglist = $sf->get_tag_values("locus_tag");
                $ltag = $ltaglist[0];
                my @translationlist = $sf->get_tag_values("translation");
                $protein = $translationlist[0];
            
                chomp $ltag;
                print ">$keggcode:$ltag\n";
                chomp $protein;
                print "$protein\n";
         }
      }
   }
}

