#!/usr/bin/env perl

use warnings;
use strict;
use Bio::SeqIO;
use Bio::DB::EUtilities;
use Getopt::Long;

my $acc_number;             # accession number of query
my $genome_file;            # genome file

GetOptions('acc=s'      =>  \$acc_number,
           'g|genome=s' =>  \$genome_file,);

&usage unless ($acc_number && $genome_file);
#die "file $genome_file does already exist\n" if (-e $genome_file);
exit(0) if (-e $genome_file); ## edit 2.0.3

# get GI for given genome accession number
my $factory = Bio::DB::EUtilities->new(-eutil => 'esearch',
                                       -db => 'nuccore',
                                       -field => 'Accession',
                                       -term => $acc_number,
                                       -tool => 'coprarna', ## edit 2.0.5.1 // changed to coprarna // also removed email
                                       -verbose => -1); # no warnings
my ($id) = $factory->get_ids;
die "no genome with accession number '$acc_number' found\n" unless ($id);

# get a summary
$factory->reset_parameters(-eutil => 'esummary',
                           -db => 'nuccore',
                           -id => $id);
# one individual DocSum objects per ID
my $ds = $factory->next_DocSum;
# flattened mode, iterate through all Item objects and read sequence length
my $seq_length;
while (my $item = $ds->next_Item('flattened'))  {
    # not all Items have content, so need to check...
    $seq_length = $item->get_content
        if ($item->get_content && $item->get_name =~ /Length/);
}

# download the full genbank file to a file (not retained in memory)
$factory = Bio::DB::EUtilities->new(-eutil => 'efetch',
                                    -db => 'nuccore',
                                    -id => $id,
                                    -tool => 'intarna_web',
                                    -email => 'intarna@informatik.uni-freiburg.de',
                                    -rettype => 'gbwithparts'); # gbwithparts to get full sequence
$factory->get_Response(-file => $genome_file); 


#####################################################################################################################################

sub usage {
   die("\nUsage: ./get_refseq_from_refid.pl -acc accession-number -g genome genome-file \n\n");
}

