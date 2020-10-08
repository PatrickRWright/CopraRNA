#!/usr/bin/env perl

use warnings;
use strict;
use Bio::SeqIO;
use Bio::DB::EUtilities;
use Getopt::Long;

use Cwd 'abs_path'; 

my $acc_number;             # accession number of query
my $genome_file;            # genome file
my $cores;
my $genomePath;

GetOptions('acc=s'      =>  \$acc_number,
           'g|genome=s' =>  \$genome_file,
	   'c=s'	=>  \$cores,
	   'gPath=s'	=>  \$genomePath,);

&usage unless ($acc_number && $genome_file);

# check if file is already present
exit(0) if (-e $genome_file);

# get absolute path
my $ABS_PATH = abs_path($0); 
# remove script name at the end
# match all non slash characters at the end of the string
$ABS_PATH =~ s|[^/]+$||g; 
my $PATH_COPRA_SUBSCRIPTS = $ABS_PATH;

# get GI for given genome accession number
my $factory = Bio::DB::EUtilities->new(-eutil => 'esearch',
                                       -db => 'nuccore',
                                       -field => 'Accession',
                                       -term => $acc_number,
                                       -tool => 'coprarna',
                                       -email => 'rna@informatik.uni-freiburg.de',
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

my $dwnCall = "python ". $PATH_COPRA_SUBSCRIPTS  . "download_genomes.py -g $acc_number -c $cores -p $genomePath/";
system $dwnCall;

sub usage {
   die("\nUsage: ./get_refseq_from_refid.pl -acc accession-number -g genome-file \n\n");
}

