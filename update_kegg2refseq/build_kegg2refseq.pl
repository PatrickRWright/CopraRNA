#!/usr/bin/env perl

use warnings;
use strict;
use List::MoreUtils qw(uniq);

# usage
# ../build_kegg2refseq.pl

system "wget ftp://ftp.ncbi.nih.gov/genomes/GENOME_REPORTS/prokaryotes.txt";

open MYDATA, "prokaryotes.txt" or die("\nERROR: Can't open prokaryotes.txt at build_kegg2refseq.pl\n\n");
    my @genomeInfo = <MYDATA>;
close MYDATA;

# open writing file handles
open (HP, '>CopraRNA_available_organisms.tmp');
open (New, '>kegg2refseqnew.csv');

my %rdmStringHash = ();
my %printedHash = ();

foreach my $line (@genomeInfo) {
    my $printLineHP = "";
    my $printLineNew = "";
    my $switch = 1;
    if ($line =~ m/NC_\d{6}|NZ_.+/) { ## edit 2.0.2
        my $rdmString = &generate_random_string(4);
        
        # make sure no duplicate letter codes are present
        while ($switch) {
            if (exists $rdmStringHash{$rdmString}) {
                $rdmString = &generate_random_string(4);
            } else {
                $switch = 0;
                $rdmStringHash{$rdmString} = "exists";
                $printLineNew = $rdmString . "\t";
            }
        }

        my @splitLine = split(/\t/, $line);
        
        # specify ID cell
        my $ID_cell = $splitLine[8]; ## edit 2.0.5.1
        my @split_ID_cell = split(/;/, $ID_cell);       
 
        foreach my $entry (@split_ID_cell) {
            if ($entry =~ m/(NC_\d{6}|NZ_.+)/) { ## edit 2.0.2
                my $RID = $1;
                chomp $RID;
                # remove .1 at the end of new ids
                $RID =~ s/\.\d+//g; ## edit 2.0.2
                # remove all after '/'
                $RID =~ s/\/.+//g; ## edit 2.0.5.1
                $printLineHP = $printLineHP . $RID . " ";
                $printLineNew = $printLineNew . $RID . " ";
            }
        }
        $printLineHP =~ s/\s+$//;
        $printLineNew =~ s/\s+$//;
        my $orgName = $splitLine[0];
        $orgName =~ s/^\s+//;
        $orgName =~ s/\s+$//;
        $orgName =~ s/\s+/_/g;
        $printLineHP = $printLineHP . "\t" . $orgName . "\n";
        $printLineNew = $printLineNew . "\n";
 
        # don't print same line twice
        print HP $printLineHP unless (exists $printedHash{$printLineHP});
        print New $printLineNew unless (exists $printedHash{$printLineHP});

        $printedHash{$printLineHP} = "exists";
    }
}

# close writing file handles
close (HP);
close (New);

system "bash ../add_date_omit_incompatible.sh";

sub generate_random_string
    {
	my $length_of_randomstring=shift;# the length of 
        my @chars=('a'..'z');
	my $random_string;
	foreach (1..$length_of_randomstring) 
	{
            $random_string.=$chars[rand @chars];
	}
	return $random_string;
    }

# usage
# my $random_string=&generate_random_string(4);

