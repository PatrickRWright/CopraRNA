#!/usr/bin/env perl

use strict;
use warnings;

# this script produces the auxilliary enrichment 
# for the organism of interest

# read in main termClusterReport
open(DATA,"termClusterReport_cop1.txt") or die ("\nError: cannot open termClusterReport_cop1.txt at find_single_specific_targets_in_termCluster.pl\n\n");
    my @mainClustering = <DATA>;
close DATA;

# read in single organism termClusterReport
open(DATA,"IntaRNA_chartReport_grepped.txt") or die ("\nError: cannot open IntaRNA_chartReport_grepped.txt at find_single_specific_targets_in_termCluster.pl\n\n");
    my @singleOrgClustering = <DATA>;
close DATA;

print "term;predicted_targets EntrezGeneID(locus_tag|gene_name|target_start|target_end)\n";

# build hashes from enrichments (term -> gene Ids) 
my %mainClusteringHash = ();
my %singleOrgClusteringHash = ();

my $currEnrichmentScore = 0; ## edit 2.0.2

foreach (@mainClustering) {
    if ($_ =~ m/Enrichment\s+Score:\s+(\d+\.\d+)/) { ## edit 2.0.2
        $currEnrichmentScore = $1;
    }
    my @tabSplit = split(/\t/, $_);
    if (exists $tabSplit[1] and exists $tabSplit[5]) { # skip "Enrichment Score:" and empty lines
        # only for enrichment scores >= 1
        $mainClusteringHash{$tabSplit[1]} = $tabSplit[5] unless ($tabSplit[1] eq "Term" or $currEnrichmentScore < 1); ## edit 2.0.2
    }
}

for (my $i=0;$i<scalar(@singleOrgClustering);$i=$i+2) { ## edit 2.0.3.1

    my $gidLine = $singleOrgClustering[$i];
    chomp $gidLine;
    my $termLine = $singleOrgClustering[$i+1];
    chomp $termLine;

    my @splitGidLine = split(/=/, $gidLine);
    my $GIDs = $splitGidLine[1];
    my @splitTermLine = split(/=/, $termLine);
    my $term = $splitTermLine[1];

    $singleOrgClusteringHash{$term} = $GIDs;
}

# compare cluster memberships
foreach my $key (keys %mainClusteringHash) {
    my $printLine = "";
    my @splitMainClustering = split(/,/, $mainClusteringHash{$key});
    my @splitSingleClustering = ();
    # only for those clusters that exist in both CopraRNA and IntaRNA enrichment
    if (exists $singleOrgClusteringHash{$key}) {
        @splitSingleClustering = split(/,/, $singleOrgClusteringHash{$key});

        my %compareHash = ();    

        foreach (@splitMainClustering) { 
            $_ =~ s/\s+//g; 
            $compareHash{$_} = "exists"; 
        }
        for (my $i=0;$i<scalar(@splitSingleClustering);$i++) { 
            $splitSingleClustering[$i] =~ s/\s+//g;
            unless (exists $compareHash{$splitSingleClustering[$i]}) {
                $printLine = $printLine . $splitSingleClustering[$i];
                my $IntaRNA_ncbi_line = `grep $splitSingleClustering[$i] intarna_websrv_table_ncbi.csv`;
                chomp $IntaRNA_ncbi_line;
                my @split = split(/;/,$IntaRNA_ncbi_line);
                my $ltag = $split[0];
                my $geneName = $split[36]; ## edit 2.0.5.1 
                my $tarStart = $split[8]; ## edit 2.0.5.1
                my $tarStop = $split[9]; ## edit 2.0.5.1
                my $added_info = "(" . $ltag . "|" . $geneName . "|" . $tarStart . "|" . $tarStop . ");";
                $printLine = $printLine . $added_info;
            }
        }
        chop $printLine; # remove trailing ;
        $printLine = "$key;" . $printLine . "\n" unless (length($printLine) == 0);
        print $printLine unless (length($printLine) == 0);
    }
}


