#!/usr/bin/env perl 

use warnings;
use strict;

my $pythonEnrichment = $ARGV[0];

my $clusterCounter = 1;
my $header = "Category\tTerm\tCount\t%\tPvalue\tGenes\tList Total\tPop Hits\tPop Total\tFold Enrichment\tBonferroni\tBenjamini\tFDR\n";
my $firstSwtich = 1;

open(MYDATA, $pythonEnrichment) or die("\nError: cannot open file $pythonEnrichment at make_enrichment_table_from_py_output.pl\n\n");
    my @pythonEnrichmentLine = ();
    @pythonEnrichmentLine = <MYDATA>;
close MYDATA;

for (my $i=0;$i<scalar(@pythonEnrichmentLine);$i++) {

    my $currLine = $pythonEnrichmentLine[$i];
    chomp $currLine;
    if($currLine =~ m/score=/) {
        my @split = split(/=/, $currLine);
        my $enrichmentScore = $split[1];
        if ($firstSwtich) { # no newline in the first iteration
            print "Annotation Cluster $clusterCounter\tEnrichment Score:  $enrichmentScore\n"; 
            $firstSwtich = 0;
        } else {
            print "\nAnnotation Cluster $clusterCounter\tEnrichment Score:  $enrichmentScore\n"; 
        }
        print $header;
        $clusterCounter++;
    } 

    if ($currLine =~ m/afdr=/) {

        my %nextPrintLineHash = (); # hash is more convenient because we can directly access
        my @split = split(/=/, $currLine);
        $nextPrintLineHash{$split[0]}=$split[1];
        for (my $j=($i + 1);$j<scalar(@pythonEnrichmentLine);$j++) { # afdr line is already take care of so we do j = i + 1
            my $currSubline = $pythonEnrichmentLine[$j];
            chomp $currSubline;
            if ($currSubline =~ m/afdr=/) { # break loop 
                last;
            } else { # fill hash
                my @split = split(/=/, $currSubline);
                $nextPrintLineHash{$split[0]}=$split[1];
            }
        }
        # print line
        my $printLine = $nextPrintLineHash{"categoryName"} . "\t" . $nextPrintLineHash{"termName"} . "\t" . $nextPrintLineHash{"listHits"} . "\t" . $nextPrintLineHash{"percent"} . "\t" . $nextPrintLineHash{"ease"} . "\t" . $nextPrintLineHash{"geneIds"} . "\t" . $nextPrintLineHash{"listTotals"} . "\t" . $nextPrintLineHash{"popHits"} . "\t" . $nextPrintLineHash{"popTotals"} . "\t" . $nextPrintLineHash{"foldEnrichment"} . "\t" . $nextPrintLineHash{"bonferroni"} . "\t" . $nextPrintLineHash{"benjamini"} . "\t" . $nextPrintLineHash{"afdr"};
        $printLine =~ s/"//g;
        print $printLine . "\n"; 
    }

}


