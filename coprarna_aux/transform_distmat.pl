#!/usr/bin/env perl 


# transforms the distmat(EMBOSS) distance matrix to a
# format readable by fneigbor(EMBOSS)

use warnings;
use strict;

my $distmat = $ARGV[0];
my $maxcluster = 0;
my $firstswitch = 1;

open(MYDATA, $distmat) or die("Error: cannot open file $distmat'\n");

my @distmatlines = ();

@distmatlines = <MYDATA>;

close MYDATA;
my $c = 0;
my %matrixhash = ();
my @orgarray = ();

foreach(@distmatlines) {
    $c++;
    if($c > 8) { # skip first lines
        my @splitline = split(/\t/,$_);
        my @cleansplit = ();
        for(my $i=0;$i<scalar(@splitline);$i++) { 
            if ($splitline[$i] =~ m/\d/) { 
                chomp $splitline[$i];
                $splitline[$i] =~ s/^\s+//; #remove leading spaces
                $splitline[$i] =~ s/\s+$//; #remove trailing spaces
                push(@cleansplit, $splitline[$i]);
            }
        }

         if ($firstswitch) {
             $maxcluster = scalar(@cleansplit) - 1;
             $firstswitch = 0;
             print "    $maxcluster\n";
         }
        my @splitorg = split(/\s/, ($cleansplit[(scalar(@cleansplit) - 1)]));
        push(@orgarray, lc($splitorg[0]));

        pop @cleansplit;  
        if(scalar(@cleansplit) < $maxcluster) {
            for(my $i=scalar(@cleansplit);$i<$maxcluster;$i++){
                unshift(@cleansplit, "0.00");
            } 
        }

        for(my $i=0;$i<scalar(@cleansplit);$i++) {
            $matrixhash{($c - 8)}{($i + 1)} = $cleansplit[$i];
        }
 
    }
}

foreach my $key (sort keys %matrixhash) {
    foreach my $key2 (sort keys %{$matrixhash{$key}}) {
        if($matrixhash{$key}{$key2} == 0) {
            $matrixhash{$key}{$key2} = $matrixhash{$key2}{$key};
        }
    }
}


foreach my $key (sort { $a <=> $b } keys %matrixhash) {
    print $orgarray[$key - 1] . "\t";
    foreach my $key2 (sort { $a <=> $b } keys %{$matrixhash{$key}}) {
        while(length($matrixhash{$key}{$key2}) < 5) {
            $matrixhash{$key}{$key2} = $matrixhash{$key}{$key2} . "0";
        }
        if ($key == $key2) { $matrixhash{$key}{$key2} = "0.000"; } 
        print "$matrixhash{$key}{$key2}  ";
    }
print "\n";
}

