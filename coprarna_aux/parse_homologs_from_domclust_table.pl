#!/usr/bin/env perl

use strict;
use warnings;

use List::MoreUtils qw(uniq);

# get additional homologs from cluster.tab for the organism of interest

my $finallist = $ARGV[0];
my $clustertab = $ARGV[1];
my @clustertablines = ();
my $line;
my %homologtaghash = ();

open(MYDATA, $clustertab) or die("Error: cannot open file $clustertab'\n");

@clustertablines = <MYDATA>;

close MYDATA;


foreach(@clustertablines) {
    if (not $_ =~ m/^#|^running/) { # skip
        my @splitlinearray = split(/\t/, $_); # split each line by tabs
        for (my $i=1;$i<scalar(@splitlinearray);$i++) {
            my @homologtaglist = ();
            my @splitthesplitted = split(/\s/, $splitlinearray[$i]); # split the tab splitted values by whitespaces
            foreach (@splitthesplitted) {
                #print "$_\n";
                if ($_ =~ m/\w{3}:(.*)/) { # filter the locus tags
                    chomp $1;
                    my $temp  = $1;
                    if ($temp =~ m/(.*?)\(\d+\)/) { # remove number at the back e.g. only keep locus tag
                        chomp $1;
                        $temp = $1;
                    }
                    push(@homologtaglist, lc($temp));
                }
            }
            foreach(@homologtaglist) {
                push(@{ $homologtaghash{ $_ }}, @homologtaglist); # make each entry in the list point to the list via hash
            }
        }
    }
}

my @finallistlines = ();

open(MYDATA, $finallist) or die("Error: cannot open file $finallist'\n");

@finallistlines = <MYDATA>;

close MYDATA;

chomp $finallistlines[0];

print $finallistlines[0];
print ",Additional homologs\n";

for(my $i=1; $i<scalar(@finallistlines); $i++) {
    chomp $finallistlines[$i];
    print $finallistlines[$i];
    my @splitline = split(/,/, $finallistlines[$i]); # split each line by comma
    if($splitline[1] =~ m/^(\S+?)\(\S+GeneID\S+\)/) { # get the locus tag for the first organism (e.g. org of interest)
        chomp $1;
        print ",";
        if (exists $homologtaghash{ $1 }) {
            my @temphomologtagarray = @{$homologtaghash{$1}};
            @temphomologtagarray = uniq(@temphomologtagarray);
            foreach(@temphomologtagarray) { 
                if($_ ne $1) {
                    print $_ . " "; 
                }
            }
        }
        print "\n";
    } else { print ",\n" ; }
}

