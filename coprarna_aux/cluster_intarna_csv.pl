#!/usr/bin/env perl

use strict;
use warnings;
use List::MoreUtils qw(uniq);
## run in a CopraRNA directory with
## *_opt.intarna.csv files and cluster.tab
## to create a clustering of IntaRNA predictions 

my @final_csv_files = <*_opt.intarna.csv>; 
my $intarna_csv_header = "";

my %ltag_line_hash = (); # locus_tag -> intarna.csv line

# create a hash entry for every line in *_opt.intarna.csv files
# this hash will be later accessed using the lines from cluster.tab
foreach ( @final_csv_files ) {

    open(MYDATA, $_) or die("\nError: cannot open file $_ in cluster_intarna_csv.pl\n\n"); ## 

        my @final_csv_lines = <MYDATA>;

		# copy header information if not already done
		$intarna_csv_header = $final_csv_lines[0] if ( $intarna_csv_header eq "" );

		# store intarna output 
        for (my $i=1;$i<scalar(@final_csv_lines);$i++) { 
            my $curr_line = lc($final_csv_lines[$i]);
            chomp $curr_line;
            my @split = split(/;/,$curr_line);
            my $ltag = $split[0];
            $ltag_line_hash{$ltag} = $curr_line;
        }
    close MYDATA;
}

# print new output header
print "$intarna_csv_header;clusternumber\n";

my $homologlist = "cluster.tab";

open(MYDATA, $homologlist) or die("\nError: cannot open file $homologlist in cluster_intarna_csv.pl\n\n");
    my @cluster_tab_lines = <MYDATA>;
close MYDATA;

# go through cluster.tab and produce tags.clustered file
for (my $i=1;$i<scalar(@cluster_tab_lines);$i++) {
    
    my $curr_line = lc($cluster_tab_lines[$i]);
    chomp $curr_line;
    # split line from cluster.tab -> one entry per organism
    my @split_line = split(/\t/, $curr_line);
    
    for (my $j=1;$j<scalar(@split_line);$j++) { # skip first column 
        my $pvalue = 100; # need this to save the homolog with the lowest pvalue 
                          # if one org. contains more than one homolog
        my $print_line = "";
        my $curr_cell = $split_line[$j];
        my @split_cell = split(/\s/,$curr_cell);
        foreach(@split_cell) { # if there is only one entry this loops only once
           $_ =~ s/\(\d+\)//g; # remove digit in brackets at the end 
           $_ =~ s/\w+://g; # remove pseudo KEGG id and colon
           if (exists $ltag_line_hash{$_}) {
               my @split_line = split(/;/, $ltag_line_hash{$_});
               # $split_line[35] is the current p-value
               if ($split_line[35] < $pvalue) {
                   $print_line = $ltag_line_hash{$_};
                   $pvalue = $split_line[35];
               } 
           } 
        }
        print $print_line . ";$i\n" if ($print_line); # only if its not empty
    }
}

