#!/usr/bin/env perl 

use warnings;
use strict;

my $inputTable = $ARGV[0]; ## edit 2.0.5.1
my $cop2 = $ARGV[1]; ## edit 2.0.5.1 // 1=yes 0=no

open(MYDATA, $inputTable) or die("\nError: cannot open $inputTable at get_amount_sampled_values_and_add_to_table.pl\n\n");
    my @resultLines = <MYDATA>;
close MYDATA;

# header
chomp $resultLines[0];
print $resultLines[0] . ",Amount sampled\n";

for (my $i=1;$i<scalar(@resultLines);$i++) {

    chomp $resultLines[$i];
    print $resultLines[$i];
    my $emptyCount = 0;
    my @split = split(/,/,$resultLines[$i], -1);
    # remove annotation and additional homologs
    pop @split;
    pop @split;

    foreach(@split) {
        $emptyCount++ if ($_ eq ""); 
    }
    # workaround for the webserver output // CopraRNA2 never samples
    $emptyCount = 0 if ($cop2); ## edit 2.0.5.1
    print ",$emptyCount\n";
}

