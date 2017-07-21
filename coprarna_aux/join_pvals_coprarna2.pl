#!/usr/bin/env perl

use warnings;
use strict;
use Cwd 'abs_path'; ## edit 2.0.5.1

# combining opt_tags.clustered_trunc

# get absolute path
my $ABS_PATH = abs_path($0); ## edit 2.0.5.1
# remove script name at the end
# match all non slash characters at the end of the string
$ABS_PATH =~ s|[^/]+$||g; ## edit 2.0.5.1
my $PATH_COPRA_SUBSCRIPTS = $ABS_PATH;

my $tags_clust_file = $ARGV[0]; ## edit 2.0.5.1 

# first read the weights
my @weight_files = <weight_permutations/*zscore*>; ## edit 2.0.5.1 // fixed dir

my %weightfile_pointer_hash = (); # cat refids -> weight file

# make the refseqid set to file pointers
foreach (@weight_files) {
    my $ref_ids_rows = `awk -F';' '{print \$1}' $_ | tr '\n' '#'`;
    chomp $ref_ids_rows;
    chop $ref_ids_rows;
    #print $ref_ids_rows . "\n";
    my @split = split(/#/, $ref_ids_rows);
    @split = sort(@split); # sorting needs to be performed to keep everything comparable
    #foreach(@split) {
    #     print $_ . "\n";
    #}
    my $hash_key = join("#", @split);
    #print $hash_key . "+++\n";
    $weightfile_pointer_hash{$hash_key} = $_;
}

#foreach my $key (keys %weightfile_pointer_hash) {
#    print "$key:$weightfile_pointer_hash{$key}\n";
#}

# we need to permute over the tags.clustered file
# to find the correct weight file for each cluster
# and pass the cluster number and its affiliated
# weight file to the R script that does the combination
open (MYDATA, $tags_clust_file) or die ("Error: cannot open file $tags_clust_file\n");
    my @tags_clustered = <MYDATA>;
close (MYDATA);

my $curr_cluster_num = 0;
my @org_array = ();

# RCMD
open(RCMD, "| R --vanilla --slave"); # start R listener
print RCMD "data <- read.table(\"$tags_clust_file\", sep=';', header=T)\n"; # feed opt_tags.clustered_trunc

# basically we fill an array with refseq IDs
# to find the weightfile and pass the clusternumber
# and the weightfile to R
for (my $i=1;$i<scalar(@tags_clustered);$i++) {

    my $curr_line = $tags_clustered[$i];
    chomp $curr_line;
    my @split_curr_line = split(/;/, $curr_line);
    $curr_cluster_num = $split_curr_line[-1]; ## edit 2.0.5.1 // changed index to -1
    my $curr_org = $split_curr_line[1]; ## edit 2.0.5.1 // changed index to 1
    $curr_org =~ s/ncrna_//g;
    push(@org_array, $curr_org);

    my $next_line = "";
    my $next_line_cl_num = "";

    if(exists($tags_clustered[($i+1)])) {
        $next_line = $tags_clustered[($i+1)];
        chomp $next_line;
        my @split_next_line = split(/;/, $next_line);
        $next_line_cl_num = $split_next_line[-1]; ## edit 2.0.5.1 // changed index to -1
    } else { # exception for the last line in the file
        @org_array = () if (scalar(@org_array) < 3);
        next if (scalar(@org_array) < 3);
        @org_array = sort(@org_array);
        my $org_hash_key = join("#",@org_array);
        if(exists $weightfile_pointer_hash{$org_hash_key}) {
            # run combination of pvals for this cluster
            #print "Run R!!\n";
            print RCMD "cl_num <- $curr_cluster_num\n";
            print RCMD "weights <- read.table(\"$weightfile_pointer_hash{$org_hash_key}\", sep=';', header=F, colClasses=c(V1='character'))\n";
            print RCMD "source(\"$PATH_COPRA_SUBSCRIPTS/Hartung.R\")\n"; ## edit 2.0.4.2
            print RCMD "source(\"$PATH_COPRA_SUBSCRIPTS/combine_pval_hartung.R\")\n"; ## edit 2.0.4.1 // changed to hartung
        } else {
            # write a protest // this should never happen
            print "PROTEST!!\n";
        }
        last;
    }

    if($curr_cluster_num ne $next_line_cl_num) { # next line is different cluster
        @org_array = () if (scalar(@org_array) < 3);
        next if (scalar(@org_array) < 3);
        @org_array = sort(@org_array);
        my $org_hash_key = join("#",@org_array);
        if(exists $weightfile_pointer_hash{$org_hash_key}) {
            # run combination of pvals for this cluster
            #print "Run R!!\n";
            print RCMD "cl_num <- $curr_cluster_num\n";
            print RCMD "weights <- read.table(\"$weightfile_pointer_hash{$org_hash_key}\", sep=';', header=F, colClasses=c(V1='character'))\n";
            print RCMD "source(\"$PATH_COPRA_SUBSCRIPTS/Hartung.R\")\n"; ## edit 2.0.4.2
            print RCMD "source(\"$PATH_COPRA_SUBSCRIPTS/combine_pval_hartung.R\")\n"; ## edit 2.0.4.1 // changed to hartung
        } else {
            # write a protest // this should never happen
            print "PROTEST!!\n";
        }
        @org_array = ();
    }
}

close(RCMD);
