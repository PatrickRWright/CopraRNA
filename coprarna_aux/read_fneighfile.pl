#!/usr/bin/env perl 

use warnings;
use strict;

## edit 2.0.4 // major changes to adjust the entire script
# to UPGMA tree -> no more triple root

my $treefile = $ARGV[0];
my $fneighfile = $ARGV[1];

open(MYDATA, $treefile) or die("Error: cannot open file $treefile'\n");
my @treefilelines = ();
@treefilelines = <MYDATA>;
close MYDATA;


my $treelength = 0;

# split the newick formated tree up to get the length of the 
# entire tree
foreach(@treefilelines) {
    chomp $_;
    my @splitarray = split(/[:\(\),]/, $_);
    foreach(@splitarray) {
        if($_ =~ m/^(\d+\.\d+)$/) {  ## edit 1.3.0 fixed regex for double digit branchlengths (i.e. 95.0381)
            $treelength = $treelength + $1;
        }
    }
}

open(MYDATA, $fneighfile) or die("Error: cannot open file $fneighfile'\n");
my @fneighfilelines = ();
    @fneighfilelines = <MYDATA>;
close MYDATA;

my $theorgcnt = 0;

# get the total count of organisms
# in the tree
foreach(@fneighfilelines) {
    if($_ =~ m/(\d+)\sPopulations/) {
        $theorgcnt = $1;
        chomp $theorgcnt;
    }
}

my @lastlines = ();
my $switch = 0;

# get the lines that contain the information
# on the branch lengths
foreach(@fneighfilelines) {
    chomp $_;
    if($_ =~ m/From\s+To\s+Length\s+Height/) { ## edit 2.0.4
        $switch = 1;  
    }
    if($switch) {
        if($_ =~ m/\d/) {
            push @lastlines, $_;
        }
    }
}

my %weighthash = ();

my $lastlinecnt = scalar(@lastlines);


# find final root in lastlines and add a pseudoline
my $final_root = 1000000000000000;

# first col in last lines
my @col1 = ();
# second col in last lines
my @col2 = ();

foreach(@lastlines) {
    $_ =~ s/^\s+//;
    $_ =~ s/\s+$//;
    my @split = split(/\s+/,$_);
    if ($split[0] =~ m/^\d+$/) {
        push(@col1,$split[0]);
    }
    if ($split[1] =~ m/^\d+$/) {
        push(@col2,$split[1]);
    }
}

my %col2_hash = ();

foreach (@col2) {
    $col2_hash{$_} = 1;
}

foreach (@col1) {
    unless (exists $col2_hash{$_}) {
        $final_root = $_;
    }
}

my $pseudoline = "999 $final_root 0 0";
push(@lastlines, $pseudoline);
#### up to here we simply prepare the data for calculating the weights


my %tree = ();

foreach(@lastlines) {
    my @split = split(/\s+/,$_);
    $tree{$split[0]}{$split[1]} = $split[2];
}

my %flippedtree = ();

foreach(@lastlines) {
    my @split = split(/\s+/,$_);
    $flippedtree{$split[1]}{$split[0]} = $split[2];
}

my @refseqarray = ();

for my $key (keys %tree) {
    for my $key2 (keys %{$tree{$key}}) {
        if($key2 =~ m/n[zc]/) { ## edit 2.0.2
            push @refseqarray, $key2;
        }
    }
}

my %subtreelengthhash = ();

# edit the substring here ## edit 1.2.3
foreach(@lastlines) {
    my @split = split(/\s+/, $_);
    my $nodenum = $split[0];
    $subtreelengthhash{$nodenum} = &subtreelength($nodenum, \%tree);
}

foreach my $refid (@refseqarray) {
    for my $root (keys %{$flippedtree{$refid}}) {
        for my $rootlength (keys %{$flippedtree{$root}}) {
            $weighthash{$refid} = &calctree($root, $flippedtree{$refid}{$root}, \%flippedtree, $flippedtree{$root}{$rootlength}, \%subtreelengthhash, $treelength, $theorgcnt);
        }
    }
}


# subtreelength does not contain the
# root branch
sub subtreelength 
{
    (my $nodenum, my $hashref) = @_;

    my $subtreelength = 0;
    my %localhash = %{$hashref};
    for my $key (keys %{$localhash{$nodenum}}) {
        if($key =~ m/n[zc]/) { ## edit 2.0.2 ## if an internal node points at a leaf we are done
            $subtreelength = $subtreelength + $localhash{$nodenum}{$key};
        } else { ## if an internal node points at an internal node we need to continue
            $subtreelength = $subtreelength + $localhash{$nodenum}{$key} + &subtreelength($key, \%localhash);
        }
    }
    return $subtreelength;
}


sub calctree
{ ## edit 1.2.8 added orgcnt to input arguments for division by zero error  
    (my $root, my $branchlength, my $treeref, my $rootlength, my $subtreelengthsref, my $treelength, my $orgcnt) = @_;
    my $result = 1;
    my %localtree = %{$treeref};
    my %localsubtree = %{$subtreelengthsref};
    for my $nextroot (keys %{$localtree{$root}}) {
            my $nextsubtree = $localsubtree{$root}+$localtree{$root}{$nextroot};
            my $nextrootlength = 0;
            for my $thenextroot (keys %{$localtree{$nextroot}}) {
                $nextrootlength = $localtree{$nextroot}{$thenextroot};
            }
            unless ($treelength == 0) { ## edit 1.2.8 fix divide by zero bug for treelength 0
                $result = ($branchlength + ($localtree{$root}{$nextroot}/2))/($localsubtree{$root}+$localtree{$root}{$nextroot}) *
                &calctree($nextroot, $nextsubtree, \%localtree, $nextrootlength, \%localsubtree, $treelength, $theorgcnt);
            } else {
                $result = 1/$orgcnt;
            }
    }
    return $result;
}

# print result
for my $entry (keys %weighthash) {
    print "$entry;$weighthash{$entry}\n";
}

