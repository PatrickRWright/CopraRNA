#!/usr/bin/env perl

use strict;
use warnings;
# file handles
use IO::File;
use IO::Uncompress::Gunzip qw(gunzip $GunzipError);
# genbank file parsing
use Bio::SeqIO;

# replaces get_genname_genid_note_from_gbk_opt.pl
# and now works with tags.clustered files as
# input to annotate IntaRNA result parts

my $finallist = $ARGV[0];
my $tags_clusterd = $ARGV[1];
my %ltaggennamehash = ();
my $ltag = "";
my $genname = "";
my $note = "";
my $genid = "";
my @ids = ();
my $columncount = 1;
my %ltagcolumnhash = ();
my %ltaggenidhash = ();
my $argofinterestswitch = 1;

print "p-value,";

for(my $i=2;$i<=(scalar(@ARGV)-1);$i++) {
        $columncount++;
        my @splitarg = split(/,/, $ARGV[$i]);
        if ($splitarg[0] =~ m/(N[ZC]_.+?)\.gb(\.gz)?/) {
            my $temp = $1;
            chomp $temp;
            # print RefSeq IDs header
            print "$temp,";
        }
        foreach my $genomeFile (@splitarg) {
          my $fileHandle = undef;
          if ($genomeFile =~ m/.+\.gz$/) {
          	$fileHandle = new IO::Uncompress::Gunzip $genomeFile or die "IO::Uncompress::Gunzip failed: $GunzipError\n";
          } else {
          	$fileHandle = IO::File->new($genomeFile, "r");
          }

        my $seqio_object = Bio::SeqIO->new(-fh => $fileHandle, -format => 'genbank');
        my $seq_object = $seqio_object->next_seq;

        for my $feat_object ($seq_object->get_SeqFeatures) {
            push @ids, $feat_object->get_tag_values("locus_tag") if ($feat_object->has_tag("locus_tag"));
            if (exists $ids[0]) {
            } else { $ids[0] = "N/A"; }
            push @ids, $feat_object->get_tag_values("gene") if ($feat_object->has_tag("gene"));
            if (exists $ids[1]) {
            } else { $ids[1] = "N/A"; }
            push @ids, $feat_object->get_tag_values("db_xref") if ($feat_object->has_tag("db_xref"));
            if (exists $ids[2]) {
            } else { $ids[2] = "N/A"; }


            if (exists $ids[0]) {
                $ltag = lc($ids[0]);
            } else {
                $ltag = "N/A";
            }
            if (exists $ids[1]) {
                $genname = $ids[1];
            } else {
                $genname = "N/A";
            }
            if (exists $ids[2]) {
                for(my $i=2; $i < scalar(@ids); $i++) {
                    if ($ids[$i] =~ m/GeneID:/) {
                        $genid = $ids[$i];
                    }
                }
            } else {
               $genid = "N/A";
            }

            chomp $ltag;
            chomp $genname;
            chomp $genid;
            if($ltag ne "N/A") {
                $ltagcolumnhash{$ltag} = $columncount; 
            }
            if ($ltag ne "N/A" and $genname ne "N/A") {
                $genname =~ s/,/;/g;
                $genname =~ s/\(//g;
                $genname =~ s/\)//g; 
                $ltaggennamehash{$ltag} = $genname;
            }
            if ($ltag ne "N/A" and $genid ne "N/A" and $genid =~ m/GeneID:/) {
                $ltaggenidhash{$ltag} = $genid;
            }
            $ltag = "";
            $genname = "";
            $genid = "";
            @ids = ();
        }
    }
$argofinterestswitch = 0;
}

print "Annotation";

print "\n";

open(MYDATA, $finallist) or die("\nError: cannot open file " . $finallist . " at annotate_raw_output.pl\n\n");
    my @finallistlines = <MYDATA>;
close MYDATA;

my %pvalenergyhash = ();
my $pvalenergy = '';

open (DATA, $tags_clusterd) or die ("\nError: cannot open file " . $tags_clusterd . " at annotate_raw_output.pl\n\n");
    my @datalines = <DATA>;
close DATA;

foreach my $line (@datalines) {
    my @splitarg = split(/;/, $line);
    foreach(@splitarg) { chomp $_; }
    $pvalenergy = "|" . $splitarg[14] . "|" . $splitarg[-2] . "|" . $splitarg [8] . "|" . $splitarg [9] . "|" . $splitarg [10] . "|" . $splitarg [11];
    chomp $pvalenergy;
    $pvalenergyhash{lc($splitarg[0])} = $pvalenergy;
}

my %annotationhash = (); #locustags -> annotation

my @splitThirdArgv = split(/,/,$ARGV[2]);

foreach my $genomeFile (@splitThirdArgv) {
    my $fileHandle = undef;
    if ($genomeFile =~ m/.+\.gz$/) {
    	$fileHandle = new IO::Uncompress::Gunzip $genomeFile or die "IO::Uncompress::Gunzip failed: $GunzipError\n";
    } else {
    	$fileHandle = IO::File->new($genomeFile, "r");
    }
    my $in  = Bio::SeqIO->new(-fh => $fileHandle, '-format' => 'genbank');
    while ( my $seq = $in->next_seq() ) {
        foreach my $sf ( $seq->get_SeqFeatures() ) {
            if( $sf->primary_tag eq 'CDS' ) {
            # add exception for "hypothetical protein" annotation
            my $product = "";
            my @products = ();
            if ($sf->has_tag("locus_tag") and $sf->has_tag("product")) {
                @products = $sf->get_tag_values("product");
                $product = $products[0];
            }

            if ($sf->has_tag("locus_tag") and $sf->has_tag("product") and (not $product =~ m/hypothetical protein/)) { # product
                my @ltaglist = $sf->get_tag_values("locus_tag");
                my $ltag = $ltaglist[0];
                my @prodlist = $sf->get_tag_values("product");
                my $prod = $prodlist[0];
                $annotationhash{$ltag} = $prod;
            } 
            elsif ($sf->has_tag("locus_tag") and $sf->has_tag("note")) { # note
                my @ltaglist = $sf->get_tag_values("locus_tag");
                my $ltag = $ltaglist[0];
                my @notelist = $sf->get_tag_values("note");
                my $note = $notelist[0];
                $annotationhash{$ltag} = $note;
            }
            elsif ($sf->has_tag("locus_tag") and $sf->has_tag("function")) { # function
                my @ltaglist = $sf->get_tag_values("locus_tag");
                my $ltag = $ltaglist[0];
                my @funclist = $sf->get_tag_values("function");
                my $func = $funclist[0];
                $annotationhash{$ltag} = $func;
            }
            }
        }
    }    
}

foreach (@finallistlines) {
    my $commacount = 1;
    my @splitlines = split(/;/, $_);
    pop(@splitlines);
    print $splitlines[0] . ",";
    shift(@splitlines);
    my @newsplit = ();
    my @newsplit2 = ();
    
    foreach (@splitlines) { 
       $newsplit[$ltagcolumnhash{$_}] = $_; 
    } 
 
    foreach (@newsplit) {
       if (defined $_) {
       push (@newsplit2, $_);
       }
    }

  

    my $nscount = 0;
    if ($ltagcolumnhash{$newsplit2[0]} ne 2) {

        for(my $k=2; $k<$ltagcolumnhash{$newsplit2[0]};$k++) {
            $commacount++;
            print ",";
        }
    }    

    foreach(@newsplit2) {
        print $newsplit2[$nscount];
        if (exists $ltaggennamehash{$newsplit2[$nscount]}) { # this is the locus tag

            print "(" . $ltaggennamehash{$newsplit2[$nscount]} . $pvalenergyhash{$newsplit2[$nscount]} . "|";
        
        } else {

            print "(N/A" . $pvalenergyhash{$newsplit2[$nscount]} . "|";
        }
       
        if (exists $ltaggenidhash{$newsplit2[$nscount]}) { # this is the gene id
            print $ltaggenidhash{$newsplit2[$nscount]} . ")";
        } else {
            print "N/A)";
        }

        if (exists $newsplit2[$nscount + 1]) {
        for (my $j=$ltagcolumnhash{$newsplit2[$nscount]}; $j<$ltagcolumnhash{$newsplit2[$nscount + 1]};$j++) {
            print ",";
            $commacount++;
        }
        }
 
$nscount++;
}
# print annotation
for(my $i=$commacount; $i < (scalar(@ARGV)-1); $i++) {
        print ",";
}

foreach my $key (keys %annotationhash) {
    my $thetag = $newsplit2[0];
    my $temp = $annotationhash{$key};
    if ($key =~ m/$thetag/i) {
        my $omitcomma = $annotationhash{$key}; 
        $omitcomma =~ s/,/ /g;
        chomp $omitcomma;
        print $omitcomma;
        last;
    } 
}
                                       
print "\n";  
}

