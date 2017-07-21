#!/usr/bin/env perl

use strict;
use warnings;

use Bio::SeqIO;

my $gbk = $ARGV[0];
my $up = $ARGV[1];
my $down = $ARGV[2];
my $region = $ARGV[3];

unless(defined $ARGV[0] and $ARGV[1] and $ARGV[2] and $ARGV[3]) { die("usage: ./parse_region_from_genome.pl N[ZC]_XXXXXXXX ntup(number) ntdown(number) region(i.e. 5utr/3utr/cds) > outfile.fas\n"); }

# comment: with up 200 donw 100 the startcodon is on position 201 202 203 of the 300 nt long string

my $in  = Bio::SeqIO->new(-file => $gbk , '-format' => 'genbank');
my $seq = $in->next_seq();

my $in2  = Bio::SeqIO->new(-file => $gbk , '-format' => 'genbank');

if ($region eq "5utr") {
    while ( my $seq2 = $in2->next_seq() ) {
        foreach my $sf2 ( $seq2->get_SeqFeatures() ) {
            if( $sf2->primary_tag eq 'CDS' ) {
            my $ltag = "";
            if ($sf2->has_tag("locus_tag")) {
                my @ltaglist = $sf2->get_tag_values("locus_tag");
                $ltag = $ltaglist[0];
            }
            my $start = $sf2->start;
            my $end = $sf2->end;
            my $loc = $sf2->location->strand;
            if($loc == 1) {
                if ( $sf2->location->isa('Bio::Location::SplitLocationI') )  {
                    my @locarray = ();
                    foreach my $location ( $sf2->location->sub_Location ) {
                        push(@locarray,$location->start);
                    }
                    print &get_region(1, $locarray[0], $up, $down, $ltag, $seq);
                } else {
                    print &get_region(1, $start, $up, $down, $ltag, $seq);
                } 
            } else {
                if ( $sf2->location->isa('Bio::Location::SplitLocationI') )  {
                    my @locarray = ();
                    foreach my $location ( $sf2->location->sub_Location ) {
                        push(@locarray,$location->end); 
                    }
                    print &get_region(-1, $locarray[-1], $up, $down, $ltag, $seq);
                } else {
                    print &get_region(-1, $end, $up, $down, $ltag, $seq);
                }
            }
            }
        }
    }
}

if ($region eq "3utr") {
     while ( my $seq2 = $in2->next_seq() ) {
        foreach my $sf2 ( $seq2->get_SeqFeatures() ) {
            if( $sf2->primary_tag eq 'CDS' ) {
            my $ltag = "";
            if ($sf2->has_tag("locus_tag")) {
                my @ltaglist = $sf2->get_tag_values("locus_tag");
                $ltag = $ltaglist[0];
            }
            my $start = $sf2->start;
            my $end = $sf2->end;
            my $loc = $sf2->location->strand;
            if($loc == 1) {
                if ( $sf2->location->isa('Bio::Location::SplitLocationI') )  {
                    my @locarray = ();
                    foreach my $location ( $sf2->location->sub_Location ) {
                        push(@locarray,$location->end); 
                    }
                    print &get_region(1, ($locarray[-1] + 1), $up, $down, $ltag, $seq);
                } else {
                    print &get_region(1, ($end + 1), $up, $down, $ltag, $seq);
                }
            } else {
                if ( $sf2->location->isa('Bio::Location::SplitLocationI') )  {
                    my @locarray = ();
                    foreach my $location ( $sf2->location->sub_Location ) {
                        push(@locarray,$location->start);
                    }
                    print &get_region(-1, ($locarray[0] - 1), $up, $down, $ltag, $seq);
                } else {
                    print &get_region(-1, ($start - 1), $up, $down, $ltag, $seq);
                }
            }
            }
        }
    }
}

if ($region eq "cds") {
     while ( my $seq2 = $in2->next_seq() ) {
        foreach my $sf2 ( $seq2->get_SeqFeatures() ) {
            if( $sf2->primary_tag eq 'CDS' ) {
            my $ltag = "";
            if ($sf2->has_tag("locus_tag")) {
                my @ltaglist = $sf2->get_tag_values("locus_tag");
                $ltag = $ltaglist[0];
            }
            my $start = $sf2->start;
            my $end = $sf2->end;
            my $loc = $sf2->location->strand;
            my $tempStart = $start;
            my $tempUp = $up;
            my $tempDown = $down;
            if($loc == 1) {
                if ( $sf2->location->isa('Bio::Location::SplitLocationI') )  {
                    my @locarray = (); # ends
                    my @locarray2 = (); # starts
                    foreach my $location ( $sf2->location->sub_Location ) {
                        push(@locarray,$location->end); 
                        push(@locarray2,$location->start); 
                    }
                    my $printseq1 = &get_region_cds(1, $locarray2[0], $up, 0, $ltag, $seq);
                    chomp $printseq1;
                    my $printseq2 = $sf2->spliced_seq->seq;
                    my $printseq3 = &get_region_cds(1, ($locarray[-1] + 1), 0, $down, $ltag, $seq);
                    $printseq3 =~ s/>.*\n//; # remove locus_tag
                    my $printseq = $printseq1 . $printseq2 . $printseq3;
                    print $printseq;
                } else {
                    if ($tempStart==1) {
                        $start++;
                        $tempUp++;
                    }
                    my $printseq1 = &get_region_cds(1, $start, $tempUp, 0, $ltag, $seq);
                    chomp $printseq1;
                    
                    chop $printseq1 if($tempStart == 1);

                    my $printseq2 = $sf2->spliced_seq->seq;
                    my $printseq3 = &get_region_cds(1, $end, 0, ($down + 1), $ltag, $seq);
                    $printseq3 =~ s/>.*\n//; # remove locus_tag
                    $printseq3 = reverse($printseq3);
                    chop $printseq3;
                    $printseq3 = reverse($printseq3);
                    my $printseq = $printseq1 . $printseq2 . $printseq3;
                    print $printseq;
                }
            } else {
                if ( $sf2->location->isa('Bio::Location::SplitLocationI') )  {
                    my @locarray = (); # starts
                    my @locarray2 = (); # ends
                    foreach my $location ( $sf2->location->sub_Location ) {
                        push(@locarray,$location->start); 
                        push(@locarray2,$location->end); 
                    }
                    my $printseq1 = &get_region_cds(-1, $locarray2[-1], $up, 0, $ltag, $seq);
                    chomp $printseq1;
                    my $printseq2 = $sf2->spliced_seq->seq;
                    my $printseq3 = &get_region_cds(-1, ($locarray[0] - 1), 0, $down, $ltag, $seq);
                    $printseq3 =~ s/>.*\n//; # remove locus_tag
                    my $printseq = $printseq1 . $printseq2 . $printseq3;
                    print $printseq;
                } else {
                    if ($tempStart==1) {
                        $start++;
                        $tempDown++;
                    }
                    my $printseq1 = &get_region_cds(-1, ($end - 1), ($up + 1), 0, $ltag, $seq);
                    chomp $printseq1;
                    chop $printseq1;
                    my $printseq2 = $sf2->spliced_seq->seq;
                    my $printseq3 = &get_region_cds(-1, ($start - 1), 0, $tempDown, $ltag, $seq);
                    $printseq3 =~ s/>.*\n//; # remove locus_tag
                    if ($tempStart==1) {
                        $printseq3 = reverse($printseq3);
                        chop $printseq3;
                        $printseq3 = reverse($printseq3);
                    }
                    my $printseq = $printseq1 . $printseq2 . $printseq3;
                    print $printseq;
                }
            }
            }
        }
    }

}


# subs

sub get_region
{
        my $strand = $_[0];
        my $pos = $_[1];
        my $upstream = $_[2];
        my $downstream = $_[3];
        my $ltag = $_[4];
        my $seq_obj = $_[5];
        my $genomelength = $seq_obj->length;
        my $circularTest = $seq_obj->is_circular;

        if ($strand == -1) {
            if ((($pos - $downstream) >= 1) and (($pos + $upstream) <= $genomelength)) { 
                my $printseq = reverse($seq_obj->subseq(($pos - $downstream + 1), ($pos + $upstream)));
                $printseq =~ tr/[g,a,t,c,G,A,T,C]/[c,t,a,g,C,T,A,G]/; 
                return ">$ltag\n$printseq\n";
            } else { # this exception deals with sequences at the sides of the genome
                unless($circularTest) { # only do this if genome is not curcular
                    if (($pos - $downstream) < 1) {
                        my $printseq = reverse($seq_obj->subseq(1, ($pos + $upstream)));
                        $printseq =~ tr/[g,a,t,c,G,A,T,C]/[c,t,a,g,C,T,A,G]/;
                        if (length($printseq) >= 150) { ## edit 2.0.2, changed this 140 -> 150
                            return ">$ltag\n$printseq\n";
                        } else { return; }
                    }
                    if (($pos + $upstream) > $genomelength) {
                        my $printseq = reverse($seq_obj->subseq(($pos - $downstream + 1), $genomelength));
                        $printseq =~ tr/[g,a,t,c,G,A,T,C]/[c,t,a,g,C,T,A,G]/;
                        if (length($printseq) >= 150) { ## edit 2.0.2, changed this 140 -> 150
                            return ">$ltag\n$printseq\n";
                        } else { return; }
                    }
                } else { # here we add the exception for the edges in circular genomes 
                    if (($pos - $downstream) < 1) {
                        my $overhang = ($pos - $downstream) * (-1);
                        my $printseq1 = reverse($seq_obj->subseq(1, ($pos + $upstream)));
                        my $printseq2 = reverse($seq_obj->subseq(($genomelength - $overhang + 1), $genomelength));
                        my $printseq = $printseq1 . $printseq2;
                        $printseq =~ tr/[g,a,t,c,G,A,T,C]/[c,t,a,g,C,T,A,G]/;
                        return ">$ltag\n$printseq\n";
                    }
                    if (($pos + $upstream) > $genomelength) {
                        my $overhang = ($pos + $upstream) - $genomelength;
                        my $printseq1 = reverse($seq_obj->subseq(1, $overhang));
                        my $printseq2 = reverse($seq_obj->subseq(($pos - $downstream + 1), $genomelength));
                        my $printseq = $printseq1 . $printseq2;
                        $printseq =~ tr/[g,a,t,c,G,A,T,C]/[c,t,a,g,C,T,A,G]/;
                        return ">$ltag\n$printseq\n";
                    }
                }
            }
        } else {
            if ((($pos - $upstream) >= 1) and (($pos + $downstream) <= $genomelength)) {
                my $printseq = $seq_obj->subseq(($pos - $upstream), ($pos + $downstream - 1));
                return ">$ltag\n$printseq\n";
            } else { # this exception deals with sequences at the sides of the genome
                unless($circularTest) { # only do this if genome is not curcular
                     if (($pos - $upstream) < 1) {
                         my $printseq = $seq_obj->subseq(1, ($pos + $downstream - 1));
                         if (length($printseq) >= 150) { ## edit 2.0.2, changed this 140 -> 150
                             return ">$ltag\n$printseq\n";
                         } else { return; }
                     }
                     if (($pos + $downstream) > $genomelength) {
                         my $printseq = $seq_obj->subseq(($pos - $upstream), $genomelength);
                         if (length($printseq) >= 150) { ## edit 2.0.2, changed this 140 -> 150
                             return ">$ltag\n$printseq\n";
                         } else { return; }
                     }
                } else { # here we add the exception for the edges in circular genomes
                     if (($pos - $upstream) < 1) {
                          my $overhang = ($pos - $upstream) * (-1);
                          my $printseq2 = $seq_obj->subseq(1, ($pos + $downstream - 1));
                          my $printseq1 = $seq_obj->subseq(($genomelength - $overhang), $genomelength);
                          my $printseq = $printseq1 . $printseq2;
                          return ">$ltag\n$printseq\n";
                     }
                     if (($pos + $downstream) > $genomelength) {
                          my $overhang = ($pos + $downstream) - $genomelength;
                          my $printseq1 = $seq_obj->subseq(($pos - $upstream), $genomelength);
                          my $printseq2 = $seq_obj->subseq(1, ($overhang - 1));
                          my $printseq = $printseq1 . $printseq2;
                          return ">$ltag\n$printseq\n";
                     }
                }
            }
        }
}

sub get_region_cds
{
        my $strand = $_[0];
        my $pos = $_[1];
        my $upstream = $_[2];
        my $downstream = $_[3];
        my $ltag = $_[4];
        my $seq_obj = $_[5];
        my $genomelength = $seq_obj->length;
        my $circularTest = $seq_obj->is_circular;

        if ($strand == -1) {
            if ((($pos - $downstream) >= 1) and (($pos + $upstream) <= $genomelength)) { 
                my $printseq = reverse($seq_obj->subseq(($pos - $downstream + 1), ($pos + $upstream)));
                $printseq =~ tr/[g,a,t,c,G,A,T,C]/[c,t,a,g,C,T,A,G]/; 
                return ">$ltag\n$printseq\n";
            } else { # this exception deals with sequences at the sides of the genome
                unless($circularTest) { # only do this if genome is not curcular
                    if (($pos - $downstream) < 1) {
                        my $printseq = reverse($seq_obj->subseq(1, ($pos + $upstream)));
                        $printseq =~ tr/[g,a,t,c,G,A,T,C]/[c,t,a,g,C,T,A,G]/;
                            return ">$ltag\n$printseq\n";
                    }
                    if (($pos + $upstream) > $genomelength) {
                        my $printseq = reverse($seq_obj->subseq(($pos - $downstream + 1), $genomelength));
                        $printseq =~ tr/[g,a,t,c,G,A,T,C]/[c,t,a,g,C,T,A,G]/;
                            return ">$ltag\n$printseq\n";
                    }
                } else { # here we add the exception for the edges in circular genomes 
                    if (($pos - $downstream) < 1) {
                        my $overhang = ($pos - $downstream) * (-1);
                        my $printseq1 = reverse($seq_obj->subseq(1, ($pos + $upstream)));
                        my $printseq2 = reverse($seq_obj->subseq(($genomelength - $overhang + 1), $genomelength));
                        my $printseq = $printseq1 . $printseq2;
                        $printseq =~ tr/[g,a,t,c,G,A,T,C]/[c,t,a,g,C,T,A,G]/;
                        return ">$ltag\n$printseq\n";
                    }
                    if (($pos + $upstream) > $genomelength) {
                        my $overhang = ($pos + $upstream) - $genomelength;
                        my $printseq1 = reverse($seq_obj->subseq(1, $overhang));
                        my $printseq2 = reverse($seq_obj->subseq(($pos - $downstream + 1), $genomelength));
                        my $printseq = $printseq1 . $printseq2;
                        $printseq =~ tr/[g,a,t,c,G,A,T,C]/[c,t,a,g,C,T,A,G]/;
                        return ">$ltag\n$printseq\n";
                    }
                }
            }
        } else {
            if ((($pos - $upstream) >= 1) and (($pos + $downstream) <= $genomelength)) {
                my $printseq = $seq_obj->subseq(($pos - $upstream), ($pos + $downstream - 1));
                return ">$ltag\n$printseq\n";
            } else { # this exception deals with sequences at the sides of the genome
                unless($circularTest) { # only do this if genome is not curcular
                     if (($pos - $upstream) < 1) {
                         my $printseq = $seq_obj->subseq(1, ($pos + $downstream - 1));
                             return ">$ltag\n$printseq\n";
                     }
                     if (($pos + $downstream) > $genomelength) {
                         my $printseq = $seq_obj->subseq(($pos - $upstream), $genomelength);
                             return ">$ltag\n$printseq\n";
                     }
                } else { # here we add the exception for the edges in circular genomes
                     if (($pos - $upstream) < 1) {
                          my $overhang = ($pos - $upstream) * (-1);
                          my $printseq2 = $seq_obj->subseq(1, ($pos + $downstream - 1));
                          my $printseq1 = $seq_obj->subseq(($genomelength - $overhang), $genomelength);
                          my $printseq = $printseq1 . $printseq2;
                          return ">$ltag\n$printseq\n";
                     }
                     if (($pos + $downstream) > $genomelength) {
                          my $overhang = ($pos + $downstream) - $genomelength;
                          my $printseq1 = $seq_obj->subseq(($pos - $upstream), $genomelength);
                          my $printseq2 = $seq_obj->subseq(1, ($overhang - 1));
                          my $printseq = $printseq1 . $printseq2;
                          return ">$ltag\n$printseq\n";
                     }
                }
            }
        }
}




