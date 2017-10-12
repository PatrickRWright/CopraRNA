#!/usr/bin/env perl

use strict;
use warnings;

# prepares an output table to display on the Freiburg RNA tools webserver

my $CopraRNA_out = $ARGV[0];
my $intarna_out = $ARGV[1];

# check CopraRNA1 switch 
my $cop1 = `grep 'CopraRNA1:' CopraRNA_option_file.txt | sed 's/CopraRNA1://g'`; ## edit 2.0.5.1
chomp $cop1;

open(MYDATA, $CopraRNA_out) or die("Error: cannot open file $CopraRNA_out at prepare_output_for_websrv_new.pl\n");
    my @CopraRNA_out_lines = <MYDATA>;
close MYDATA;

open(MYDATA, $intarna_out) or die("Error: cannot open file $intarna_out at prepare_output_for_websrv_new.pl\n");
    my @intarna_out_lines = <MYDATA>;
close MYDATA;

if (scalar(@CopraRNA_out_lines) <= 1) {
    die("\nNo data in $CopraRNA_out at prepare_output_for_websrv_new.pl\n\n");
}

open(WRITEINTERNAL, ">coprarna_internal_table.csv");

my $printedcounter = 1; # when this is on 101 we stop because we have then printed 100 lines

print WRITEINTERNAL "Rank,CopraRNA p-value,CopraRNA fdr,Locus Tag,Gene Name,Energy kcal/mol,IntaRNA p-value,Position mRNA,Position ncRNA,Annotation,Additional homologs,Entrez GeneID,Interaction,Position Seed - mRNA,Position Seed - ncRNA,Hybridization Energy kcal/mol,Unfolding Energy - mRNA kcal/mol,Unfolding Energy - ncRNA kcal/mol,Amount sampled\n"; # header ## edit 2.0.1

for (my $i=1;$i<scalar(@CopraRNA_out_lines);$i++) {
    my @split = split(/,/,$CopraRNA_out_lines[$i]);
    if($split[2]) { ## edit 1.1.0 // check if organism of interest has an entry
        
        # print Rank and CopraRNA p-value and fdr
        print WRITEINTERNAL "$printedcounter,";
        printf WRITEINTERNAL ("%.4g", $split[1]); ## edit 1.1.0
        print WRITEINTERNAL ",";                  ## edit 1.1.0
        printf WRITEINTERNAL ("%.4g", $split[0]); ## edit 1.1.0
        print WRITEINTERNAL ",";
        $printedcounter++;

        # print the IntaRNA outputs for the organism of interest
        # split up organism of interest cell
        my @split_ooi = split(/[\(,|]/, $split[2]); ## edit 2.0.6 changed var name to split_ooi
        # locus tag
        my $ltag = $split_ooi[0]; ## edit 2.0.6
        print WRITEINTERNAL $ltag . ",";
        # gene name
        print WRITEINTERNAL $split_ooi[1] . ",";
        # IntaRNA energy score
        printf WRITEINTERNAL ("%.2f", $split_ooi[2]);
        print WRITEINTERNAL ",";
        # IntaRNA p-value
        printf WRITEINTERNAL ("%.6f", $split_ooi[3]);
        print WRITEINTERNAL ",";
        # interacting region mRNA
        my $intmRNA = $split_ooi[4] . " -- " . $split_ooi[5];
        print WRITEINTERNAL $split_ooi[4] . " -- " . $split_ooi[5] . ",";
        # interacting region ncRNA
        my $intncRNA = $split_ooi[6] . " -- " . $split_ooi[7];
        print WRITEINTERNAL $split_ooi[6] . " -- " . $split_ooi[7] . ",";

        if ($cop1) {
            # annotation
            print WRITEINTERNAL $split[-3] . ","; ## edit 2.0.1
            # additional homologs
            my $temp = $split[-2]; ## edit 2.0.1
            chomp $temp;
            print WRITEINTERNAL $temp . ",";
        } else { ## edit 2.0.6
            # annotation
            print WRITEINTERNAL $split[-2] . ",";
            # additional homologs
            my $temp = $split[-1];
            chomp $temp;
            print WRITEINTERNAL $temp . ",";
        }
        my $GID = "";
        if ($split_ooi[8] =~ m/GeneID:(\d+)\)/) {
            $GID = $1;
        }
        print WRITEINTERNAL $GID . ",";
        my @intarna_array = [];
        for (my $j=0; $j<scalar(@intarna_out_lines);$j++) {
            if ($intarna_out_lines[$j] =~ m/$ltag/i) {
                for (my $line=($j+5);$line<=($j+18);$line++) {
                    push(@intarna_array, $intarna_out_lines[$line]);
                }
                last;
            }
        }
        chomp $intarna_array[1];
        chomp $intarna_array[2];
        chomp $intarna_array[3];
        chomp $intarna_array[4];
        chomp $intarna_array[7];
        $intarna_array[7] =~ s/.*: //;
        chomp $intarna_array[10];
        $intarna_array[10] =~ s/.*: //;
        chomp $intarna_array[12];
        $intarna_array[12] =~ s/.*: //;
        chomp $intarna_array[13];
        $intarna_array[13] =~ s/.*: //;
        chomp $intarna_array[14];
        $intarna_array[14] =~ s/.*: //;
        # this is the interaction
        #print WRITEINTERNAL $intarna_array[1] . "\\n" . $intarna_array[2] . "\\n" . $intarna_array[3] . "\\n" . $intarna_array[4] .",";       
 
        # call the subroutine here
        print WRITEINTERNAL &prune_interaction($intarna_array[1],$intarna_array[2],$intarna_array[3],$intarna_array[4],$intmRNA,$intncRNA); 
        print WRITEINTERNAL ",";       
 
        # these are the interaction properties
        if ($cop1) {
            print WRITEINTERNAL $intarna_array[7] . "," . $intarna_array[10] . "," . $intarna_array[14] . "," . $intarna_array[12] . "," . $intarna_array[13] . "," . $split[-1]; ## edit 2.0.1
        } else { ## edit 2.0.6
            print WRITEINTERNAL $intarna_array[7] . "," . $intarna_array[10] . "," . $intarna_array[14] . "," . $intarna_array[12] . "," . $intarna_array[13] . "," . "0\n";
        }       
    }
}

close(WRITEINTERNAL);




#&prune_interaction("5'-AAGUAUUUUUCAGCUUUUCAUUCUGACUGCAACGGGCAAUAUGUCUCUGUGUGGAUUAAAAAAAGAGUGUCUGAUAGCAGCUUCUGAACUGGUU             U AAUUAAAAUUUUAUU    UAGG                 AACCAAUAUAGGCAUAGCGCACAGACAGAUAAAAAUUACAGAGUACACAACAUCCAUGAAACGCAUUAGCACCACCAUUACCACCACCAUCACCAUUACCACAGGUAACGGUGCGGGCUGACGCGUACAGGAAACACAGAAAAAAGCCCGCACCU-3'\n","                                                                                                 ACCUGC CGUGAG A               GACU    UCACUAA   AUACUUU\n","                                                                                                 UGGACG GUACUC U               CUGA    AGUGGUU   UAUGGAG\n","                                                                                       3'-UUUUUUU      C                                      AGU       ACCC-5'\n", "101 -- 222", "74 -- 92");


#&prune_interaction("5'-CAGCUUCUGAACUGGUU             U AAUUAAAAUUUUAUU    UAGG                 AACCAAUAACAGAAAAAAGCCCGCACCU-3'\n","ACCUGC CGUGAG A               GACU    UCACUAA   AUACUUU\n","UGGACG GUACUC U               CUGA    AGUGGUU   UAUGGAG\n","3'-AAAAAAAAUUUUUUU      C                                      AGU       ACCCAAAAAAAA-5'\n", "1 -- 12", "16 -- 20");

#&prune_interaction("5'-U             U AAUUAAAAUUUUAUU    UAGG                 CGCU-3'\n","ACCUGC CGUGAG A               GACU    UCACUAA   AUACUUU\n","UGGACG GUACUC U               CUGA    AGUGGUU   UAUGGAG\n","3'-AUUUU      C                                      AGU       -5'\n", "77 -- 98", "1000 -- 2000");
#print "##########################\n";
#&prune_interaction("5'-AA          A-3'\n","CUCUCUCUCU\n","GGGGGGGGGG\n","3'-          GG-5'\n", "3 -- 12", "3 -- 12");
#print "##########################\n";
#&prune_interaction("5'-AA          -3'\n","CUCUCUCUCU\n","GGGGGGGGGG\n","3'-          -5'\n", "3 -- 12", "1 -- 10");

#&prune_interaction("5'-AAAAA          GGGG-3'\n","CUCUCUCUCU\n","GGGGGGGGGG\n","3'-GUCG          U-5'\n", "6 -- 15", "2 -- 11");

#&prune_interaction("5'-          -3'\n","CUCUCUCUCU\n","GGGGGGGGGG\n","3'-          -5'\n", "1 -- 10", "1 -- 10");

#&prune_interaction("5'-          AAGGG-3'\n","CUCUCUCUCU\n","GGGGGGGGGG\n","3'-AAAA          -5'\n", "1 -- 10", "1 -- 10");

sub prune_interaction
{
    ###
    ### instead of printing, append all final parts to one string which is then returned at the end of the sub
    ###
    
    my $line1 = $_[0];
    my $line2 = $_[1];
    my $line3 = $_[2];
    my $line4 = $_[3];
    my $regionTarget = $_[4];
    my $regionNCRNA = $_[5];

    $line1 =~ s/^\s+//;
    $line1 =~ s/\s+$//;
    $line4 =~ s/^\s+//;
    $line4 =~ s/\s+$//;
    $line2 =~ s/^\s+//;
    $line2 =~ s/\s+$//;
    $line3 =~ s/^\s+//;
    $line3 =~ s/\s+$//;   
    my $finalOut = "";

    my $flankingswitch5primeTop = 0;    
    my $flankingswitch3primeTop = 0;    

    # print position
    my @splitposiTarget = split(/\s/,$regionTarget);
    my @splitposiNCRNA = split(/\s/, $regionNCRNA);

 
    my @split5to3prime = split(/\s/, $line1);
    my $fivePrimeTop = shift @split5to3prime;
    my $threePrimeTop = pop @split5to3prime;

    my $posTargetUp = $splitposiTarget[0] - 1;
    ########disallow exceeding the boundaries
    $posTargetUp = $posTargetUp + 1 if(length($fivePrimeTop) == 3);
    
    if(length($fivePrimeTop) == 3) {
        $finalOut = $finalOut . "              " . $posTargetUp;
    } else {    
        $finalOut = $finalOut . "             " . $posTargetUp;
    }
    
    if(length($fivePrimeTop) == 3) { $flankingswitch5primeTop = 1; }
    if(length($threePrimeTop) == 3) { $flankingswitch3primeTop = 1; }
    my $sumswitchTop = $flankingswitch5primeTop + $flankingswitch3primeTop;
 
    for (my $spaces=1;$spaces<=(1 + length($line2) - length($posTargetUp) - $sumswitchTop);$spaces++) {
        $finalOut = $finalOut . " ";
    }


    my $posTargetDown = $splitposiTarget[-1] + 1;
    ########disallow exceeding the boundaries
    $posTargetDown = $posTargetDown - 1 if(length($threePrimeTop) == 3);
    $finalOut = $finalOut . $posTargetDown . "\\n";

    if(length($fivePrimeTop) == 3) {
        $finalOut = $finalOut . "              |";
    } else {
        $finalOut = $finalOut . "             |";
    }

    for (my $spaces=1;$spaces<=length($line2) - $sumswitchTop;$spaces++) {
        $finalOut = $finalOut . " ";
    }

    $finalOut = $finalOut . "|\\n";

    # make first line of interaction
    @split5to3prime = split(/\s/, $line1);
    $fivePrimeTop = shift @split5to3prime;
    $threePrimeTop = pop @split5to3prime;
    # remove the flanks
    $line1 =~ s/$fivePrimeTop//;
    $line1 =~ s/$threePrimeTop//;
    # specify top 5' end
    my $topleft = "";
    if(length($fivePrimeTop) > 11) {
        $topleft = substr($fivePrimeTop,0,6);
        $fivePrimeTop = reverse($fivePrimeTop);
        $topleft = $topleft . "..." . reverse(substr($fivePrimeTop,0,5));
        $finalOut = $finalOut . $topleft;
    } else {
        $topleft = $fivePrimeTop;
        my $temp = length($topleft);
        for (my $spaces=1;$spaces<=(14-$temp);$spaces++) { $topleft = " " . $topleft; }
        $finalOut = $finalOut . $topleft;
    }

    chomp $line1; 
    $finalOut = $finalOut . $line1;
    
    my $topright = "";
    if(length($threePrimeTop) > 11) {
        $topright = substr($threePrimeTop,0,5);
        $threePrimeTop = reverse($threePrimeTop);
        $topright = $topright . "..." . reverse(substr($threePrimeTop,0,6));
        $finalOut = $finalOut . $topright;
    } else {
        $topright = $threePrimeTop;
        $finalOut = $finalOut . $topright;
    }

    $finalOut = $finalOut . "\\n";    

    # make second and third line of interaction

    for (my $i=1;$i<=length($topleft);$i++) {
        $finalOut = $finalOut . " ";
    }
    $finalOut = $finalOut . $line2 . "\\n";

    # add a line with | and : for the interacting nucleotides here
    for (my $i=1;$i<=length($topleft);$i++) {
        $finalOut = $finalOut . " ";
    }

    for (my $j=0;$j<length($line2);$j++) {
        if(substr($line2,$j,1) eq "A" and substr($line3,$j,1) eq "U") {
            $finalOut = $finalOut . "|";
        }
        elsif(substr($line2,$j,1) eq "U" and substr($line3,$j,1) eq "A") {
            $finalOut = $finalOut . "|";
        }
        elsif(substr($line2,$j,1) eq "G" and substr($line3,$j,1) eq "C") {
            $finalOut = $finalOut . "|";
        } 
        elsif(substr($line2,$j,1) eq "C" and substr($line3,$j,1) eq "G") {
            $finalOut = $finalOut . "|";
        }
        elsif(substr($line2,$j,1) eq "G" and substr($line3,$j,1) eq "U") {
            $finalOut = $finalOut . ":";
        }
        elsif(substr($line2,$j,1) eq "U" and substr($line3,$j,1) eq "G") {
            $finalOut = $finalOut . ":";
        } else {
            $finalOut = $finalOut . " ";
        }
    }
    $finalOut = $finalOut . "\\n";


    for (my $i=1;$i<=length($topleft);$i++) {
        $finalOut = $finalOut . " ";
    }
    $finalOut = $finalOut . $line3 . "\\n";

    # make fourth line of interaction
    my @split3to5prime = split(/\s/,$line4);
    my $fivePrimeTail = shift @split3to5prime;
    my $threePrimeTail = pop @split3to5prime;
    # remove flanks
    $line4 =~ s/$fivePrimeTail//;
    $line4 =~ s/$threePrimeTail//;
    my $bottomLeft = "";
    
    if(length($fivePrimeTail) > 11) {
        $bottomLeft = substr($fivePrimeTail,0,6);
        $fivePrimeTail = reverse($fivePrimeTail);
        $bottomLeft = $bottomLeft . "..." . reverse(substr($fivePrimeTail,0,5));
        $finalOut = $finalOut . $bottomLeft;
    } else {
        $bottomLeft = $fivePrimeTail;
        my $temp = length($bottomLeft);
        for (my $spaces=1;$spaces<=(14-$temp);$spaces++) { $bottomLeft = " " . $bottomLeft; }
        $finalOut = $finalOut . $bottomLeft;
    }

    chomp $line4;
    $finalOut = $finalOut . $line4;

    my $bottomRight = "";
    if(length($threePrimeTail) > 11) {
        $bottomRight = substr($threePrimeTail,0,5);
        $threePrimeTail = reverse($threePrimeTail);
        $bottomRight = $bottomRight . "..." . reverse(substr($threePrimeTail,0,6));
        $finalOut = $finalOut . $bottomRight;
    } else {
        $bottomRight = $threePrimeTail;
        $finalOut = $finalOut . $bottomRight;
    }
    $finalOut = $finalOut . "\\n";
 
    # switches
    my $flankingswitch5primeBottom = 0;
    my $flankingswitch3primeBottom = 0;
    my $sumswitchBottom = 0;       

    if(length($threePrimeTail) == 3) { $flankingswitch5primeBottom = 1; }
    if(length($fivePrimeTail) == 3) { $flankingswitch3primeBottom = 1; }
    $sumswitchBottom = $flankingswitch5primeBottom + $flankingswitch3primeBottom;

    if(length($fivePrimeTail) == 3) {
        $finalOut = $finalOut . "              |";
    } else {
        $finalOut = $finalOut . "             |";
    }


    for (my $spaces=1;$spaces<=length($line3) - $sumswitchBottom;$spaces++) {
        $finalOut = $finalOut . " ";
    }

    $finalOut = $finalOut . "|\\n";

    my $posTargetUp2 = $splitposiNCRNA[0] - 1;
    my $posTargetDown2 = $splitposiNCRNA[-1] + 1;

    ############ not letting indices out of bounds
    $posTargetUp2 = $posTargetUp2 + 1 if(length($threePrimeTail) == 3);
    $posTargetDown2 = $posTargetDown2 - 1 if(length($fivePrimeTail) == 3);

    if(length($fivePrimeTail) == 3) {
        $finalOut = $finalOut . "              " . $posTargetDown2;
    } else {
        $finalOut = $finalOut . "             " . $posTargetDown2;
    }

    for (my $spaces=1;$spaces<=(1 + length($line3) - length($posTargetDown2) - $sumswitchBottom);$spaces++) {
        $finalOut = $finalOut . " ";
    }
    $finalOut = $finalOut . $posTargetUp2 . "\\n"; 
    

    # spaces -> "_"    
    $finalOut =~ s/\s/_/g;
    return $finalOut;

}

