#!/usr/bin/env perl

use strict;
use warnings;
use Cwd 'abs_path'; ## edit 2.0.5.1

# get absolute path
my $ABS_PATH = abs_path($0); ## edit 2.0.5.1
# remove script name at the end
# match all non slash characters at the end of the string
$ABS_PATH =~ s|[^/]+$||g; ## edit 2.0.5.1
my $PATH_COPRA_SUBSCRIPTS = $ABS_PATH;

my $ncrnas = $ARGV[0]; # input_sRNA.fa
my $upfromstartpos = $ARGV[1]; # 200
my $down = $ARGV[2]; # 100
my $mrnapart = $ARGV[3]; # cds or 5utr or 3utr
my $GenBankFiles = "";
my $orgcount = 0;

my $cores = `grep 'core count:' CopraRNA_option_file.txt | grep -oP '\\d+'`; ## edit 2.0.4
chomp $cores; ## edit 2.0.4

# check if CopraRNA1 prediction should be made
my $cop1 = `grep 'CopraRNA1:' CopraRNA_option_file.txt | sed 's/CopraRNA1://g'`; ## edit 2.0.5.1
chomp $cop1;

# check nooi switch
my $nooi = `grep 'nooi:' CopraRNA_option_file.txt | sed 's/nooi://g'`; ## edit 2.0.6
chomp $nooi;

# check for verbose printing
my $verbose = `grep 'verbose:' CopraRNA_option_file.txt | sed 's/verbose://g'`; ## edit 2.0.5.1
chomp $verbose;

# get amount of top predictions to return
my $topcount = `grep 'top count:' CopraRNA_option_file.txt | grep -oP '\\d+'`; ## edit 2.0.5.1
chomp $topcount;
$topcount++; # need this to include the header

# check for websrv output printing
my $websrv = `grep 'websrv:' CopraRNA_option_file.txt | sed 's/websrv://g'`; ## edit 2.0.5.1
chomp $websrv;

# check for enrichment on/off and count
my $enrich = `grep 'enrich:' CopraRNA_option_file.txt | sed 's/enrich://g'`; ## edit 2.0.5.1
chomp $enrich;

# get window size option
my $winsize = `grep 'win size:' CopraRNA_option_file.txt | sed 's/win size://g'`; ## edit 2.0.5.1
chomp $winsize;

# get maximum base pair distance 
my $maxbpdist = `grep 'max bp dist:' CopraRNA_option_file.txt | sed 's/max bp dist://g'`; ## edit 2.0.5.1
chomp $maxbpdist;

# get consensus prediction option
my $cons = `grep 'cons:' CopraRNA_option_file.txt | sed 's/cons://g'`; ## edit 2.0.6
chomp $cons;

# get ooifilt
my $ooi_filt = `grep 'ooifilt:' CopraRNA_option_file.txt | sed 's/ooifilt://g'`;
chomp $ooi_filt;

open ERRORLOG, ">>err.log" or die("\nError: cannot open file err.log in homology_intaRNA.pl\n\n"); ## edit 2.0.2 

my $keggtorefseqnewfile = $PATH_COPRA_SUBSCRIPTS . "kegg2refseqnew.csv";
# RefSeqID -> space separated RefSeqIDs // 'NC_005140' -> 'NC_005139 NC_005140 NC_005128'
my %refseqaffiliations = ();

# read kegg2refseqnew.csv
open(MYDATA, $keggtorefseqnewfile) or die("\nError: cannot open file $keggtorefseqnewfile in homology_intaRNA.pl\n\n");
    my @keggtorefseqnew = <MYDATA>;
close MYDATA;

# get the refseq affiliations
foreach(@keggtorefseqnew) {
    # split off quadruplecode (pseudokegg id)
    my @split = split("\t", $_);
    my $all_refseqs = $split[1];
    chomp $all_refseqs;
    # split up refseq ids
    my @split_refseqs = split(/\s/, $all_refseqs);
    foreach(@split_refseqs) {
        $refseqaffiliations{$_} = $all_refseqs;
    }
}

# add "ncRNA_" to fasta headers
system "sed 's/>/>ncRNA_/g' $ncrnas > ncrna.fa"; ## edit 2.0.5.1 // replaced put_ncRNA_fasta_together.pl with this statement

# assign correct refseq IDs for each sequence
for (my $i=4;$i<scalar(@ARGV);$i++) {
    # split up the refseq list for one organism
    my @split = split(/\s/, $refseqaffiliations{$ARGV[$i]});
    # get the first id entry
    my $first_refseq_id = $split[0];
    # override in ncrna.fa
    system "sed -i 's/$ARGV[$i]/$first_refseq_id/g' ncrna.fa";
}

# override $ncrnas variable
$ncrnas = "ncrna.fa";

# get Orgcount
$orgcount = (scalar(@ARGV) - 4);

## prepare input for combine_clusters.pl
## Download Refseq files by Refseq ID 
my $RefSeqIDs = `grep ">" input_sRNA.fa | tr '\n' ' ' | sed 's/>//g'`; ## edit 2.0.5.1
my @split_RefIds = split(/\s+/, $RefSeqIDs);

foreach(@split_RefIds) {
    my $currRefSeqID = $_;

    my $presplitreplicons = $refseqaffiliations{$currRefSeqID};
    my @replikons = split(/\s/, $presplitreplicons); # added this
    
    foreach(@replikons) {
        my $refseqoutputfile = $_ . ".gb"; # added .gb
        $GenBankFiles = $GenBankFiles . $refseqoutputfile . ",";
        my $accessionnumber = $_;
        print $PATH_COPRA_SUBSCRIPTS  . "get_refseq_from_refid.pl -acc $accessionnumber -g $accessionnumber.gb \n" if ($verbose); ## edit 1.2.1 ## edit 2.0.2
        system $PATH_COPRA_SUBSCRIPTS . "get_refseq_from_refid.pl -acc $accessionnumber -g $accessionnumber.gb"; ## edit 1.2.1 ## edit 2.0.2
    }
    chop $GenBankFiles;
    $GenBankFiles = $GenBankFiles . " ";
}

## RefSeq correct download check for 2nd try ## edit 1.2.5
my @files = ();
@files = <*gb>;

foreach(@files) {
    open(GBDATA, $_) or die("\nError: cannot open file $_ in homology_intaRNA.pl\n\n");
        my @gblines = <GBDATA>;
    close GBDATA;

    my $lastLine = $gblines[-2]; ## edit 2.0.2
    my $lastLine_new = $gblines[-1]; ## edit 2.0.3, because of new file donwload the bottom differs
    if ($lastLine =~ m/^\/\//) {
        # all is good
    } elsif ($lastLine_new =~ m/^\/\//) {
        # all is good
    } else {
        system "rm $_"; # remove file to try download again later
    }
}

## refseq availability check
@files = ();

my @totalrefseqFiles = split(/\s|,/, $GenBankFiles);
my $consistencyswitch = 1;

my $limitloops = 0;

my $sleeptimer = 30; ## edit 1.2.0
while($consistencyswitch) {
    @files = ();
    @files = <*gb>;
    foreach(@totalrefseqFiles) {
        chomp $_;
        my $value = $_;
        if(grep( /^$value$/, @files )) { 
            $consistencyswitch = 0;
        } else {
             $limitloops++;
             $consistencyswitch = 1;
 
             if($limitloops > 100) { 
                 $consistencyswitch = 0;
                 print ERRORLOG "Not all RefSeq *gb files downloaded correctly. Restart your job.\n"; 
                 last;
             }
             my $accNr = $_;
             chop $accNr;
             chop $accNr;
             chop $accNr;
             sleep $sleeptimer; ## edit 1.2.0
             $sleeptimer = $sleeptimer * 1.1; ## edit 1.2.0
             print "next try: " . $PATH_COPRA_SUBSCRIPTS . "get_refseq_from_refid.pl -acc $accNr -g $accNr.gb\n" if ($verbose); ## edit 1.2.0 ## edit 1.2.1 ## edit 2.0.2
             system $PATH_COPRA_SUBSCRIPTS . "get_refseq_from_refid.pl -acc $accNr -g $accNr.gb"; ## edit 1.2.1 ## edit 2.0.2
             last;
        }
    }
}

### end availability check


### refseq correct DL check kill job ## edit 1.2.5
@files = <*gb>;

foreach(@files) {
    open(GBDATA, $_) or die("\nError: cannot open file $_ in homology_intaRNA.pl\n\n");
        my @gblines = <GBDATA>;
    close GBDATA;

    my $lastLine = $gblines[-2]; ## edit 2.0.2
    my $lastLine_new = $gblines[-1]; ## edit 2.0.3, because of new file donwload the bottom differs
    if ($lastLine =~ m/^\/\//) {
        # all is good
    } elsif ($lastLine_new =~ m/^\/\//) {
        # all is good
    } else {
        print ERRORLOG "File $_ did not download correctly. This is probably due to a connectivity issue on your or the NCBI's side. Please try to resubmit your job later (~2h.).\n"; # kill ## edit 2.0.2
    }
}


## fixing issue with CONTIG and ORIGIN both in gbk file (can't parse without this) ## edit 1.2.4

@files = <*gb>;

foreach (@files) {
    system "sed -i '/^CONTIG/d' $_"; ## d stands for delete
}

#### end quickfix

## edit 1.2.2 adding new exception check
@files = <*gb>;

foreach (@files) {
    system $PATH_COPRA_SUBSCRIPTS . "check_for_gene_CDS_features.pl $_ >> gene_CDS_exception.txt";
}

open(MYDATA, "gene_CDS_exception.txt") or die("\nError: cannot open file gene_CDS_exception.txt at homology_intaRNA.pl\n\n");
    my @exception_lines = <MYDATA>;
close MYDATA;


if (scalar(@exception_lines) >= 1) {
    my $exceptionRefSeqs = "";
    foreach(@exception_lines) {
        my @split = split(/\s+/,$_);
        $exceptionRefSeqs = $exceptionRefSeqs . $split[-1] . " ";
    }
    print ERRORLOG "Error: gene but no CDS features present in $exceptionRefSeqs.\n This is most likely connected to currently corrupted RefSeq record(s) at the NCBI.\nPlease resubmit your job without the currently errorous organism(s) or wait some time with your resubmission.\nUsually the files are fixed within ~1 week.\n"; ## edit 1.2.2 added \n ## edit 2.0.2
}
## end CDS gene exception check


## get cluster.tab with DomClust
unless (-e "cluster.tab") { # only do if cluster.tab has not been imported ## edit 2.0.4 changed this to -e

    ### get AA fasta for homolog clustering

    @files = <*gb>;

    foreach(@files) {
        system $PATH_COPRA_SUBSCRIPTS . "get_CDS_from_gbk.pl $_ >> all.fas"; ## edit 2.0.5.1 // removed unless 
    }

    # prep for DomClust
    system "formatdb -i all.fas" unless (-e "all.fas.blast"); ## edit 2.0.1
    # blast sequences
    my $blastallErrorFile = "blastall.error";
    system "blastall -a $cores -p blastp -d all.fas -e 0.001 -i all.fas -Y 1e9 -v 30000 -b 30000 -m 8 -o all.fas.blast 2> $blastallErrorFile" unless (-e "all.fas.blast"); # change the -a parameter to qdjust core usage ## edit 2.0.1 // ## edit 2.0.5.1 // added 2> /dev/null to prevent output to the terminal
    # remove empty error file
    system("rm -f $blastallErrorFile") if ( -z $blastallErrorFile );
    system $PATH_COPRA_SUBSCRIPTS . "blast2homfile.pl all.fas.blast > all.fas.hom"; ## edit 2.0.5.1 // removed -distconv this is now fixed within the script
    system $PATH_COPRA_SUBSCRIPTS . "fasta2genefile.pl all.fas";
    # DomClust
    my $domclustErrorFile = "domclust.error";
    my $domclustExitStatus = system "domclust all.fas.hom all.fas.gene -HO -S -c60 -p0.5 -V0.6 -C80 -o5 > cluster.tab 2> $domclustErrorFile"; ## edit 2.0.5.1 // changed to conda domclust
    $domclustExitStatus /= 256; # get original exit value
    # ensure domclust went fine
    if ($domclustExitStatus != 0) {
    	# restart domclust with --nobreak option
	    my $domclustExitStatus = system "domclust all.fas.hom all.fas.gene -HO -S -c60 -p0.5 -V0.6 -C80 -o5 --nobreak > cluster.tab 2> $domclustErrorFile"; ## edit 2.0.5.1 // changed to conda domclust
	    $domclustExitStatus /= 256; # get original exit value
	    # check if second run was successful
	    if ($domclustExitStatus != 0) {
	    	die("\nERROR: 'domclust' returned with non-zero exit status $domclustExitStatus.\n\n");
	    }
    }
    # remove empty error file
    system("rm -f $domclustErrorFile") if ( -z $domclustErrorFile );

    # edit 2.0.2
    system "grep '>' all.fas | uniq -d > N_chars_in_CDS.txt";
    if (-s "N_chars_in_CDS.txt") {
        print ERRORLOG "'N' characters found in nucleotide CDS. Please remove organism(s) with locus tags:\n";
        system "cat err.log N_chars_in_CDS.txt >> err.log";
    }

}

# 16s sequence parsing 
print $PATH_COPRA_SUBSCRIPTS . "parse_16s_from_gbk.pl $GenBankFiles > 16s_sequences.fa\n" if ($verbose);
system $PATH_COPRA_SUBSCRIPTS . "parse_16s_from_gbk.pl $GenBankFiles > 16s_sequences.fa" unless (-e "16s_sequences.fa");

# check 16s
open(MYDATA, "16s_sequences.fa") or die("\nError: cannot open file 16s_sequences.fa in homology_intaRNA.pl\n\n");
    my @sixteenSseqs = <MYDATA>;
close MYDATA;

my $sixteenScounter = 0;
my $temp_16s_ID = ""; ## edit 2.0.2
foreach (@sixteenSseqs) {
    if ($_ =~ m/>/) {
        $temp_16s_ID = $_; ## edit 2.0.2
        chomp $temp_16s_ID; ## edit 2.0.2
        $sixteenScounter++;
    } else {
        if ($_ =~ m/N/) { print ERRORLOG "\nError: 'N' characters present in 16s_sequences.fa. Remove $temp_16s_ID from the input for the job to execute correctly.\n"; } ## edit 2.0.2
    }
}

if ($sixteenScounter ne $orgcount) {
    my $no16sOrgs = `(grep ">" 16s_sequences.fa && grep ">" input_sRNA.fa) | sort | uniq -u | tr '\n' ' '`; ## edit 2.0.3
    chomp $no16sOrgs; ## edit 2.0.3
    print ERRORLOG "\nError: wrong number of sequences in 16s_sequences.fa.\nOne (or more) of your entered organisms does not contain a correctly annotated 16s rRNA sequence and needs to be removed.\nPlease remove $no16sOrgs\n"; ## edit 2.0.2 
}

## prepare single organism whole genome target predictions 
system "echo $GenBankFiles > merged_refseq_ids.txt"; ## edit 2.0.2 # need this for iterative region plot construction

print $PATH_COPRA_SUBSCRIPTS . "prepare_intarna_out.pl $ncrnas $upfromstartpos $down $mrnapart $GenBankFiles\n" if ($verbose);
system $PATH_COPRA_SUBSCRIPTS . "prepare_intarna_out.pl $ncrnas $upfromstartpos $down $mrnapart $GenBankFiles";
## end  edit 2.0.0

# re-cluster based on 5'UTRs
system "R --slave -f " . $PATH_COPRA_SUBSCRIPTS . "refine_clustertab.r"; #edit jens


# do CopraRNA combination 
## edit 2.0.4 // removed all N*final.csv files as input to combine_clusters.pl
print $PATH_COPRA_SUBSCRIPTS . "combine_clusters.pl $orgcount\n" if ($verbose);
system $PATH_COPRA_SUBSCRIPTS . "combine_clusters.pl $orgcount";

# make annotations
system $PATH_COPRA_SUBSCRIPTS . "annotate_raw_output.pl CopraRNA1_with_pvsample_sorted.csv opt_tags.clustered_rcsize $GenBankFiles > CopraRNA1_anno.csv" if ($cop1); ## edit 2.0.6
system $PATH_COPRA_SUBSCRIPTS . "annotate_raw_output.pl CopraRNA2_prep_sorted.csv opt_tags.clustered $GenBankFiles > CopraRNA2_prep_anno.csv"; ## edit 2.0.6

# get additional homologs in cluster.tab
system $PATH_COPRA_SUBSCRIPTS . "parse_homologs_from_domclust_table.pl CopraRNA1_anno.csv cluster.tab > CopraRNA1_anno_addhomologs.csv" if ($cop1); ## edit 2.0.6
system $PATH_COPRA_SUBSCRIPTS . "parse_homologs_from_domclust_table.pl CopraRNA2_prep_anno.csv cluster.tab > CopraRNA2_prep_anno_addhomologs.csv"; ## edit 2.0.6

# add corrected p-values (padj) - first column
system "awk -F',' '{ print \$1 }' CopraRNA1_anno_addhomologs.csv > CopraRNA1_pvalues.txt" if ($cop1); ## edit 2.0.6
# just for formatting
system "awk -F',' '{ print \$1 }' CopraRNA2_prep_anno_addhomologs.csv > CopraRNA2_pvalues.txt"; ## edit 2.0.6

system "R --slave -f $PATH_COPRA_SUBSCRIPTS/calc_padj.R --args CopraRNA1_pvalues.txt" if ($cop1); ## edit 2.0.6
system "paste padj.csv CopraRNA1_anno_addhomologs.csv -d ',' > CopraRNA1_anno_addhomologs_padj.csv" if ($cop1); ## edit 2.0.6

# just for formatting
system "R --slave -f $PATH_COPRA_SUBSCRIPTS/calc_padj.R --args CopraRNA2_pvalues.txt";
system "paste padj.csv CopraRNA2_prep_anno_addhomologs.csv -d ',' > CopraRNA2_prep_anno_addhomologs_padj.csv"; ## edit 2.0.6

# add amount sampled values CopraRNA 1 // CopraRNA 2 has no sampling
system $PATH_COPRA_SUBSCRIPTS . "get_amount_sampled_values_and_add_to_table.pl CopraRNA1_anno_addhomologs_padj.csv 0 > CopraRNA1_anno_addhomologs_padj_amountsamp.csv" if ($cop1); ## edit 2.0.6
# make consistent names
system "mv CopraRNA1_anno_addhomologs_padj_amountsamp.csv CopraRNA1_final_all.csv" if ($cop1); ## edit 2.0.6
system $PATH_COPRA_SUBSCRIPTS . "get_amount_sampled_values_and_add_to_table.pl CopraRNA2_prep_anno_addhomologs_padj.csv 1 > CopraRNA2_prep_anno_addhomologs_padj_amountsamp.csv"; ## edit 2.0.6

# get ooi refseq id
my @split = split(/\s/, $refseqaffiliations{$ARGV[4]});
# get the first id entry
my $ooi_refseq_id = $split[0];


### edit jens
unless ($cop1) {
    # align homologous targets
    system $PATH_COPRA_SUBSCRIPTS . "parallelize_target_alignments.pl CopraRNA2_prep_anno_addhomologs_padj_amountsamp.csv";
    # run position script
    system "cp " . $PATH_COPRA_SUBSCRIPTS . "CopraRNA_available_organisms.txt ."; ## edit 2.0.6
    #system "R --slave -f " . $PATH_COPRA_SUBSCRIPTS . "copraRNA2_position_script_for_evo_precalculated_alignments_w_ooi.R --args $ooi_refseq_id 2> /dev/null > /dev/null"; ## edit 2.0.6
	system "R --slave -f " . $PATH_COPRA_SUBSCRIPTS . "copraRNA2_phylogenetic_sorting.r > /dev/null > /dev/null"; #edit jens
    # perform actual CopraRNA 2 p-value combination
   # system "R --slave -f " . $PATH_COPRA_SUBSCRIPTS . "join_pvals_coprarna2.R --args $ooi_refseq_id ooi_consensus overall_consensus 2> /dev/null > /dev/null"; ## edit 2.0.6
	system "R --slave -f " . $PATH_COPRA_SUBSCRIPTS . "join_pvals_coprarna_2.r > /dev/null > /dev/null"; #edit jens
    
}

#edit jens
# truncate final output // ## edit 2.0.5.1
system "head -n $topcount CopraRNA1_final_all.csv > CopraRNA1_final.csv" if ($cop1); ## edit 2.0.6
unless ($cop1) {
    system "head -n $topcount CopraRNA_result_all.csv > CopraRNA_result.csv"; ## edit 2.0.6
    #system "head -n $topcount CopraRNA2_final_all_balanced.csv > CopraRNA2_final_balanced.csv"; ## edit 2.0.6
    #system "head -n $topcount CopraRNA2_final_all_balanced_consensus.csv > CopraRNA2_final_balanced_consensus.csv"; ## edit 2.0.6
    #system "head -n $topcount CopraRNA2_final_all_ooi_consensus.csv > CopraRNA2_final_ooi_consensus.csv"; ## edit 2.0.6
    #system "head -n $topcount CopraRNA2_final_all_ooi_ooiconsensus.csv > CopraRNA2_final_ooi_ooiconsensus.csv"; ## edit 2.0.6
}

#edit jens
# # figure out which result is the primary result ## edit 2.0.6
# if ($cop1) { # CopraRNA 1 is the primary requested result
    # system "cp CopraRNA1_final.csv CopraRNA_result.csv";    
    # system "cp CopraRNA1_final_all.csv CopraRNA_result_all.csv";    
# } elsif ($nooi and (not $cons)) { # CopraRNA 2 with balanced mode is the requested result
    # system "cp  CopraRNA2_final_balanced.csv CopraRNA_result.csv";
    # system "cp  CopraRNA2_final_all_balanced.csv CopraRNA_result_all.csv";
# } elsif ($nooi and ($cons eq 2)) { # CopraRNA 2 balanced prediction with overall consensus
    # system "cp CopraRNA2_final_balanced_consensus.csv CopraRNA_result.csv"; 
    # system "cp CopraRNA2_final_all_balanced_consensus.csv CopraRNA_result_all.csv";
# } elsif ($cons eq 1) { # CopraRNA 2 ooi prediction with ooi consensus
    # system "cp CopraRNA2_final_ooi_ooiconsensus.csv CopraRNA_result.csv";
    # system "cp CopraRNA2_final_all_ooi_ooiconsensus.csv CopraRNA_result_all.csv"; 
# } elsif ($cons eq 2) { # CopraRNA 2 ooi prediction with overall consensus
    # system "cp CopraRNA2_final_ooi_consensus.csv CopraRNA_result.csv";
    # system "cp CopraRNA2_final_all_ooi_consensus.csv CopraRNA_result_all.csv"; 
# } else { # CopraRNA 2 with org of interest focus (standard)
    # system "cp CopraRNA2_final_ooi.csv CopraRNA_result.csv";
    # system "cp CopraRNA2_final_all_ooi.csv CopraRNA_result_all.csv";
# }

# filtering for ooi single p-value
if ($ooi_filt) {

    my @not_filtered_list = (); # values below the p-value threshold
    my @filtered_list = ();     # values empty or above the p-value threshold

    open(MYDATA, "CopraRNA_result_all.csv") or die("\nError: cannot open file CopraRNA_result_all.csv at homology_intaRNA.pl\n\n");
        my @CopraRNA_all_out_lines = <MYDATA>;
    close MYDATA;

    push(@not_filtered_list, $CopraRNA_all_out_lines[0]); # header

    for (my $i=1;$i<scalar(@CopraRNA_all_out_lines);$i++) {
        my $curr_line = $CopraRNA_all_out_lines[$i];
        my @split = split(/,/,$curr_line);
        my $curr_ooi_cell = $split[2];
        if ($curr_ooi_cell) {
            my @split_ooi_cell = split(/\|/,$curr_ooi_cell);
            my $curr_ooi_pv = $split_ooi_cell[2];
            if($curr_ooi_pv<=$ooi_filt) { # smaller or eq to the set ooi_filt threshold
                push(@not_filtered_list, $curr_line);
            } else { # bigger tahn the set ooi_filt threshold
                push(@filtered_list, $curr_line);
            }
        } else { # empty cell
            push(@filtered_list, $curr_line);
        }
    }
    # print
    open WRITEFILT, ">", "CopraRNA_result_all_filt.csv";

    foreach(@not_filtered_list) {
        print WRITEFILT $_;
    }
    foreach(@filtered_list) {
        print WRITEFILT $_;
    }
    close WRITEFILT;
    system "cp CopraRNA_result_all_filt.csv CopraRNA_result_all.csv";
    system "head -n $topcount CopraRNA_result_all.csv > CopraRNA_result.csv";
}

# plot CopraRNA 2 evo heatmap and jalview files for selection
unless ($cop1) {


    system "R --slave -f " . $PATH_COPRA_SUBSCRIPTS . "copraRNA2_find_conserved_sites.r > /dev/null > /dev/null"; ## edit jens
	system "R --slave -f " . $PATH_COPRA_SUBSCRIPTS . "copraRNA2_conservation_heatmaps.r > /dev/null > /dev/null"; ## edit jens
    system "rm CopraRNA_available_organisms.txt"; ## edit 2.0.6
}

# check for run fail CopraRNA
open(MYDATA, "CopraRNA_result.csv") or die("\nError: cannot open file CopraRNA_result.csv at homology_intaRNA.pl\n\n");
    my @CopraRNA_out_lines = <MYDATA>;
close MYDATA;

if (scalar(@CopraRNA_out_lines) <= 1) { ## edit 2.0.6
    print ERRORLOG "Error: No predictions in CopraRNA_result.csv. CopraRNA run failed.\n"; ## edit 2.0.2
}

# trim off last column (initial_sorting) if CopraRNA 2 prediction mode
unless ($cop1) {
     system "awk -F',' '{ print \$NF }' CopraRNA_result.csv > CopraRNA_result.map_evo_align" if ($websrv);
     system "awk -F, -vOFS=, '{NF-=1;print}' CopraRNA_result.csv > CopraRNA_result_temp.csv";
     system "mv CopraRNA_result_temp.csv CopraRNA_result.csv";
     system "awk -F, -vOFS=, '{NF-=1;print}' CopraRNA_result_all.csv > CopraRNA_result_all_temp.csv";
     system "mv CopraRNA_result_all_temp.csv CopraRNA_result_all.csv";
     # change header
     system "sed -i 's/,Additional.homologs,/,Additional homologs,/g' CopraRNA_result.csv";
     system "sed -i 's/,Amount.sampled/,Amount sampled/g' CopraRNA_result.csv";
     system "sed -i 's/p.value/p-value/g' CopraRNA_result.csv";
     system "sed -i 's/,Additional.homologs,/,Additional homologs,/g' CopraRNA_result_all.csv";
     system "sed -i 's/,Amount.sampled/,Amount sampled/g' CopraRNA_result_all.csv";
     system "sed -i 's/p.value/p-value/g' CopraRNA_result_all.csv";
}

if ($websrv) { # only if webserver output is requested via -websrv ## edit 2.0.5.1

    my $allrefs = $refseqaffiliations{$ARGV[4]};
    my @splitallrefs = split(/\s/,$allrefs);

    my $themainrefid = $splitallrefs[0]; # organism of interest RefSeq ID
    my $orgofintTargets = $themainrefid . "_upfromstartpos_" . $upfromstartpos . "_down_" . $down . ".fa";
    my $orgofintsRNA = "ncRNA_" . $themainrefid . ".fa";

    # returns comma separated locus tags (first is always refseq ID). Example: NC_000913,b0681,b1737,b1048,b4175,b0526,b1093,b1951,,b3831,b3133,b0886,,b3176 
    my $top_predictons_locus_tags = `awk -F',' '{print \$3}' CopraRNA_result.csv | sed 's/(.*)//g' | tr '\n' ','`; ## edit 2.0.6 switched to generic output file

    # split
    my @split = split(/,/, $top_predictons_locus_tags);
    
    # remove RefSeqID
    shift @split;

    foreach (@split) {
        if ($_) {
            system "grep -iA1 '$_' $orgofintTargets >> CopraRNA_top_targets.fa";
        }
    }

    system "IntaRNA_1ui.pl -t CopraRNA_top_targets.fa -m $orgofintsRNA -o -w $winsize -L $maxbpdist > Cop_IntaRNA1_ui.intarna";
    # fix for ambiguous nt in intarna output
    system "sed -i '/contains ambiguous IUPAC nucleotide encodings/d' Cop_IntaRNA1_ui.intarna";

    system $PATH_COPRA_SUBSCRIPTS . "prepare_output_for_websrv_new.pl CopraRNA_result.csv Cop_IntaRNA1_ui.intarna";
    system "mv coprarna_internal_table.csv coprarna_websrv_table.csv";

    system "cp $orgofintTargets target_sequences_orgofint.fa";
}

system $PATH_COPRA_SUBSCRIPTS . "print_archive_README.pl > README.txt";

if ($enrich) { ## edit 2.0.5.1 // ## edit 2.0.6 changes in file names

    ##### create DAVID enrichment table
    ## this has all been changed to python in version 2.0.3.1 because the DAVID-WS perl client was flawed
    system $PATH_COPRA_SUBSCRIPTS . "DAVIDWebService_CopraRNA.py CopraRNA_result_all.csv $enrich > DAVID_enrichment_temp.txt"; ## edit 2.0.5.1 // added $enrich as input 
    system "grep -P 'termName\\s=|categoryName\\s=|score\\s=|listHits\\s=|percent\\s=|ease\\s=|geneIds\\s=|listTotals\\s=|popHits\\s=|popTotals\\s=|foldEnrichment\\s=|bonferroni\\s=|benjamini\\s=|afdr\\s=' DAVID_enrichment_temp.txt | sed 's/^[ ]*//g' | sed 's/ = /=/g' | sed 's/, /,/g' > DAVID_enrichment_grepped_temp.txt"; ## edit 2.0.6 // only removing obsolete spaces and keeping others
    system $PATH_COPRA_SUBSCRIPTS . "make_enrichment_table_from_py_output.pl DAVID_enrichment_grepped_temp.txt > termClusterReport.txt"; ## edit 2.0.3.1

    open(MYDATA, "termClusterReport.txt") or system "echo 'If you are reading this, then your prediction did not return an enrichment, your organism of interest is not in the DAVID database\nor the DAVID webservice is/was termporarily down. You can either rerun your CopraRNA\nprediction or create your enrichment manually at the DAVID homepage.' > termClusterReport.txt";
        my @enrichment_lines = <MYDATA>;
    close MYDATA;

    unless($enrichment_lines[0]) {
        system "echo -e 'If you are reading this, then your prediction did not return an enrichment, your organism of interest is not in the DAVID database\nor the DAVID webservice is/was termporarily down. You can either rerun your CopraRNA\nprediction or create your enrichment manually at the DAVID homepage.' > termClusterReport.txt";
    }

    ##### end DAVID enrichment

    ## enrichment visualization ## edit 1.2.5
    system "cp $PATH_COPRA_SUBSCRIPTS" . "copra_heatmap.html ."; ## edit 1.2.5 ## edit 1.2.7 (edited html file)
    system "R --slave -f " . $PATH_COPRA_SUBSCRIPTS . "extract_functional_enriched.R --args CopraRNA_result_all.csv termClusterReport.txt enrichment.txt"; ## edit 1.2.5 ## edit 1.2.7 (edited R code) ## edit 2.0.5.1 added args  # edit jens R - script
    system $PATH_COPRA_SUBSCRIPTS . "make_heatmap_json.pl enrichment.txt"; ## edit 1.2.5
    system "cp $PATH_COPRA_SUBSCRIPTS" . "index-thumb.html ."; ## edit 1.2.5
    system "cp $PATH_COPRA_SUBSCRIPTS" . "index-pdf.html ."; ## edit 1.2.6
    system "phantomjs " . $PATH_COPRA_SUBSCRIPTS . "rasterize.js " . "./index-thumb.html enriched_heatmap_big.png"; ## edit 2.0.6 // phantomjs now via conda and in path
    system "phantomjs " . $PATH_COPRA_SUBSCRIPTS . "rasterize.js " . "./index-pdf.html enriched_heatmap_big.pdf"; ## edit 2.0.6
    system "rm index-thumb.html"; ## edit 1.2.5
    system "rm index-pdf.html"; ## edit 1.2.6
    ## end add enrichment vis
}


close ERRORLOG;

